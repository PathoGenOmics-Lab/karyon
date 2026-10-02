//! Just enough XML to read back what [`SvgWriter`](crate::SvgWriter) wrote.
//!
//! The writer emits one shape of markup: double-quoted attributes, the five
//! named entities, no namespaces beyond the root's `xmlns`, no comments and no
//! doctype. Reading exactly that would be fifty lines. This reads a little
//! more, single quotes, numeric references, comments, processing
//! instructions, a doctype and CDATA, because [`Pdf::from_svg`] is public and
//! an SVG that came out of the writer and through an editor that added a
//! comment should still convert. Nothing is validated: a tag that never closes
//! is closed by the reader, as [`SvgWriter::finish`](crate::SvgWriter::finish)
//! closes the groups a track left open, and an end tag with no start is let go.
//!
//! [`Pdf::from_svg`]: crate::pdf::Pdf::from_svg

use std::borrow::Cow;

/// One piece of the document.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Token<'a> {
    /// `<name ...>`, or `<name .../>` when `empty` is true.
    Start {
        name: &'a str,
        attributes: Attributes<'a>,
        empty: bool,
    },
    /// `</name>`.
    End(&'a str),
    /// Character data, still escaped.
    Text(&'a str),
    /// The inside of a CDATA section, which is not escaped.
    Raw(&'a str),
}

/// The attributes of a start tag, read when they are asked for.
///
/// Kept as the slice of the tag they came in, so the hundreds of thousands of
/// rectangles in a large figure cost no allocation to walk past.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Attributes<'a>(&'a str);

impl<'a> Attributes<'a> {
    /// The raw value of `name`, still escaped, if the tag has one.
    pub(crate) fn get(&self, name: &str) -> Option<&'a str> {
        self.iter()
            .find(|(given, _)| *given == name)
            .map(|(_, value)| value)
    }

    /// Every attribute, in the order it was written, as name and raw value.
    pub(crate) fn iter(&self) -> AttributeIter<'a> {
        AttributeIter { rest: self.0 }
    }
}

/// Walks the attributes of one tag.
pub(crate) struct AttributeIter<'a> {
    rest: &'a str,
}

impl<'a> Iterator for AttributeIter<'a> {
    type Item = (&'a str, &'a str);

    fn next(&mut self) -> Option<Self::Item> {
        let rest = self.rest.trim_start();
        let name_end = rest.find(|c: char| c == '=' || c.is_whitespace())?;
        let name = &rest[..name_end];
        let after = rest[name_end..]
            .trim_start()
            .strip_prefix('=')?
            .trim_start();
        let quote = after.chars().next().filter(|c| *c == '"' || *c == '\'')?;
        let body = &after[1..];
        let close = body.find(quote)?;
        self.rest = &body[close + 1..];
        Some((name, &body[..close]))
    }
}

/// The tokens of a document, in order.
pub(crate) struct Tokens<'a> {
    rest: &'a str,
}

/// Reads `document` a token at a time.
pub(crate) fn tokens(document: &str) -> Tokens<'_> {
    Tokens { rest: document }
}

impl<'a> Iterator for Tokens<'a> {
    type Item = Token<'a>;

    fn next(&mut self) -> Option<Token<'a>> {
        loop {
            if self.rest.is_empty() {
                return None;
            }
            if !self.rest.starts_with('<') {
                let end = self.rest.find('<').unwrap_or(self.rest.len());
                let text = &self.rest[..end];
                self.rest = &self.rest[end..];
                return Some(Token::Text(text));
            }
            if let Some(after) = self.rest.strip_prefix("<!--") {
                self.rest = skip_past(after, "-->");
                continue;
            }
            if let Some(after) = self.rest.strip_prefix("<?") {
                self.rest = skip_past(after, "?>");
                continue;
            }
            if let Some(after) = self.rest.strip_prefix("<![CDATA[") {
                let end = after.find("]]>").unwrap_or(after.len());
                self.rest = skip_past(after, "]]>");
                return Some(Token::Raw(&after[..end]));
            }
            if let Some(after) = self.rest.strip_prefix("<!") {
                // A doctype, whose internal subset can hold a `>` of its own
                // inside brackets.
                let close = after.find('>').unwrap_or(after.len());
                self.rest = match after.find('[') {
                    Some(open) if open < close => skip_past(after, "]>"),
                    _ => skip_past(after, ">"),
                };
                continue;
            }
            if let Some(after) = self.rest.strip_prefix("</") {
                let end = after.find('>').unwrap_or(after.len());
                let name = after[..end].trim();
                self.rest = skip_past(after, ">");
                return Some(Token::End(name));
            }
            let after = &self.rest[1..];
            let close = tag_end(after);
            let inside = &after[..close];
            self.rest = after.get(close + 1..).unwrap_or("");
            let (inside, empty) = match inside.trim_end().strip_suffix('/') {
                Some(open) => (open, true),
                None => (inside, false),
            };
            let name_end = inside
                .find(|c: char| c.is_whitespace() || c == '/')
                .unwrap_or(inside.len());
            let name = &inside[..name_end];
            if name.is_empty() {
                // A `<` that opens nothing, as in text written unescaped.
                continue;
            }
            return Some(Token::Start {
                name,
                attributes: Attributes(&inside[name_end..]),
                empty,
            });
        }
    }
}

/// Where the `>` closing a start tag is, past any quoted `>` inside it.
fn tag_end(tag: &str) -> usize {
    let mut quote: Option<u8> = None;
    for (at, byte) in tag.bytes().enumerate() {
        match (quote, byte) {
            (None, b'"' | b'\'') => quote = Some(byte),
            (Some(open), _) if open == byte => quote = None,
            (None, b'>') => return at,
            _ => {}
        }
    }
    tag.len()
}

/// Everything after the first `marker` in `text`, or nothing if there is none.
fn skip_past<'a>(text: &'a str, marker: &str) -> &'a str {
    match text.find(marker) {
        Some(at) => &text[at + marker.len()..],
        None => "",
    }
}

/// Text with its entities replaced by the characters they stand for.
///
/// Borrowed when there is nothing to replace, which is nearly always: the
/// writer escapes only the five metacharacters, and only in names and labels.
/// A reference to something that is no character is left as it was written.
pub(crate) fn unescape(text: &str) -> Cow<'_, str> {
    if !text.contains('&') {
        return Cow::Borrowed(text);
    }
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.find('&') {
        out.push_str(&rest[..at]);
        let after = &rest[at + 1..];
        let replaced = after.find(';').and_then(|end| {
            let name = &after[..end];
            let c = match name {
                "amp" => Some('&'),
                "lt" => Some('<'),
                "gt" => Some('>'),
                "quot" => Some('"'),
                "apos" => Some('\''),
                _ => {
                    let number = name.strip_prefix('#')?;
                    let hex = number
                        .strip_prefix('x')
                        .or_else(|| number.strip_prefix('X'));
                    let code = match hex {
                        Some(hex) => u32::from_str_radix(hex, 16).ok()?,
                        None => number.parse::<u32>().ok()?,
                    };
                    char::from_u32(code)
                }
            }?;
            Some((c, end))
        });
        match replaced {
            Some((c, end)) => {
                out.push(c);
                rest = &after[end + 1..];
            }
            None => {
                out.push('&');
                rest = after;
            }
        }
    }
    out.push_str(rest);
    Cow::Owned(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(document: &str) -> Vec<String> {
        tokens(document)
            .map(|token| match token {
                Token::Start { name, empty, .. } => {
                    format!("<{name}{}>", if empty { "/" } else { "" })
                }
                Token::End(name) => format!("</{name}>"),
                Token::Text(text) => format!("text {text}"),
                Token::Raw(text) => format!("raw {text}"),
            })
            .collect()
    }

    #[test]
    fn the_writer_s_markup_comes_back_as_its_tags() {
        let read = names(
            r#"<svg width="10"><g clip-path="url(#a)"><rect x="1"/></g><text x="0">a &amp; b</text></svg>"#,
        );
        assert_eq!(
            read,
            [
                "<svg>",
                "<g>",
                "<rect/>",
                "</g>",
                "<text>",
                "text a &amp; b",
                "</text>",
                "</svg>"
            ]
        );
    }

    #[test]
    fn comments_declarations_and_processing_instructions_are_passed_over() {
        let read = names(
            "<?xml version=\"1.0\"?><!DOCTYPE svg [ <!ENTITY a \"b\"> ]><!-- a > b --><svg/><![CDATA[x < y]]>",
        );
        assert_eq!(read, ["<svg/>", "raw x < y"]);
    }

    #[test]
    fn a_quoted_closing_bracket_does_not_end_the_tag() {
        let mut read = tokens(r#"<text a="1 > 0" b='it said "x"'/>"#);
        let Some(Token::Start {
            name,
            attributes,
            empty,
        }) = read.next()
        else {
            panic!("no start tag");
        };
        assert_eq!(name, "text");
        assert!(empty);
        assert_eq!(attributes.get("a"), Some("1 > 0"));
        assert_eq!(attributes.get("b"), Some(r#"it said "x""#));
        assert_eq!(attributes.get("c"), None);
        assert!(read.next().is_none());
    }

    #[test]
    fn entities_are_replaced_and_unknown_ones_are_kept() {
        assert_eq!(unescape("a &lt;b&gt; &amp; &quot;c&apos;"), "a <b> & \"c'");
        assert_eq!(unescape("&#969; &#x2264;"), "\u{3c9} \u{2264}");
        assert_eq!(unescape("&nbsp; & &#xZZ;"), "&nbsp; & &#xZZ;");
        assert!(matches!(unescape("plain"), Cow::Borrowed("plain")));
    }
}
