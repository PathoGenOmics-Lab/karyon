//! The page `karyon man` prints: `karyon help all` set in roff, for man.
//!
//! The page is the help and not a second text kept beside it. Each heading
//! of the help is a section of the page, each entry under one a tagged
//! paragraph, and the prose and the commands to type are set as they are
//! written. Only what a manual page has and the help has not is written
//! here: the name line, a synopsis made of the usage and the commands, the
//! exit status and where the guide is. So a flag written up in the help is
//! in the page with no edit here, and a test reads the whole help back out
//! of it.

use karyon::cli::args;

use super::{commands, split_entry, GUIDE, HELP};

/// The day the version the page describes was released, which a test holds
/// to the changelog's heading for that version.
const DATE: &str = "2026-07-31";

/// The page, in roff for the `man` macros.
pub(crate) fn page() -> String {
    let mut page = Page {
        out: format!(
            ".TH KARYON 1 {DATE} \"karyon {}\" \"User Commands\"\n",
            karyon::VERSION
        ),
        fresh: false,
        flags: args::FLAGS.iter().copied().chain(["-h", "-V"]).collect(),
    };
    let mut blocks = blocks(HELP).into_iter();
    // `karyon, genomic track plots on one shared coordinate axis`.
    let title = blocks
        .next()
        .and_then(|block| block.first().copied())
        .unwrap_or_default();
    let what = title.split_once(", ").map_or(title, |(_, what)| what);
    page.heading("NAME");
    page.line(format!("karyon \\- {}", escape(what)));

    let mut described = false;
    for block in blocks {
        let (first, rest) = block.split_first().expect("a block holds a line");
        if heading(first) {
            if *first == "USAGE" {
                page.synopsis(rest);
            } else {
                page.heading(first);
                page.body(rest);
            }
        } else if !first.starts_with(' ') {
            // The help's paragraph at the margin, under its usage, which is
            // what the program does.
            if !described {
                page.heading("DESCRIPTION");
                described = true;
            }
            page.prose(&block);
        } else {
            page.body(&block);
        }
    }

    page.heading("EXIT STATUS");
    page.prose(&[
        "A failing command prints one line to standard error, starting with",
        "karyon:, and exits with status 1. Success exits with 0, and so do the",
        "help, a completion script and this page.",
    ]);
    page.heading("SEE ALSO");
    let guide = format!("The guide, with every track drawn and its options at work: {GUIDE}");
    page.prose(&[guide.as_str()]);
    page.out
}

/// The page as it is written, and what it needs to know to go on.
struct Page {
    out: String,
    /// Whether a section heading was the last thing written, after which a
    /// paragraph needs no macro of its own: `.PP` there is one mandoc warns
    /// it skips.
    fresh: bool,
    /// Every spelling of an option, set in bold where the prose names one.
    flags: Vec<&'static str>,
}

impl Page {
    fn heading(&mut self, text: &str) {
        self.out.push_str(&format!(".SH {}\n", escape(text)));
        self.fresh = true;
    }

    fn paragraph(&mut self) {
        if !self.fresh {
            self.out.push_str(".PP\n");
        }
        self.fresh = false;
    }

    /// A line of text, guarded against reading as a request: roff takes a
    /// line that starts with a full stop or an apostrophe for one.
    fn line(&mut self, text: String) {
        if text.starts_with(['.', '\'']) {
            self.out.push_str("\\&");
        }
        self.out.push_str(&text);
        self.out.push('\n');
        self.fresh = false;
    }

    /// The usage line, and each command after it in the form it is typed.
    fn synopsis(&mut self, usage: &[&str]) {
        self.heading("SYNOPSIS");
        let mut forms: Vec<(String, String)> = usage
            .iter()
            .map(|line| {
                let line = line.trim();
                let (name, rest) = line.split_once(' ').unwrap_or((line, ""));
                (name.to_string(), escape(rest))
            })
            .collect();
        for (_, entry) in commands() {
            let (names, value, _) = split_entry(&entry);
            let rest = value.map_or(String::new(), |value| {
                format!("\\fI<{}>\\fR", escape(value))
            });
            forms.push((names.join(" "), rest));
        }
        for (at, (name, rest)) in forms.iter().enumerate() {
            if at > 0 {
                self.out.push_str(".br\n");
            }
            let space = if rest.is_empty() { "" } else { " " };
            self.line(format!("\\fB{}\\fR{space}{rest}", escape(name)));
        }
    }

    /// What follows a heading: commands to type, entries, or prose.
    fn body(&mut self, lines: &[&str]) {
        let Some(first) = lines.first() else {
            return;
        };
        if first.trim_end().ends_with('\\') {
            // A command to type, kept as written, with its indentation.
            self.paragraph();
            self.out.push_str(".RS 4\n.nf\n");
            for &line in lines {
                self.line(escape(line.strip_prefix("    ").unwrap_or(line)));
            }
            self.out.push_str(".fi\n.RE\n");
        } else if first.starts_with("    -") || first.starts_with("    karyon ") {
            // An entry starts at the help's indent, and its lines go on
            // deeper, from where its text starts.
            let mut entry: Vec<&str> = Vec::new();
            for line in lines {
                if !line.starts_with("     ") && !entry.is_empty() {
                    self.entry(&entry);
                    entry.clear();
                }
                entry.push(line);
            }
            self.entry(&entry);
        } else {
            self.prose(lines);
        }
    }

    /// One entry as a tagged paragraph: what is typed, in bold, and the
    /// value it takes, in italics, over what it does.
    fn entry(&mut self, lines: &[&str]) {
        let Some((first, more)) = lines.split_first() else {
            return;
        };
        let (names, value, said) = split_entry(first);
        // A command is typed as one phrase, `karyon man`; the spellings of
        // an option are each a word of their own, `-o, --output`.
        let mut tag = if names.first() == Some(&"karyon") {
            format!("\\fB{}\\fR", escape(&names.join(" ")))
        } else {
            names
                .iter()
                .map(|name| format!("\\fB{}\\fR", escape(name)))
                .collect::<Vec<_>>()
                .join(", ")
        };
        if let Some(value) = value {
            tag.push_str(&format!(" \\fI<{}>\\fR", escape(value)));
        }
        self.out.push_str(".TP\n");
        self.line(tag);
        let mut text = said;
        for line in more {
            text.push(' ');
            text.push_str(line.trim());
        }
        self.fill(&text);
    }

    /// A paragraph of the help's prose.
    fn prose(&mut self, lines: &[&str]) {
        self.paragraph();
        self.fill(&lines.join(" "));
    }

    /// Text as lines of roff no longer than mandoc asks a page's source to
    /// be, 80 bytes, escapes and all, with room for the guard a line may need.
    fn fill(&mut self, text: &str) {
        let mut line = String::new();
        for word in text.split_whitespace() {
            let set = self.set(word);
            if !line.is_empty() && line.len() + 1 + set.len() > 78 {
                self.line(std::mem::take(&mut line));
            }
            if !line.is_empty() {
                line.push(' ');
            }
            line.push_str(&set);
        }
        if !line.is_empty() {
            self.line(line);
        }
    }

    /// A word of prose, escaped, and in bold where it is a flag, the
    /// punctuation around it left as it was.
    fn set(&self, word: &str) -> String {
        let flag = word
            .trim_start_matches(['(', '\''])
            .trim_end_matches([',', '.', ';', ':', ')', '\'']);
        match word.find(flag) {
            Some(at) if !flag.is_empty() && self.flags.contains(&flag) => format!(
                "{}\\fB{}\\fR{}",
                escape(&word[..at]),
                escape(flag),
                escape(&word[at + flag.len()..])
            ),
            _ => escape(word),
        }
    }
}

/// The help in blocks, as its blank lines part it.
fn blocks(text: &str) -> Vec<Vec<&str>> {
    let mut blocks: Vec<Vec<&str>> = vec![Vec::new()];
    for line in text.lines() {
        match blocks.last_mut() {
            Some(block) if !line.trim().is_empty() => block.push(line),
            Some(block) if !block.is_empty() => blocks.push(Vec::new()),
            _ => {}
        }
    }
    blocks.retain(|block| !block.is_empty());
    blocks
}

/// Whether a line of the help is a heading: at the margin, opening with a
/// word in capitals, as `TRACKS` and `TRACK OPTIONS, each describing...` do
/// and the prose at the margin does not.
fn heading(line: &str) -> bool {
    let word = line.split([' ', ',']).next().unwrap_or_default();
    word.len() > 1 && word.chars().all(|c| c.is_ascii_uppercase())
}

/// `text` as roff sets it. A backslash is written `\e`, a hyphen-minus `\-`,
/// so a flag copied off the page is the flag, and a character beyond ASCII
/// its Unicode escape, which groff and mandoc both read.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\e"),
            '-' => out.push_str("\\-"),
            c if c.is_ascii() => out.push(c),
            c => out.push_str(&format!("\\[u{:04X}]", u32::from(c))),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The page as its reader sees its words: requests dropped but for a
    /// heading's, and fonts and escapes undone.
    fn plain(page: &str) -> String {
        let mut text = String::new();
        for line in page.lines() {
            let line = match line.strip_prefix(".SH ") {
                Some(heading) => heading,
                None if line.starts_with('.') => continue,
                None => line,
            };
            let mut line = ["\\fB", "\\fI", "\\fR", "\\&"]
                .iter()
                .fold(line.to_string(), |line, font| line.replace(font, ""))
                .replace("\\-", "-");
            while let Some(at) = line.find("\\[u") {
                let end = at + line[at..].find(']').expect("an escape ends");
                let code = u32::from_str_radix(&line[at + 3..end], 16).expect("a code");
                let character = char::from_u32(code).expect("a character");
                line.replace_range(at..=end, &character.to_string());
            }
            text.push_str(&line.replace("\\e", "\\"));
            text.push('\n');
        }
        text
    }

    /// Against the parser's list: a flag written up in the help is an entry
    /// of the page, set in bold as its tag.
    #[test]
    fn every_flag_the_parser_answers_to_has_an_entry() {
        let page = page();
        let tags: Vec<&str> = page
            .split(".TP\n")
            .skip(1)
            .filter_map(|entry| entry.lines().next())
            .collect();
        for flag in args::FLAGS.iter().copied().chain(["-h", "-V"]) {
            let bold = format!("\\fB{}\\fR", escape(flag));
            assert!(
                tags.iter().any(|tag| tag
                    .split(", ")
                    .any(|name| name.split(' ').next() == Some(bold.as_str()))),
                "{flag} has no entry in the page"
            );
        }
    }

    /// Every word of `karyon help all`, in its order, so nothing the help
    /// says is dropped on the way to the page. The usage is the page's
    /// synopsis, under the heading man pages give it.
    #[test]
    fn the_page_holds_every_word_of_the_help_in_its_order() {
        let page = plain(&page());
        let bare = |word: &str| word.trim_matches(|c: char| ",.;:".contains(c)).to_string();
        let mut said = page.split_whitespace().map(bare);
        let mut matched = 0;
        for word in HELP.split_whitespace().filter(|word| *word != "USAGE") {
            let word = bare(word);
            assert!(
                said.any(|there| there == word),
                "{word}, the help's word {matched}, is not in the page where the help has it"
            );
            matched += 1;
        }
        assert!(matched > 3000, "{matched} words");
    }

    /// What mandoc's lint asks of a page, held here where CI has no mandoc:
    /// only the requests this writes, no line roff would take for a request,
    /// nothing at the end of a line, no line of filled text over 80 bytes,
    /// and no backslash but the escapes written on purpose.
    #[test]
    fn the_page_is_roff_that_man_reads() {
        let page = page();
        let mut filled = true;
        for line in page.lines() {
            assert_eq!(line, line.trim_end(), "space at the end of a line");
            if let Some(request) = line.strip_prefix('.') {
                let name = request.split(' ').next().unwrap_or_default();
                let known = ["TH", "SH", "TP", "PP", "RS", "RE", "nf", "fi", "br"];
                assert!(known.contains(&name), "{line}");
                filled = match name {
                    "nf" => false,
                    "fi" => true,
                    _ => filled,
                };
                continue;
            }
            assert!(!line.starts_with('\''), "{line}");
            assert!(!filled || line.len() <= 80, "{} bytes: {line}", line.len());
            let mut rest = line;
            while let Some(at) = rest.find('\\') {
                rest = &rest[at + 1..];
                let escapes = ["-", "e", "fB", "fI", "fR", "&", "[u"];
                assert!(
                    escapes.iter().any(|escape| rest.starts_with(escape)),
                    "{line}"
                );
            }
        }
        assert!(page.starts_with(".TH KARYON 1 "), "{page}");
        assert!(page.contains("\n.SH NAME\nkaryon \\- genomic track plots"));
        for heading in [
            "SYNOPSIS",
            "DESCRIPTION",
            "TRACKS",
            "COMMANDS",
            "EXIT STATUS",
        ] {
            assert!(page.contains(&format!("\n.SH {heading}")), "{heading}");
        }
    }

    /// The date on the page is the day its version was released, which the
    /// changelog's heading for the version gives: a release that moves the
    /// version moves this too.
    #[test]
    fn the_page_is_dated_as_its_version_was_released() {
        const CHANGELOG: &str = include_str!("../../../CHANGELOG.md");
        let heading = format!("## [{}] - ", karyon::VERSION);
        let released = CHANGELOG
            .lines()
            .find_map(|line| line.strip_prefix(heading.as_str()));
        assert_eq!(
            released,
            Some(DATE),
            "CHANGELOG.md dates karyon {} {released:?}, and DATE in src/bin/karyon/manual.rs \
             is {DATE}",
            karyon::VERSION
        );
    }

    #[test]
    fn the_synopsis_has_the_usage_and_every_command() {
        let page = page();
        let synopsis = page
            .split_once(".SH SYNOPSIS\n")
            .and_then(|(_, rest)| rest.split_once("\n.SH "))
            .expect("a synopsis")
            .0;
        assert!(synopsis.starts_with("\\fBkaryon\\fR <PLACE> [FILE | TRACK FILE]..."));
        for (word, _) in commands() {
            assert!(
                synopsis.contains(&format!("\\fBkaryon {word}\\fR")),
                "{word}"
            );
        }
        assert_eq!(synopsis.matches("\n.br\n").count(), commands().len());
    }

    #[test]
    fn text_is_escaped_as_roff_sets_it() {
        assert_eq!(escape("--with-sequence"), "\\-\\-with\\-sequence");
        assert_eq!(escape("aln.cram \\"), "aln.cram \\e");
        assert_eq!(escape("r² at 1×"), "r\\[u00B2] at 1\\[u00D7]");
        assert!(heading(
            "TRACK OPTIONS, each describing the track before it, once"
        ));
        assert!(heading("USAGE"));
        assert!(!heading(
            "The place comes first: a 1-based inclusive locus string"
        ));
        assert!(!heading("A file"));
    }
}
