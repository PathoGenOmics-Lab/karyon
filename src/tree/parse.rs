//! Reading a phylogeny: Newick, annotated Newick and NEXUS.
//!
//! Turning text into a [`Tree`](super::Tree) is a different job from anything
//! done to one afterwards, and it is the job with all the ambiguity in it. The
//! format writes support values and internal names in the same place, brackets
//! carry BEAST and NHX annotations in one dialect and comments in another, and
//! a quoted label may contain any of the punctuation the grammar uses. Every
//! one of those decisions is here rather than spread through the operations
//! that assume they were already made.

use super::*;

pub(super) fn parse_newick_impl(input: &str, preserve_annotations: bool) -> Result<Tree, Error> {
    let text = input.trim().trim_end_matches(';').trim();
    if text.is_empty() {
        return Err(Error::InvalidNewick {
            reason: "empty tree",
        });
    }

    let mut nodes: Vec<Clade> = Vec::new();
    let mut annotations: Vec<Annotations> = Vec::new();
    let mut tree_annotations = Annotations::new();
    let mut rooted = None;
    let mut stack: Vec<usize> = Vec::new();
    let mut current: Option<usize> = None;
    let mut chars = text.chars().peekable();
    // One buffer each for the pieces that are read and thrown away, rather
    // than one per node. A million tip tree is two million nodes, so a `String`
    // per branch length and per bracket is two million allocations that live
    // for the length of a number.
    let mut number = String::new();
    let mut comment = String::new();

    while let Some(c) = chars.next() {
        match c {
            '(' => {
                let parent = stack.last().copied();
                if parent.is_none() && !nodes.is_empty() {
                    return Err(Error::InvalidNewick {
                        reason: "more than one root",
                    });
                }
                let index = add_node(&mut nodes, &mut annotations, parent);
                stack.push(index);
                current = None;
            }
            ')' => {
                if let (None, Some(parent)) = (current, stack.last().copied()) {
                    add_node(&mut nodes, &mut annotations, Some(parent));
                }
                let closed = stack.pop().ok_or(Error::InvalidNewick {
                    reason: "unbalanced parentheses",
                })?;
                current = Some(closed);
            }
            ',' => {
                if stack.is_empty() {
                    return Err(Error::InvalidNewick {
                        reason: "comma outside any clade",
                    });
                }
                if current.is_none() {
                    add_node(&mut nodes, &mut annotations, stack.last().copied());
                }
                current = None;
            }
            ':' => {
                number.clear();
                while let Some(next) = chars.peek() {
                    if next.is_ascii_digit() || matches!(next, '.' | '-' | '+' | 'e' | 'E') {
                        number.push(*next);
                        chars.next();
                    } else {
                        break;
                    }
                }
                let length = number.parse::<f64>().map_err(|_| Error::InvalidNewick {
                    reason: "branch length is not a number",
                })?;
                let target = match current {
                    Some(index) => index,
                    None => {
                        let parent = stack.last().copied().ok_or(Error::InvalidNewick {
                            reason: "branch length with nothing to attach to",
                        })?;
                        let index = add_node(&mut nodes, &mut annotations, Some(parent));
                        current = Some(index);
                        index
                    }
                };
                nodes[target].branch_length = Some(length);
            }
            '[' => {
                comment.clear();
                for next in chars.by_ref() {
                    if next == ']' {
                        break;
                    }
                    comment.push(next);
                }
                if preserve_annotations {
                    let (root_marker, fields) = parse_comment(&comment);
                    if let Some(value) = root_marker {
                        rooted = Some(value);
                    }
                    if let Some(node) = current {
                        annotations[node].extend(fields);
                    } else {
                        tree_annotations.extend(fields);
                    }
                }
            }
            c if c.is_whitespace() => {}
            _ => {
                let mut name = String::new();
                let quoted = c == '\'' || c == '"';
                if !quoted {
                    name.push(c);
                }
                let quote = c;
                while let Some(next) = chars.peek() {
                    if quoted {
                        if *next == quote {
                            chars.next();
                            if chars.peek() == Some(&quote) {
                                name.push(quote);
                                chars.next();
                                continue;
                            }
                            break;
                        }
                        name.push(*next);
                        chars.next();
                    } else if matches!(*next, '(' | ')' | ',' | ':' | ';' | '[') {
                        break;
                    } else {
                        name.push(*next);
                        chars.next();
                    }
                }
                let name = name.trim().to_string();

                match current {
                    Some(index) if !nodes[index].is_leaf() => match name.parse::<f64>() {
                        Ok(support) => nodes[index].support = Some(support),
                        Err(_) => nodes[index].name = Some(name),
                    },
                    _ => {
                        let parent = stack.last().copied();
                        if parent.is_none() && !nodes.is_empty() {
                            return Err(Error::InvalidNewick {
                                reason: "more than one root",
                            });
                        }
                        let index = add_node(&mut nodes, &mut annotations, parent);
                        nodes[index].name = Some(name);
                        current = Some(index);
                    }
                }
            }
        }
    }

    if !stack.is_empty() {
        return Err(Error::InvalidNewick {
            reason: "unbalanced parentheses",
        });
    }
    if nodes.is_empty() {
        return Err(Error::InvalidNewick {
            reason: "empty tree",
        });
    }
    Ok(Tree {
        nodes,
        root: 0,
        annotations,
        tree_annotations,
        rooted,
    })
}

pub(super) fn add_node(
    nodes: &mut Vec<Clade>,
    annotations: &mut Vec<Annotations>,
    parent: Option<usize>,
) -> usize {
    nodes.push(Clade {
        name: None,
        branch_length: None,
        support: None,
        children: Vec::new(),
        parent,
    });
    annotations.push(Annotations::new());
    let index = nodes.len() - 1;
    if let Some(parent) = parent {
        nodes[parent].children.push(index);
    }
    index
}

pub(super) fn parse_comment(comment: &str) -> (Option<bool>, Annotations) {
    let text = comment.trim();
    if text.eq_ignore_ascii_case("&R") {
        return (Some(true), Annotations::new());
    }
    if text.eq_ignore_ascii_case("&U") {
        return (Some(false), Annotations::new());
    }

    let mut annotations = Annotations::new();
    if let Some(body) = text.strip_prefix("&&NHX:") {
        for field in split_delimited(body, ':') {
            insert_annotation(&mut annotations, &field, '=');
        }
    } else if let Some(body) = text.strip_prefix('&') {
        for field in split_delimited(body, ',') {
            insert_annotation(&mut annotations, &field, '=');
        }
    } else if !text.is_empty() {
        annotations.insert(
            "comment".to_string(),
            AnnotationValue::Text(text.to_string()),
        );
    }
    (None, annotations)
}

pub(super) fn insert_annotation(annotations: &mut Annotations, field: &str, separator: char) {
    let field = field.trim();
    if field.is_empty() {
        return;
    }
    let (key, value) = field
        .split_once(separator)
        .map_or((field, "true"), |(key, value)| (key.trim(), value.trim()));
    if !key.is_empty() {
        annotations.insert(key.to_string(), parse_annotation_value(value));
    }
}

pub(super) fn parse_annotation_value(value: &str) -> AnnotationValue {
    let value = value.trim();
    if value.len() >= 2 {
        let first = value.as_bytes()[0];
        let last = value.as_bytes()[value.len() - 1];
        if (first == b'\'' && last == b'\'') || (first == b'"' && last == b'"') {
            return AnnotationValue::Text(value[1..value.len() - 1].to_string());
        }
        if first == b'{' && last == b'}' {
            return AnnotationValue::List(
                split_delimited(&value[1..value.len() - 1], ',')
                    .into_iter()
                    .map(|item| parse_annotation_value(&item))
                    .collect(),
            );
        }
    }
    if value.eq_ignore_ascii_case("true") {
        return AnnotationValue::Boolean(true);
    }
    if value.eq_ignore_ascii_case("false") {
        return AnnotationValue::Boolean(false);
    }
    if let Ok(number) = value.parse::<f64>() {
        if number.is_finite() {
            return AnnotationValue::Number(number);
        }
    }
    AnnotationValue::Text(value.to_string())
}

pub(super) fn split_delimited(input: &str, delimiter: char) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quote = None;
    let mut braces = 0usize;
    for character in input.chars() {
        if let Some(active) = quote {
            field.push(character);
            if character == active {
                quote = None;
            }
            continue;
        }
        match character {
            '\'' | '"' => {
                quote = Some(character);
                field.push(character);
            }
            '{' => {
                braces += 1;
                field.push(character);
            }
            '}' => {
                braces = braces.saturating_sub(1);
                field.push(character);
            }
            value if value == delimiter && braces == 0 => {
                fields.push(field.trim().to_string());
                field.clear();
            }
            _ => field.push(character),
        }
    }
    if !field.trim().is_empty() {
        fields.push(field.trim().to_string());
    }
    fields
}

pub(super) fn parse_nexus(input: &str) -> Result<Tree, Error> {
    let statements = nexus_statements(input);
    let mut translation = BTreeMap::new();
    let mut expression = None;

    for statement in &statements {
        let trimmed = statement.trim();
        let lower = trimmed.to_ascii_lowercase();
        if lower.starts_with("translate") {
            let body = trimmed.get("translate".len()..).unwrap_or_default().trim();
            for entry in split_delimited(body, ',') {
                let split = entry.find(char::is_whitespace).ok_or(Error::InvalidNexus {
                    reason: "a translate entry has no taxon name",
                })?;
                let key = entry[..split].trim();
                let name = unquote(entry[split..].trim());
                if key.is_empty() || name.is_empty() {
                    return Err(Error::InvalidNexus {
                        reason: "an empty translate entry",
                    });
                }
                translation.insert(key.to_string(), name);
            }
        } else if lower.starts_with("tree ") || lower.starts_with("utree ") {
            expression = trimmed
                .split_once('=')
                .map(|(_, tree)| tree.trim().to_string());
            if expression.is_none() {
                return Err(Error::InvalidNexus {
                    reason: "a tree statement has no equals sign",
                });
            }
            break;
        }
    }

    let expression = expression.ok_or(Error::InvalidNexus {
        reason: "no tree statement",
    })?;
    let mut tree = Tree::parse_annotated_newick(&expression)?;
    for leaf in tree.leaves() {
        let Some(name) = tree.nodes[leaf].name.as_deref() else {
            continue;
        };
        if let Some(translated) = translation.get(name) {
            tree.nodes[leaf].name = Some(translated.clone());
        }
    }
    Ok(tree)
}

pub(super) fn nexus_statements(input: &str) -> Vec<String> {
    let mut statements = Vec::new();
    let mut statement = String::new();
    let mut quote = None;
    let mut bracket_depth = 0usize;
    for character in input.chars() {
        if let Some(active) = quote {
            statement.push(character);
            if character == active {
                quote = None;
            }
            continue;
        }
        match character {
            '\'' | '"' => {
                quote = Some(character);
                statement.push(character);
            }
            '[' => {
                bracket_depth += 1;
                statement.push(character);
            }
            ']' => {
                bracket_depth = bracket_depth.saturating_sub(1);
                statement.push(character);
            }
            ';' if bracket_depth == 0 => {
                if !statement.trim().is_empty() {
                    statements.push(statement.trim().to_string());
                }
                statement.clear();
            }
            _ => statement.push(character),
        }
    }
    if !statement.trim().is_empty() {
        statements.push(statement.trim().to_string());
    }
    statements
}

pub(super) fn unquote(value: &str) -> String {
    let value = value.trim();
    if value.len() >= 2 {
        let first = value.as_bytes()[0];
        let last = value.as_bytes()[value.len() - 1];
        if (first == b'\'' && last == b'\'') || (first == b'"' && last == b'"') {
            return value[1..value.len() - 1].replace("''", "'");
        }
    }
    value.to_string()
}
