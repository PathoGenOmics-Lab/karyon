//! The scripts `karyon completions` prints, which teach bash, zsh and fish
//! the words of a karyon command line.
//!
//! No flag is written in this file. Each script is built from the entries of
//! the help text, which tests hold to the flags the parser answers to, and
//! from the parser itself: whether a flag reads a file is
//! [`args::takes_a_file`], and the words a flag that takes one of a few
//! accepts are the ones it names when it refuses another. So a flag wired
//! into the parser and written up in the help is completed with no edit
//! here, and one the help leaves out fails the test that already holds the
//! help to the parser.

use karyon::cli::args;

use super::{commands, entries, section, split_entry};

/// The script for `shell`.
pub(crate) fn script(shell: args::Shell) -> String {
    match shell {
        args::Shell::Bash => bash(),
        args::Shell::Zsh => zsh(),
        args::Shell::Fish => fish(),
    }
}

/// What a shell offers after an option.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum Offer {
    /// Nothing more: the option takes no value.
    Nothing,
    /// The names of files, for an option that reads one or writes one.
    Files,
    /// These words, the only ones the option takes.
    Words(Vec<&'static str>),
    /// A value of the reader's own, a name, a number or a colour, which
    /// there is no list of.
    Anything,
}

/// An option of the help, as a shell offers it.
pub(crate) struct Opt {
    /// Every spelling of it, as `-o` and `--output`.
    pub(crate) names: Vec<String>,
    /// The value it takes as its entry names it, `FILE` for `<FILE>`.
    pub(crate) value: Option<String>,
    /// What a shell offers after it.
    pub(crate) offer: Offer,
    /// What it is: the first clause of its entry.
    pub(crate) summary: String,
}

/// Every option the help lists, in its order: the tracks, the options of a
/// track and the figure's.
pub(crate) fn options() -> Vec<Opt> {
    ["TRACKS", "TRACK OPTIONS", "FIGURE OPTIONS"]
        .into_iter()
        .flat_map(|heading| entries(section(heading)))
        .map(|(_, entry)| {
            let (names, value, said) = split_entry(&entry);
            // The parser is asked by the long spelling, which every option has.
            let flag = names.last().copied().unwrap_or_default();
            let offer = match value {
                None => Offer::Nothing,
                Some(_) if args::takes_a_file(flag) => Offer::Files,
                Some(_) => match choices(flag) {
                    words if words.is_empty() => Offer::Anything,
                    words => Offer::Words(words),
                },
            };
            Opt {
                names: names.iter().map(|name| name.to_string()).collect(),
                value: value.map(str::to_string),
                offer,
                summary: summary(&said).to_string(),
            }
        })
        .collect()
}

/// The first clause of what an entry says, short enough to stand beside its
/// word in a list: up to its first colon, semicolon or full stop.
pub(crate) fn summary(said: &str) -> &str {
    let end = [": ", "; ", ". "]
        .iter()
        .filter_map(|stop| said.find(stop))
        .min()
        .unwrap_or(said.len());
    said[..end].trim_end_matches('.')
}

/// The words `flag` takes where it takes one of a few, as the parser says
/// them.
///
/// Given a word it does not know, such an option is refused with the words
/// it does know, as `light or dark`, and each of those the parser then
/// accepts after the option is offered: `or` is refused and `light` is not.
/// An option that takes any word, or only a number, offers none.
pub(crate) fn choices(flag: &str) -> Vec<&'static str> {
    let Err(args::ArgError::BadValue { expected, .. }) = args::parse(&line(None, flag, "?")) else {
        return Vec::new();
    };
    let mut words: Vec<&'static str> = Vec::new();
    for word in expected.split(|c: char| c.is_whitespace() || c == ',') {
        let a_word = word.starts_with(|c: char| c.is_ascii_lowercase())
            && word
                .chars()
                .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-');
        if a_word && !words.contains(&word) && accepts(flag, word) {
            words.push(word);
        }
    }
    words
}

/// Whether the parser takes `word` after `flag`: anywhere, for an option of
/// the figure's, and after some track that has the option, for a track's.
fn accepts(flag: &str, word: &str) -> bool {
    match args::parse(&line(None, flag, word)) {
        // Refused for what it is, which is read before the track is.
        Err(args::ArgError::BadValue { .. }) => false,
        Err(args::ArgError::NoTrackYet(_)) => args::Kind::ALL.iter().any(|kind| {
            !matches!(
                args::parse(&line(Some(*kind), flag, word)),
                Err(args::ArgError::BadValue { .. }
                    | args::ArgError::WrongTrack { .. }
                    | args::ArgError::Misplaced { .. })
            )
        }),
        _ => true,
    }
}

/// A command line giving `word` to `flag`, after a track of `kind` where
/// there is one.
fn line(kind: Option<args::Kind>, flag: &str, word: &str) -> Vec<String> {
    let mut line = vec!["chr1:1-10".to_string()];
    if let Some(kind) = kind {
        line.push(kind.dashed().to_string());
        if !kind.reads_nothing() {
            line.push("x.txt".to_string());
        }
    }
    line.extend([flag.to_string(), word.to_string()]);
    line
}

/// What `karyon help` takes, each with what it is: every track, as its entry
/// in the help says, and all of it.
fn topics() -> Vec<(&'static str, String)> {
    let tracks = entries(section("TRACKS"));
    let mut topics = vec![("all", "every track and every option, in full".to_string())];
    for kind in args::Kind::ALL {
        let said = tracks
            .iter()
            .find(|(flag, _)| *flag == kind.dashed())
            .map(|(_, entry)| split_entry(entry).2)
            .unwrap_or_default();
        topics.push((kind.flag(), summary(&said).to_string()));
    }
    topics
}

/// A command, what it is, and the value after it: what that value is called
/// and the words it is chosen from.
pub(crate) struct Command {
    pub(crate) word: &'static str,
    pub(crate) summary: String,
    pub(crate) value: String,
    pub(crate) next: Vec<(&'static str, String)>,
}

/// The commands of the help, each with what may follow it, by the value its
/// entry names: a topic, as `karyon help` takes, or a shell.
pub(crate) fn command_list() -> Vec<Command> {
    commands()
        .into_iter()
        .map(|(word, entry)| {
            let (_, value, said) = split_entry(&entry);
            let next = match value {
                Some("TOPIC") => topics(),
                Some("SHELL") => args::Shell::ALL
                    .iter()
                    .map(|(_, shell)| (*shell, String::new()))
                    .collect(),
                _ => Vec::new(),
            };
            Command {
                word,
                summary: summary(&said).to_string(),
                value: value.unwrap_or_default().to_ascii_lowercase(),
                next,
            }
        })
        .collect()
}

/// `words` joined by `between` into lines no wider than a terminal, each
/// after `indent` and every one but the last ending in `end`.
fn wrap(words: &[&str], between: &str, indent: &str, end: &str) -> String {
    let mut lines: Vec<String> = Vec::new();
    let mut line = String::new();
    for word in words {
        if !line.is_empty() && indent.len() + line.len() + between.len() + word.len() > 76 {
            lines.push(std::mem::take(&mut line));
        }
        if !line.is_empty() {
            line.push_str(between);
        }
        line.push_str(word);
    }
    lines.push(line);
    let lines: Vec<String> = lines.iter().map(|line| format!("{indent}{line}")).collect();
    lines.join(&format!("{end}\n"))
}

/// Words for `compgen -W`, quoted, and over several lines where one would be
/// wider than a terminal: each line after `indent`, one step further in.
fn listed(words: &[&str], indent: &str) -> String {
    let joined = words.join(" ");
    if indent.len() + joined.len() < 48 {
        return format!("'{joined}'");
    }
    format!(
        "'\n{}\n{indent}'",
        wrap(words, " ", &format!("{indent}    "), "")
    )
}

/// Every spelling of the options whose offer is `wanted`.
fn spellings(options: &[Opt], wanted: impl Fn(&Offer) -> bool) -> Vec<&str> {
    options
        .iter()
        .filter(|option| wanted(&option.offer))
        .flat_map(|option| option.names.iter().map(String::as_str))
        .collect()
}

const BASH: &str = r#"# bash completion for karyon, written by `karyon completions bash` from the
# tables its parser and its help are built from. Load it from ~/.bashrc with
#
#     eval "$(karyon completions bash)"
#
# or save it as ~/.local/share/bash-completion/completions/karyon, where the
# bash-completion package finds it the first time karyon is completed.

# Adds each line compgen prints to the answers, spaces and all.
_karyon_add() {
    local line
    while IFS= read -r line; do
        COMPREPLY+=("$line")
    done < <(compgen "$@")
}

_karyon() {
    COMPREPLY=()
    local cur=${COMP_WORDS[COMP_CWORD]} prev=${COMP_WORDS[COMP_CWORD - 1]}
    # bash splits --theme=dark at the sign, into --theme, = and dark.
    if [[ $cur == "=" ]]; then
        cur=
    elif [[ $prev == "=" && $COMP_CWORD -gt 2 ]]; then
        prev=${COMP_WORDS[COMP_CWORD - 2]}
    fi

    # The value of the option before.
    case $prev in
@VALUES@
    esac

    # A command is the first word, and takes one more at most.
    if [[ $COMP_CWORD -gt 1 ]]; then
        case ${COMP_WORDS[1]} in
@COMMANDS@
        esac
    fi

    if [[ $cur == -* ]]; then
        _karyon_add -W @FLAGS@ -- "$cur"
        return
    fi
    # A command first of all, and anywhere the place or a file named on its
    # own, which is a track.
    if [[ $COMP_CWORD -eq 1 ]]; then
        _karyon_add -W '@FIRST@' -- "$cur"
    fi
    _karyon_add -f -- "$cur"
}

complete -o filenames -F _karyon karyon
"#;

/// The script for bash, a function `complete -F` calls.
fn bash() -> String {
    let options = options();
    // A pattern of several flags is continued onto the next line after its
    // bar, and an arm with no flag is left out, since `)` alone is no pattern.
    let arm = |patterns: &[&str], action: &str| {
        if patterns.is_empty() {
            return String::new();
        }
        let mut arm = wrap(patterns, " | ", "        ", " | \\");
        arm.push_str(")\n");
        if !action.is_empty() {
            arm.push_str(&format!("            {action}\n"));
        }
        arm.push_str("            return\n            ;;\n");
        arm
    };
    let mut values = arm(
        &spellings(&options, |offer| *offer == Offer::Files),
        "_karyon_add -f -- \"$cur\"",
    );
    for option in &options {
        if let Offer::Words(words) = &option.offer {
            let names: Vec<&str> = option.names.iter().map(String::as_str).collect();
            let action = format!(
                "_karyon_add -W {} -- \"$cur\"",
                listed(words, "            ")
            );
            values.push_str(&arm(&names, &action));
        }
    }
    // A value of the reader's own: nothing to offer, and not a file either.
    values.push_str(&arm(
        &spellings(&options, |offer| *offer == Offer::Anything),
        "",
    ));

    let commands = command_list();
    let mut after = String::new();
    for command in &commands {
        after.push_str(&format!("            {})\n", command.word));
        if !command.next.is_empty() {
            let words: Vec<&str> = command.next.iter().map(|(word, _)| *word).collect();
            after.push_str(&format!(
                "                [[ $COMP_CWORD -eq 2 ]] &&\n                    \
                 _karyon_add -W {} -- \"$cur\"\n",
                listed(&words, "                    ")
            ));
        }
        after.push_str("                return\n                ;;\n");
    }

    let first: Vec<&str> = commands.iter().map(|command| command.word).collect();
    BASH.replace("@VALUES@\n", &values)
        .replace("@COMMANDS@\n", &after)
        .replace(
            "@FLAGS@",
            &listed(&spellings(&options, |_| true), "        "),
        )
        .replace("@FIRST@", &first.join(" "))
}

const ZSH: &str = r#"#compdef karyon

# zsh completion for karyon, written by `karyon completions zsh` from the
# tables its parser and its help are built from. Save it as _karyon in a
# folder on $fpath before compinit runs, or load it after compinit with
#
#     source <(karyon completions zsh)

_karyon() {
    # A command is the first word, and takes one more at most.
    if (( CURRENT > 2 )); then
        case $words[2] in
@COMMANDS@
        esac
    fi
    _arguments \
@SPECS@
        '1: :_karyon_first' \
        '*:file:_files'
}

# The first word: a command, or the place, or a file named on its own.
_karyon_first() {
    if (( CURRENT == 2 )); then
        local -a commands=(
@FIRST@
        )
        _describe -t commands command commands
    fi
    _files
}

if [[ $funcstack[1] == _karyon ]]; then
    _karyon "$@"
else
    compdef _karyon karyon
fi
"#;

/// The script for zsh, a completion function for `compinit`, or for
/// `compdef` where it is sourced.
fn zsh() -> String {
    // Inside single quotes a single quote is ended, escaped and reopened.
    let quoted = |text: &str| text.replace('\'', r"'\''");
    let mut specs = String::new();
    for option in options() {
        let explained = if option.summary.is_empty() {
            String::new()
        } else {
            let said = option
                .summary
                .replace('\\', r"\\")
                .replace('[', r"\[")
                .replace(']', r"\]");
            format!("[{said}]")
        };
        for name in &option.names {
            // A long option's value may follow its sign as well as a space.
            let sign = if name.starts_with("--") && option.value.is_some() {
                "="
            } else {
                ""
            };
            let takes = match (&option.value, &option.offer) {
                (None, _) | (_, Offer::Nothing) => String::new(),
                (Some(value), offer) => {
                    let message = value.replace('\\', r"\\").replace(':', r"\:");
                    let action = match offer {
                        Offer::Files => "_files".to_string(),
                        Offer::Words(words) => format!("({})", words.join(" ")),
                        Offer::Anything | Offer::Nothing => " ".to_string(),
                    };
                    format!(":{message}:{action}")
                }
            };
            specs.push_str(&format!(
                "        '*{name}{sign}{}{}' \\\n",
                quoted(&explained),
                quoted(&takes)
            ));
        }
    }

    let described = |word: &str, said: &str| {
        let word = word.replace(':', r"\:");
        if said.is_empty() {
            format!("'{}'", quoted(&word))
        } else {
            format!("'{}:{}'", quoted(&word), quoted(said))
        }
    };
    let mut after = String::new();
    for command in command_list() {
        after.push_str(&format!("            ({})\n", command.word));
        if command.next.is_empty() {
            after.push_str("                return 1\n                ;;\n");
            continue;
        }
        after.push_str("                (( CURRENT == 3 )) || return 1\n");
        after.push_str("                local -a choices=(\n");
        for (word, said) in &command.next {
            after.push_str(&format!("                    {}\n", described(word, said)));
        }
        after.push_str("                )\n");
        after.push_str(&format!(
            "                _describe -t {0} {0} choices\n",
            command.value
        ));
        after.push_str("                return\n                ;;\n");
    }

    let first: String = command_list()
        .iter()
        .map(|command| {
            format!(
                "            {}\n",
                described(command.word, &command.summary)
            )
        })
        .collect();
    ZSH.replace("@COMMANDS@\n", &after)
        .replace("@SPECS@\n", &specs)
        .replace("@FIRST@\n", &first)
}

const FISH: &str = r#"# fish completion for karyon, written by `karyon completions fish` from the
# tables its parser and its help are built from. Save it as
# ~/.config/fish/completions/karyon.fish, where fish finds it the first time
# karyon is completed.

# Whether karyon is the only word so far, so the next one is the first.
function __karyon_first
    test (count (commandline -opc)) -eq 1
end

# Whether the words so far are karyon and one of these.
function __karyon_after
    set -l words (commandline -opc)
    test (count $words) -eq 2
    and contains -- $words[2] $argv
end

# The first word: a command, or the place, or a file named on its own.
@FIRST@
# The word after a command.
@COMMANDS@
# The options, and what each one takes.
@OPTIONS@"#;

/// The script for fish, a list of `complete` commands.
fn fish() -> String {
    // Inside single quotes fish reads a backslash before a quote or another
    // backslash as an escape, and nothing else.
    let quoted = |text: &str| format!("'{}'", text.replace('\\', r"\\").replace('\'', r"\'"));
    let described = |said: &str| {
        if said.is_empty() {
            String::new()
        } else {
            format!(" -d {}", quoted(said))
        }
    };
    let mut first = String::new();
    let mut after = String::new();
    for command in command_list() {
        first.push_str(&format!(
            "complete -c karyon -n __karyon_first -a {}{}\n",
            command.word,
            described(&command.summary)
        ));
        let condition = quoted(&format!("__karyon_after {}", command.word));
        if command.next.is_empty() {
            after.push_str(&format!("complete -c karyon -n {condition} -f\n"));
        }
        for (word, said) in &command.next {
            after.push_str(&format!(
                "complete -c karyon -n {condition} -f -a {word}{}\n",
                described(said)
            ));
        }
    }

    let mut lines = String::new();
    for option in options() {
        let mut line = String::from("complete -c karyon");
        for name in &option.names {
            match name.strip_prefix("--") {
                Some(long) => line.push_str(&format!(" -l {long}")),
                None => line.push_str(&format!(" -s {}", &name[1..])),
            }
        }
        match &option.offer {
            Offer::Nothing => {}
            Offer::Files => line.push_str(" -r -F"),
            Offer::Words(words) => line.push_str(&format!(" -x -a {}", quoted(&words.join(" ")))),
            Offer::Anything => line.push_str(" -x"),
        }
        line.push_str(&described(&option.summary));
        lines.push_str(&line);
        lines.push('\n');
    }
    FISH.replace("@FIRST@\n", &first)
        .replace("@COMMANDS@\n", &after)
        .replace("@OPTIONS@", &lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every flag the parser answers to, short spellings included.
    fn flags() -> Vec<&'static str> {
        args::FLAGS.iter().copied().chain(["-h", "-V"]).collect()
    }

    /// The arms of the bash script's `case $prev`, each flag with what its
    /// arm does before it returns, on one line: nothing, where it offers
    /// nothing.
    fn bash_arms(script: &str) -> Vec<(String, String)> {
        let body = script
            .split_once("    case $prev in\n")
            .expect("the bash script reads the option before")
            .1
            .split_once("\n    esac")
            .expect("and ends its case")
            .0
            .replace(" | \\\n", " | ");
        let mut arms = Vec::new();
        for arm in body.split(";;").filter(|arm| arm.contains(')')) {
            let (patterns, action) = arm.split_once(")\n").expect("an arm has a pattern");
            let action = action
                .split_once("return")
                .expect("an arm returns")
                .0
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ");
            for flag in patterns.split(" | ") {
                arms.push((flag.trim().to_string(), action.clone()));
            }
        }
        arms
    }

    /// The zsh script's spec for one spelling of a flag, without the line's
    /// continuation.
    fn zsh_spec<'a>(script: &'a str, flag: &str) -> Option<&'a str> {
        script
            .lines()
            .map(|line| line.trim().trim_end_matches(" \\"))
            .find(|line| {
                line.strip_prefix(&format!("'*{flag}"))
                    .is_some_and(|rest| rest.starts_with(['=', '[', ':', '\'']))
            })
    }

    /// The fish script's line for one spelling of a flag.
    fn fish_line<'a>(script: &'a str, flag: &str) -> Option<&'a str> {
        let spelled = match flag.strip_prefix("--") {
            Some(long) => format!(" -l {long}"),
            None => format!(" -s {}", &flag[1..]),
        };
        script.lines().find(|line| {
            line.starts_with("complete -c karyon ")
                && line.match_indices(&spelled).any(|(at, _)| {
                    let after = &line[at + spelled.len()..];
                    after.is_empty() || after.starts_with(' ')
                })
        })
    }

    /// Against the parser's own list, so a flag wired up and written up in
    /// the help is completed in every shell, and one left out of the help is
    /// a failing test here as well as there.
    #[test]
    fn every_flag_the_parser_answers_to_is_in_each_script() {
        let bash = script(args::Shell::Bash);
        let zsh = script(args::Shell::Zsh);
        let fish = script(args::Shell::Fish);
        let dashed = bash
            .split_once("if [[ $cur == -* ]]; then")
            .expect("bash offers the flags after a dash")
            .1
            .split_once(" -- \"$cur\"")
            .expect("from one list")
            .0;
        let offered: Vec<&str> = dashed.split_whitespace().collect();
        for flag in flags() {
            assert!(offered.contains(&flag), "bash does not offer {flag}");
            assert!(zsh_spec(&zsh, flag).is_some(), "zsh has no spec for {flag}");
            assert!(
                fish_line(&fish, flag).is_some(),
                "fish does not offer {flag}"
            );
        }
        assert!(flags().len() > 100, "{} flags", flags().len());
    }

    /// After a flag that reads a file, the names of files; after one of a
    /// few words, those words; after a value of the reader's own, nothing,
    /// and not the names of files either.
    #[test]
    fn each_flag_is_followed_by_what_it_takes() {
        let bash = script(args::Shell::Bash);
        let zsh = script(args::Shell::Zsh);
        let fish = script(args::Shell::Fish);
        let arms = bash_arms(&bash);
        let mut checked = (0, 0, 0);
        for option in options() {
            for flag in &option.names {
                let arm = arms
                    .iter()
                    .find(|(name, _)| name == flag)
                    .map(|(_, a)| a.as_str());
                let spec = zsh_spec(&zsh, flag).unwrap_or_else(|| panic!("zsh: {flag}"));
                let line = fish_line(&fish, flag).unwrap_or_else(|| panic!("fish: {flag}"));
                let files = args::takes_a_file(flag);
                assert_eq!(files, option.offer == Offer::Files, "{flag}");
                match &option.offer {
                    Offer::Nothing => {
                        assert_eq!(arm, None, "{flag}");
                        assert!(
                            spec.ends_with(&format!("{flag}'")) || spec.ends_with("]'"),
                            "{spec}"
                        );
                        assert!(!line.contains(" -r") && !line.contains(" -x"), "{line}");
                    }
                    Offer::Files => {
                        assert_eq!(arm, Some("_karyon_add -f -- \"$cur\""), "{flag}");
                        assert!(spec.ends_with(":_files'"), "{spec}");
                        assert!(line.contains(" -r -F"), "{line}");
                        checked.0 += 1;
                    }
                    Offer::Words(words) => {
                        let listed = words.join(" ");
                        // A list too long for its line goes on inside its
                        // quotes, over the lines after it.
                        let flat = |text: &str| {
                            let spaced = text.replace('\'', " ' ");
                            spaced.split_whitespace().collect::<Vec<_>>().join(" ")
                        };
                        let action = format!("_karyon_add -W '{listed}' -- \"$cur\"");
                        assert_eq!(arm.map(flat), Some(flat(&action)), "{flag}");
                        assert!(spec.ends_with(&format!(":({listed})'")), "{spec}");
                        assert!(line.contains(&format!(" -x -a '{listed}'")), "{line}");
                        checked.1 += 1;
                    }
                    Offer::Anything => {
                        assert_eq!(arm, Some(""), "{flag}");
                        assert!(spec.ends_with(": '"), "{spec}");
                        assert!(line.contains(" -x") && !line.contains(" -a "), "{line}");
                        checked.2 += 1;
                    }
                }
            }
        }
        assert!(
            checked.0 >= 40 && checked.1 >= 10 && checked.2 >= 20,
            "{checked:?}"
        );
        // And the words themselves, written out here rather than read from
        // the parser that the scripts read them from.
        let words = |flag: &str| match options()
            .into_iter()
            .find(|option| option.names.iter().any(|name| name == flag))
            .map(|option| option.offer)
        {
            Some(Offer::Words(words)) => words,
            other => panic!("{flag}: {other:?}"),
        };
        assert_eq!(words("--theme"), ["light", "dark"]);
        assert_eq!(words("--aggregate"), ["max", "mean", "min"]);
        assert_eq!(
            words("--format"),
            ["bedgraph", "depth", "values", "bed", "gff3"]
        );
        assert_eq!(words("--threshold"), ["genome-wide"]);
        assert_eq!(words("--max-rows"), ["all"]);
        let styles: Vec<&str> = args::Style::ALL.iter().map(|(_, word)| *word).collect();
        assert_eq!(words("--style"), styles);
        assert!(
            bash.contains("        --theme)\n            _karyon_add -W 'light dark' -- \"$cur\"")
        );
        assert!(zsh.contains("'*--theme=[light or dark]:NAME:(light dark)' \\"));
        assert!(fish.contains("complete -c karyon -l theme -x -a 'light dark' -d 'light or dark'"));
    }

    /// A word an option is offered is one the parser takes after it, and a
    /// word that only shows what a value looks like, as `'#fbfaff'` after
    /// `--background` or `FROM=TO` after `--rename`, is not offered.
    #[test]
    fn an_option_is_offered_the_words_the_parser_takes_and_no_others() {
        for option in options() {
            let flag = option.names.last().unwrap();
            if let Offer::Words(words) = &option.offer {
                for word in words {
                    assert!(accepts(flag, word), "{flag} {word}");
                }
            }
        }
        for flag in [
            "--background",
            "--rename",
            "--colors",
            "--color",
            "--width",
            "--label",
        ] {
            assert!(choices(flag).is_empty(), "{flag}: {:?}", choices(flag));
        }
        assert!(!accepts("--theme", "or"));
        assert!(accepts("--threshold", "genome-wide"));
        assert!(!accepts("--threshold", "genome"));
    }

    /// The first word offers the commands the help lists, and each is one
    /// the parser takes as a command, with the word after it the scripts
    /// offer: a track's name after help, a shell after completions.
    #[test]
    fn the_commands_offered_are_the_parser_s() {
        let commands = command_list();
        let words: Vec<&str> = commands.iter().map(|command| command.word).collect();
        assert_eq!(words, ["help", "completions", "man"]);
        for command in &commands {
            let mut line = vec![command.word.to_string()];
            line.extend(command.next.first().map(|(word, _)| word.to_string()));
            assert!(
                !matches!(args::parse(&line), Ok(args::Request::Draw(_)) | Err(_)),
                "{line:?}"
            );
            for (word, _) in &command.next {
                let line = [command.word.to_string(), word.to_string()];
                assert!(args::parse(&line).is_ok(), "{line:?}");
            }
        }
        let shells: Vec<&str> = args::Shell::ALL.iter().map(|(_, word)| *word).collect();
        let offered: Vec<&str> = commands[1].next.iter().map(|(word, _)| *word).collect();
        assert_eq!(offered, shells);
        assert_eq!(commands[0].next.len(), args::Kind::ALL.len() + 1);
        assert!(commands[2].next.is_empty());
        let zsh = script(args::Shell::Zsh);
        assert!(zsh.contains("            'man:this help as a manual page, set in roff for man'"));
        assert!(zsh.starts_with("#compdef karyon\n"));
    }

    #[test]
    fn a_summary_is_the_first_clause() {
        assert_eq!(
            summary("per-base signal: bedGraph, bigWig"),
            "per-base signal"
        );
        assert_eq!(
            summary("the reference, FASTA or 2bit. A dynseq"),
            "the reference, FASTA or 2bit"
        );
        assert_eq!(summary("genes; a gene is drawn once"), "genes");
        assert_eq!(summary("900 by default"), "900 by default");
        assert_eq!(summary("SJ.out.tab"), "SJ.out.tab");
        assert_eq!(summary(""), "");
    }
}
