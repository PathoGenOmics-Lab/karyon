#!/usr/bin/env python3
"""Fail if a tracked file matches no rule in .github/labeler.yml.

The labels only answer "which part of the project does this pull request touch"
if every part of the project is in the file, and it has twice not been.
`.github/labeler.yml` itself was uncovered until the rule at the bottom of it
was swept, and `src/cli/`, which is 93 per cent of the lines of the command
line, was uncovered from the day it was split out of the binary on 2026-08-26
until 2026-08-29. Both were found by somebody noticing an unlabelled pull
request, which is a slow way to find out and not one that happens every time.

So the coverage is checked by running the globs rather than by reading them.
Two files are allowed to match nothing, and they are named below with the
reason; anything else is a failure that names the file.

The shape is checked too, and strictly. actions/labeler took a bare list of
globs under the label name in v4 and has wanted `changed-files` since v5, and
the old shape does not fail: it matches nothing, quietly, for as long as nobody
looks. Anything that is not exactly the shape below is a parse error here.

    label:
      - changed-files:
          - any-glob-to-any-file:
              - "some/glob/**"

No YAML library, because the crate has no dependencies and a lint for it should
not need one either. The parser accepts that shape and refuses everything else,
which is the check rather than a shortcut around it.

    python3 .github/scripts/labeler-coverage.py
"""

import re
import subprocess
import sys
from pathlib import Path

CONFIG = Path(".github/labeler.yml")

# Tracked files that are meant to match nothing. Neither is a part of the
# project anybody filters the pull request list for, and a label that lied
# about them would be worse than no label.
ALLOWED_BARE = {".gitignore", "LICENSE"}

# The matchers v5 and later understand. Only the first is used here, and the
# rest are listed so that reaching for one is not mistaken for a typo.
MATCHERS = (
    "any-glob-to-any-file",
    "all-globs-to-all-files",
    "any-glob-to-all-files",
    "all-globs-to-any-file",
)


def parse(text):
    """The globs each label carries, refusing anything but the expected shape."""
    rules = {}
    label = None
    seen_changed_files = False
    for number, line in enumerate(text.split("\n"), 1):
        bare = line.split("#")[0].rstrip() if not line.lstrip().startswith("#") else ""
        if not bare.strip():
            continue
        indent = len(bare) - len(bare.lstrip())
        body = bare.strip()

        if indent == 0:
            if not body.endswith(":"):
                raise SystemExit(
                    f"{CONFIG}:{number}: expected a label name, found {body!r}"
                )
            label = body[:-1]
            rules[label] = []
            seen_changed_files = False
        elif body == "- changed-files:":
            if label is None:
                raise SystemExit(f"{CONFIG}:{number}: changed-files outside a label")
            seen_changed_files = True
        elif body.startswith("- ") and body.endswith(":"):
            matcher = body[2:-1]
            if not seen_changed_files:
                raise SystemExit(
                    f"{CONFIG}:{number}: {matcher!r} is not under a changed-files key. "
                    "That is the v4 shape and v5 and later match nothing with it."
                )
            if matcher not in MATCHERS:
                raise SystemExit(f"{CONFIG}:{number}: unknown matcher {matcher!r}")
        elif body.startswith('- "') and body.endswith('"'):
            if not seen_changed_files:
                raise SystemExit(
                    f"{CONFIG}:{number}: {body[2:]} is a bare glob under {label!r}. "
                    "That is the v4 shape and v5 and later match nothing with it."
                )
            rules[label].append(body[3:-1])
        else:
            raise SystemExit(f"{CONFIG}:{number}: cannot read {body!r}")
    return rules


def matcher(glob):
    """One glob as a regular expression, with minimatch's reading of the stars.

    Two stars cross a directory separator and one does not, which is the whole
    of why `src/track/*.rs` leaves `src/track/tree/` to another label and
    `src/read/**` takes all of its module.
    """
    out, i, size = "", 0, len(glob)
    while i < size:
        if glob.startswith("**/", i):
            out += "(?:.*/)?"
            i += 3
        elif glob.startswith("**", i):
            out += ".*"
            i += 2
        elif glob[i] == "*":
            out += "[^/]*"
            i += 1
        elif glob[i] == "?":
            out += "[^/]"
            i += 1
        else:
            out += re.escape(glob[i])
            i += 1
    return re.compile("^" + out + "$")


def main():
    if not CONFIG.exists():
        raise SystemExit(f"{CONFIG} is not there, so nothing labels a pull request")
    rules = parse(CONFIG.read_text(encoding="utf-8"))
    patterns = {label: [matcher(g) for g in globs] for label, globs in rules.items()}

    tracked = subprocess.run(
        ["git", "ls-files"], capture_output=True, text=True, check=True
    ).stdout.split()

    bare = [
        path
        for path in tracked
        if path not in ALLOWED_BARE
        and not any(p.match(path) for ps in patterns.values() for p in ps)
    ]

    for label in rules:
        hit = [path for path in tracked if any(p.match(path) for p in patterns[label])]
        print(f"  {label:15s} {len(hit):4d} files")
    print(f"  {len(tracked)} tracked, {len(ALLOWED_BARE)} deliberately bare")

    if bare:
        head = f"\n{len(bare)} tracked file(s) match no rule in {CONFIG}:"
        print(head, file=sys.stderr)
        for path in sorted(bare):
            print(f"  {path}", file=sys.stderr)
        print(
            "\nGive them a label, or add them to ALLOWED_BARE in this script with the\n"
            "reason. A pull request that touches only these arrives unlabelled.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
