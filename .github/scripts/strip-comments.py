#!/usr/bin/env python3
"""Take the comments out of the site's own CSS and JavaScript, at publish time.

This project writes its reasons down. `docs/stylesheets/extra.css` carries
twenty-two kilobytes of block comments saying why each rule is the way it is,
and that is the point of it: the next person to touch a rule should find the
measurement that put it there. But a comment is worth nothing to a browser, and
extra.css is render blocking on every page, so every reader was waiting on eight
kilobytes over the wire of prose addressed to somebody else.

So the source keeps its comments and the published copy does not. Nothing here
minifies: no renaming, no reordering, no whitespace games, no attempt to be
clever about what a rule means. It removes comments and collapses the blank runs
they leave behind, and that is all, because anything more is a change to the
site that nobody reviewed.

Only files this project wrote are touched. The theme's own assets are already
minified and are not ours to rewrite.

Run after `mkdocs build`, against the built directory:

    python3 .github/scripts/strip-comments.py site
"""

import re
import sys
from pathlib import Path

# Written by this project, and therefore ours to rewrite. The theme's own
# bundles live under assets/stylesheets and assets/javascripts and are left
# exactly as they arrived.
OURS = ("stylesheets/extra.css", "stylesheets/landing.css", "stylesheets/playground.css",
        "assets/figure-viewer.js", "assets/karyon-live.js", "assets/karyon-wasm.js",
        "assets/playground.js")


def strip_css(text):
    """Remove /* */ comments, skipping over strings so a comment marker inside
    one is left alone."""
    out = []
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        if c in "'\"":
            j = i + 1
            while j < n:
                if text[j] == "\\":
                    j += 2
                    continue
                if text[j] == c:
                    j += 1
                    break
                j += 1
            out.append(text[i:j])
            i = j
        elif text.startswith("/*", i):
            end = text.find("*/", i + 2)
            i = n if end < 0 else end + 2
        else:
            out.append(c)
            i += 1
    return "".join(out)


def strip_js(text):
    """The same, plus line comments, skipping strings, template literals and
    regular expression literals so nothing that looks like a comment inside one
    is removed."""
    out = []
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        if c in "'\"`":
            j = i + 1
            while j < n:
                if text[j] == "\\":
                    j += 2
                    continue
                if text[j] == c:
                    j += 1
                    break
                j += 1
            out.append(text[i:j])
            i = j
        elif text.startswith("/*", i):
            end = text.find("*/", i + 2)
            i = n if end < 0 else end + 2
        elif text.startswith("//", i):
            end = text.find("\n", i)
            i = n if end < 0 else end
        elif c == "/":
            # A slash is a regular expression only where a value may begin. The
            # previous non-space character decides, and a division follows a
            # value while a literal follows an operator or a bracket.
            prev = "".join(out).rstrip()
            last = prev[-1] if prev else ""
            if last and (last.isalnum() or last in ")]}_$"):
                out.append(c)
                i += 1
                continue
            j = i + 1
            klass = False
            while j < n:
                if text[j] == "\\":
                    j += 2
                    continue
                if text[j] == "[":
                    klass = True
                elif text[j] == "]":
                    klass = False
                elif text[j] == "/" and not klass:
                    j += 1
                    break
                elif text[j] == "\n":
                    break
                j += 1
            out.append(text[i:j])
            i = j
        else:
            out.append(c)
            i += 1
    return "".join(out)


def tidy(text):
    text = re.sub(r"[ \t]+$", "", text, flags=re.M)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip() + "\n"


def main(root):
    root = Path(root)
    saved = 0
    for rel in OURS:
        path = root / rel
        if not path.exists():
            print(f"  skipped, not built: {rel}")
            continue
        before = path.read_text(encoding="utf-8")
        after = tidy(strip_js(before) if rel.endswith(".js") else strip_css(before))
        path.write_text(after, encoding="utf-8")
        saved += len(before) - len(after)
        print(f"  {rel}: {len(before)} -> {len(after)} bytes")
    print(f"  {saved} bytes of comments removed from the published copy")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "site")
