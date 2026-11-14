# The picture a shared link shows

`card.html` is the picture a link to the site shows when it is pasted into
Slack, Mastodon, LinkedIn, X or a chat app, and
`docs/assets/social-card.png` is that page photographed at 1200 by 630
pixels. `overrides/main.html` names it in the link preview tags of every page.
The page is made of the site itself: the logo at night, the home page's title
and a strip of `docs/assets/start/reads-dark.svg`, all on the dark page's
ground. Draw it again when any of them, or the look of the site, changes.

MkDocs leaves this folder out of the site: a theme folder's `.html` and
`.md` files are templates and notes, not pages to copy.

From the root of a checkout, with Chrome and a connection for the web fonts
the page loads, as the site does:

```bash
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless \
  --hide-scrollbars --force-device-scale-factor=1 --window-size=1200,630 \
  --virtual-time-budget=10000 \
  --screenshot="$PWD/docs/assets/social-card.png" \
  "file://$PWD/overrides/social/card.html"
```

On Linux the program is `google-chrome` or `chromium`. Then look at the PNG
at its size and at a third of it, about 400 pixels wide, which is how most
previews show it: the title and the line under it have to read there. Keep
the file under 300 KB, since some apps drop a larger picture; the one drawn
on 2026-11-14 is 258 KB.
