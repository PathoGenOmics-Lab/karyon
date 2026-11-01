# Fonts

The three faces a PDF from karyon sets its text in, embedded in the file as a
subset of the letters each figure uses:

| File | Sets | Drawn to the widths of |
|:--|:--|:--|
| `LiberationSans-Regular.ttf` | words | Helvetica and Arial |
| `LiberationSans-Bold.ttf` | titles, track names and other bold text | Helvetica-Bold and Arial Bold |
| `LiberationMono-Regular.ttf` | coordinates, sequence and accessions | Courier and Courier New |

They are the Liberation Fonts, version 2.00.1, designed by Steve Matteson and
maintained by Red Hat; `AUTHORS` names them. They are copied byte for byte from
the `fontLiberation` package for R, version 0.1.0, which ships that release,
with its `LICENSE` and `AUTHORS` beside them, and nothing here changes them.
Their SHA-256 sums, to compare with any other copy of the release:

```text
e5b0af421ea2bfbc1ac8d251d647268087ae82786234c57f757d1f0b90fa8b49  LiberationSans-Regular.ttf
d723d5a272970aedf296ef6fc628180df6074bce7769701ea9e0d222c052668c  LiberationSans-Bold.ttf
a9b21391536aec1c8fad37f2d5f24750e7f2d63cd86bccbb2463b6e17005f52e  LiberationMono-Regular.ttf
```

## Licence

The fonts are under the SIL Open Font License, version 1.1, which is in
`LICENSE` with the copyright notices of Google and Red Hat and the Reserved
Font Names Arimo, Tinos, Cousine and Liberation. The rest of karyon is MIT, so
the crate as a whole is `MIT AND OFL-1.1`. The OFL lets a font be bundled with
any software as long as its notice and licence go with it, which this folder
does, and it says that its terms do not reach a document made with the font:
a PDF karyon writes is yours, under whatever terms you choose.
