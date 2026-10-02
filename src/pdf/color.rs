//! CSS colours, in every spelling a person or another tool writes one.
//!
//! The tracks write `#rrggbb` and nothing else, but the colours themselves come
//! from outside: a theme is public, and the command line takes `--color` as any
//! colour a browser would paint, `blue` or `rgb(0, 0, 255)` as much as
//! `#0000ff`. An SVG carries the string as it was given and the browser reads
//! it; a PDF has only numbers, so this module has to be the browser. It reads
//! the four hexadecimal forms, `rgb()` and `hsl()` with or without their alpha
//! and with commas or without, and the 148 names of CSS Color 4.
//!
//! [`theme::parse_hex`](crate::theme) stays as it is: it answers a narrower
//! question, which shades a theme can mix from a colour, and a theme mixes only
//! from six digits.

/// A colour as a PDF writes it: three channels and an opacity, each 0 to 1.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Color {
    pub(crate) rgb: [f64; 3],
    pub(crate) alpha: f64,
}

impl Color {
    pub(crate) const BLACK: Color = Color {
        rgb: [0.0, 0.0, 0.0],
        alpha: 1.0,
    };

    fn opaque(hex: u32) -> Color {
        Color {
            rgb: [
                f64::from((hex >> 16) & 0xff) / 255.0,
                f64::from((hex >> 8) & 0xff) / 255.0,
                f64::from(hex & 0xff) / 255.0,
            ],
            alpha: 1.0,
        }
    }
}

/// The colour `text` names, or `None` when it names none.
pub(crate) fn parse(text: &str) -> Option<Color> {
    let text = text.trim();
    if let Some(hex) = text.strip_prefix('#') {
        return parse_hex(hex);
    }
    let lower = text.to_ascii_lowercase();
    if let Some(open) = lower.find('(') {
        let inside = lower[open + 1..].strip_suffix(')')?;
        return match lower[..open].trim() {
            "rgb" | "rgba" => parse_rgb(inside),
            "hsl" | "hsla" => parse_hsl(inside),
            _ => None,
        };
    }
    if lower == "transparent" {
        return Some(Color {
            rgb: [0.0, 0.0, 0.0],
            alpha: 0.0,
        });
    }
    NAMED
        .binary_search_by(|(name, _)| (*name).cmp(lower.as_str()))
        .ok()
        .map(|at| Color::opaque(NAMED[at].1))
}

fn parse_hex(hex: &str) -> Option<Color> {
    if !hex.bytes().all(|b| b.is_ascii_hexdigit()) {
        return None;
    }
    let digit = |at: usize| u32::from_str_radix(&hex[at..=at], 16).ok();
    let pair = |at: usize| u32::from_str_radix(&hex[at..at + 2], 16).ok();
    let (channels, alpha) = match hex.len() {
        3 | 4 => {
            let short = |at: usize| digit(at).map(|d| d * 17);
            let channels = [short(0)?, short(1)?, short(2)?];
            let alpha = if hex.len() == 4 { short(3)? } else { 255 };
            (channels, alpha)
        }
        6 | 8 => {
            let channels = [pair(0)?, pair(2)?, pair(4)?];
            let alpha = if hex.len() == 8 { pair(6)? } else { 255 };
            (channels, alpha)
        }
        _ => return None,
    };
    Some(Color {
        rgb: channels.map(|c| f64::from(c) / 255.0),
        alpha: f64::from(alpha) / 255.0,
    })
}

/// The arguments of a colour function: three and an optional alpha, written
/// `a, b, c, d` in the old syntax and `a b c / d` in the new one.
fn arguments(inside: &str) -> Option<([&str; 3], Option<&str>)> {
    let pieces: Vec<&str> = if inside.contains(',') {
        inside.split(',').map(str::trim).collect()
    } else {
        let (colour, alpha) = match inside.split_once('/') {
            Some((colour, alpha)) => (colour, Some(alpha.trim())),
            None => (inside, None),
        };
        let mut pieces: Vec<&str> = colour.split_whitespace().collect();
        pieces.extend(alpha);
        pieces
    };
    match pieces.as_slice() {
        [a, b, c] => Some(([a, b, c], None)),
        [a, b, c, d] => Some(([a, b, c], Some(d))),
        _ => None,
    }
}

/// A number, or a percentage of `whole`.
fn amount(text: &str, whole: f64) -> Option<f64> {
    let value = match text.strip_suffix('%') {
        Some(percent) => percent.trim().parse::<f64>().ok()? / 100.0 * whole,
        None => text.parse::<f64>().ok()?,
    };
    value.is_finite().then_some(value)
}

fn alpha(text: Option<&str>) -> Option<f64> {
    match text {
        Some(text) => Some(amount(text, 1.0)?.clamp(0.0, 1.0)),
        None => Some(1.0),
    }
}

fn parse_rgb(inside: &str) -> Option<Color> {
    let ([r, g, b], a) = arguments(inside)?;
    let channel = |text: &str| Some((amount(text, 255.0)? / 255.0).clamp(0.0, 1.0));
    Some(Color {
        rgb: [channel(r)?, channel(g)?, channel(b)?],
        alpha: alpha(a)?,
    })
}

fn parse_hsl(inside: &str) -> Option<Color> {
    let ([h, s, l], a) = arguments(inside)?;
    let hue = {
        let (number, per_turn) = if let Some(n) = h.strip_suffix("deg") {
            (n, 360.0)
        } else if let Some(n) = h.strip_suffix("grad") {
            (n, 400.0)
        } else if let Some(n) = h.strip_suffix("rad") {
            (n, std::f64::consts::TAU)
        } else if let Some(n) = h.strip_suffix("turn") {
            (n, 1.0)
        } else {
            (h, 360.0)
        };
        let turns = number.trim().parse::<f64>().ok()? / per_turn;
        if !turns.is_finite() {
            return None;
        }
        turns.rem_euclid(1.0) * 6.0
    };
    // CSS Color 4 lets saturation and lightness go without their percent sign.
    let part = |text: &str| -> Option<f64> {
        let value = match text.strip_suffix('%') {
            Some(percent) => percent.trim().parse::<f64>().ok()?,
            None => text.parse::<f64>().ok()?,
        };
        value.is_finite().then_some((value / 100.0).clamp(0.0, 1.0))
    };
    let (s, l) = (part(s)?, part(l)?);
    // The conversion CSS Color 4 gives in its sample code.
    let k = |n: f64| (n + hue * 2.0) % 12.0;
    let chroma = s * l.min(1.0 - l);
    let f = |n: f64| l - chroma * (k(n) - 3.0).min(9.0 - k(n)).clamp(-1.0, 1.0);
    Some(Color {
        rgb: [f(0.0), f(8.0), f(4.0)],
        alpha: alpha(a)?,
    })
}

/// The named colours of CSS Color 4, sorted for a binary search.
const NAMED: [(&str, u32); 148] = [
    ("aliceblue", 0xf0f8ff),
    ("antiquewhite", 0xfaebd7),
    ("aqua", 0x00ffff),
    ("aquamarine", 0x7fffd4),
    ("azure", 0xf0ffff),
    ("beige", 0xf5f5dc),
    ("bisque", 0xffe4c4),
    ("black", 0x000000),
    ("blanchedalmond", 0xffebcd),
    ("blue", 0x0000ff),
    ("blueviolet", 0x8a2be2),
    ("brown", 0xa52a2a),
    ("burlywood", 0xdeb887),
    ("cadetblue", 0x5f9ea0),
    ("chartreuse", 0x7fff00),
    ("chocolate", 0xd2691e),
    ("coral", 0xff7f50),
    ("cornflowerblue", 0x6495ed),
    ("cornsilk", 0xfff8dc),
    ("crimson", 0xdc143c),
    ("cyan", 0x00ffff),
    ("darkblue", 0x00008b),
    ("darkcyan", 0x008b8b),
    ("darkgoldenrod", 0xb8860b),
    ("darkgray", 0xa9a9a9),
    ("darkgreen", 0x006400),
    ("darkgrey", 0xa9a9a9),
    ("darkkhaki", 0xbdb76b),
    ("darkmagenta", 0x8b008b),
    ("darkolivegreen", 0x556b2f),
    ("darkorange", 0xff8c00),
    ("darkorchid", 0x9932cc),
    ("darkred", 0x8b0000),
    ("darksalmon", 0xe9967a),
    ("darkseagreen", 0x8fbc8f),
    ("darkslateblue", 0x483d8b),
    ("darkslategray", 0x2f4f4f),
    ("darkslategrey", 0x2f4f4f),
    ("darkturquoise", 0x00ced1),
    ("darkviolet", 0x9400d3),
    ("deeppink", 0xff1493),
    ("deepskyblue", 0x00bfff),
    ("dimgray", 0x696969),
    ("dimgrey", 0x696969),
    ("dodgerblue", 0x1e90ff),
    ("firebrick", 0xb22222),
    ("floralwhite", 0xfffaf0),
    ("forestgreen", 0x228b22),
    ("fuchsia", 0xff00ff),
    ("gainsboro", 0xdcdcdc),
    ("ghostwhite", 0xf8f8ff),
    ("gold", 0xffd700),
    ("goldenrod", 0xdaa520),
    ("gray", 0x808080),
    ("green", 0x008000),
    ("greenyellow", 0xadff2f),
    ("grey", 0x808080),
    ("honeydew", 0xf0fff0),
    ("hotpink", 0xff69b4),
    ("indianred", 0xcd5c5c),
    ("indigo", 0x4b0082),
    ("ivory", 0xfffff0),
    ("khaki", 0xf0e68c),
    ("lavender", 0xe6e6fa),
    ("lavenderblush", 0xfff0f5),
    ("lawngreen", 0x7cfc00),
    ("lemonchiffon", 0xfffacd),
    ("lightblue", 0xadd8e6),
    ("lightcoral", 0xf08080),
    ("lightcyan", 0xe0ffff),
    ("lightgoldenrodyellow", 0xfafad2),
    ("lightgray", 0xd3d3d3),
    ("lightgreen", 0x90ee90),
    ("lightgrey", 0xd3d3d3),
    ("lightpink", 0xffb6c1),
    ("lightsalmon", 0xffa07a),
    ("lightseagreen", 0x20b2aa),
    ("lightskyblue", 0x87cefa),
    ("lightslategray", 0x778899),
    ("lightslategrey", 0x778899),
    ("lightsteelblue", 0xb0c4de),
    ("lightyellow", 0xffffe0),
    ("lime", 0x00ff00),
    ("limegreen", 0x32cd32),
    ("linen", 0xfaf0e6),
    ("magenta", 0xff00ff),
    ("maroon", 0x800000),
    ("mediumaquamarine", 0x66cdaa),
    ("mediumblue", 0x0000cd),
    ("mediumorchid", 0xba55d3),
    ("mediumpurple", 0x9370db),
    ("mediumseagreen", 0x3cb371),
    ("mediumslateblue", 0x7b68ee),
    ("mediumspringgreen", 0x00fa9a),
    ("mediumturquoise", 0x48d1cc),
    ("mediumvioletred", 0xc71585),
    ("midnightblue", 0x191970),
    ("mintcream", 0xf5fffa),
    ("mistyrose", 0xffe4e1),
    ("moccasin", 0xffe4b5),
    ("navajowhite", 0xffdead),
    ("navy", 0x000080),
    ("oldlace", 0xfdf5e6),
    ("olive", 0x808000),
    ("olivedrab", 0x6b8e23),
    ("orange", 0xffa500),
    ("orangered", 0xff4500),
    ("orchid", 0xda70d6),
    ("palegoldenrod", 0xeee8aa),
    ("palegreen", 0x98fb98),
    ("paleturquoise", 0xafeeee),
    ("palevioletred", 0xdb7093),
    ("papayawhip", 0xffefd5),
    ("peachpuff", 0xffdab9),
    ("peru", 0xcd853f),
    ("pink", 0xffc0cb),
    ("plum", 0xdda0dd),
    ("powderblue", 0xb0e0e6),
    ("purple", 0x800080),
    ("rebeccapurple", 0x663399),
    ("red", 0xff0000),
    ("rosybrown", 0xbc8f8f),
    ("royalblue", 0x4169e1),
    ("saddlebrown", 0x8b4513),
    ("salmon", 0xfa8072),
    ("sandybrown", 0xf4a460),
    ("seagreen", 0x2e8b57),
    ("seashell", 0xfff5ee),
    ("sienna", 0xa0522d),
    ("silver", 0xc0c0c0),
    ("skyblue", 0x87ceeb),
    ("slateblue", 0x6a5acd),
    ("slategray", 0x708090),
    ("slategrey", 0x708090),
    ("snow", 0xfffafa),
    ("springgreen", 0x00ff7f),
    ("steelblue", 0x4682b4),
    ("tan", 0xd2b48c),
    ("teal", 0x008080),
    ("thistle", 0xd8bfd8),
    ("tomato", 0xff6347),
    ("turquoise", 0x40e0d0),
    ("violet", 0xee82ee),
    ("wheat", 0xf5deb3),
    ("white", 0xffffff),
    ("whitesmoke", 0xf5f5f5),
    ("yellow", 0xffff00),
    ("yellowgreen", 0x9acd32),
];

#[cfg(test)]
mod tests {
    use super::*;

    fn rgb255(text: &str) -> Option<[u8; 3]> {
        parse(text).map(|c| c.rgb.map(|v| (v * 255.0).round() as u8))
    }

    #[test]
    fn every_spelling_of_a_colour_is_the_same_colour() {
        for spelling in [
            "blue",
            "BLUE",
            "#00f",
            "#0000ff",
            "#0000FF",
            "rgb(0,0,255)",
            "rgb(0, 0, 255)",
            "rgb(0 0 255)",
            "rgb(0%, 0%, 100%)",
            "rgba(0, 0, 255, 1)",
            "hsl(240, 100%, 50%)",
            "hsl(240deg 100% 50%)",
            "hsl(0.6667turn 100% 50%)",
        ] {
            assert_eq!(rgb255(spelling), Some([0, 0, 255]), "{spelling}");
        }
    }

    #[test]
    fn an_alpha_is_read_in_each_form_that_carries_one() {
        for spelling in [
            "#0000ff80",
            "#00f8",
            "rgba(0, 0, 255, 0.5)",
            "rgb(0 0 255 / 50%)",
            "hsla(240, 100%, 50%, 0.5)",
        ] {
            let alpha = parse(spelling).unwrap().alpha;
            assert!((alpha - 0.5).abs() < 0.04, "{spelling}: {alpha}");
        }
        assert_eq!(parse("transparent").unwrap().alpha, 0.0);
        assert_eq!(parse("#d55e00").unwrap().alpha, 1.0);
    }

    #[test]
    fn what_is_not_a_colour_is_none() {
        for text in [
            "",
            "none",
            "url(#g)",
            "#12",
            "#ggg",
            "rgb(1, 2)",
            "notacolour",
            "lab(50 0 0)",
        ] {
            assert!(parse(text).is_none(), "{text}");
        }
        // The table is sorted, or the search would miss names.
        assert!(NAMED.windows(2).all(|pair| pair[0].0 < pair[1].0));
        assert_eq!(rgb255("rebeccapurple"), Some([0x66, 0x33, 0x99]));
        assert_eq!(rgb255("grey"), rgb255("gray"));
    }
}
