//! Calendar dates as decimal years, the unit a dated tree is drawn in.
//!
//! A tip dated `2020-03-15` is a point in time a time axis can place, and a
//! tree annotated that way, as Nextstrain and TreeTime write one, was drawn as
//! a phylogram, since the date was text and not a number. A date becomes the
//! year plus the part of it that has gone by at the start of that day, so
//! `2020-01-01` is `2020.0`, and a date known only to the month is placed in
//! the middle of the month, which is the least wrong place for a date that
//! could be any day of it.

/// Days in each month of a year that is not a leap year.
const DAYS: [u32; 12] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

/// Whether a year has a 29th of February.
fn leap(year: i64) -> bool {
    (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
}

/// How many days `month` of `year` has, the month counted from 1.
fn days_in(year: i64, month: u32) -> u32 {
    if month == 2 && leap(year) {
        29
    } else {
        DAYS[(month - 1) as usize]
    }
}

/// A date as a decimal year: `2020-03-15`, `2020/03/15` or `2020-03`, a
/// month given alone placed in its middle. `None` for text that is not one,
/// or a month or a day that does not exist, as `2021-02-29`.
///
/// ```
/// use karyon::read::date::decimal_year;
///
/// assert_eq!(decimal_year("2020-01-01"), Some(2020.0));
/// assert_eq!(decimal_year("2021-07-02"), Some(2021.0 + 182.0 / 365.0));
/// assert!(decimal_year("2021-02-29").is_none());
/// ```
pub fn decimal_year(text: &str) -> Option<f64> {
    let text = text.trim().trim_matches(|c| c == '"' || c == '\'');
    let parts: Vec<&str> = text.split(['-', '/']).collect();
    let (year, month, day) = match parts.as_slice() {
        [year, month] => (*year, *month, None),
        [year, month, day] => (*year, *month, Some(*day)),
        _ => return None,
    };
    if year.len() != 4 || !year.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let year: i64 = year.parse().ok()?;
    let month: u32 = month
        .parse()
        .ok()
        .filter(|month| (1..=12).contains(month))?;
    let length = if leap(year) { 366.0 } else { 365.0 };
    let before: u32 = (1..month).map(|earlier| days_in(year, earlier)).sum();
    let into = match day {
        Some(day) => {
            let day: u32 = day
                .parse()
                .ok()
                .filter(|day| (1..=days_in(year, month)).contains(day))?;
            f64::from(day - 1)
        }
        None => f64::from(days_in(year, month)) / 2.0,
    };
    Some(year as f64 + (f64::from(before) + into) / length)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_date_is_a_year_and_the_part_of_it_gone_by() {
        assert_eq!(decimal_year("2020-01-01"), Some(2020.0));
        assert_eq!(decimal_year("2020/01/01"), Some(2020.0));
        // A leap year is 366 days long, and its 29th of February is a day.
        assert_eq!(decimal_year("2020-03-01"), Some(2020.0 + 60.0 / 366.0));
        assert!(decimal_year("2020-02-29").is_some());
        assert!(decimal_year("2021-02-29").is_none());
        // A month alone is its middle.
        assert_eq!(decimal_year("2021-02"), Some(2021.0 + 45.0 / 365.0));
        for text in [
            "",
            "2020",
            "20-01-01",
            "2020-13-01",
            "2020-00-10",
            "2020-01-32",
            "May 2020",
        ] {
            assert!(decimal_year(text).is_none(), "{text:?}");
        }
    }
}
