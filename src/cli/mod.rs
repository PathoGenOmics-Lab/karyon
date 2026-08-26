//! The command line, as a library: the grammar and the walk from it to a figure.
//!
//! This is the same code the `karyon` binary runs, and it lives here rather
//! than beside `main` for two reasons.
//!
//! # Somebody other than a shell may want to drive it
//!
//! The grammar is the crate's one declarative surface. Everything else is Rust
//! calls, which is right for a library and is no use at all to a browser, an
//! editor or a notebook that would like to hand karyon a line of text and get
//! a figure back. The library compiles to `wasm32-unknown-unknown` unchanged,
//! and what stopped a browser reaching it was that the grammar was locked
//! inside a binary target.
//!
//! # Nothing here opens a file
//!
//! [`build`] takes the text, rather than fetching it. Every `--flag <PATH>` in
//! an [`args::Invocation`] is resolved by a closure the caller
//! supplies, so a shell hands it `fs::read_to_string`, a browser hands it a
//! lookup into a map of editor buffers, and a test hands it a literal. That
//! keeps the rule the rest of the crate already follows: the library takes
//! text and returns text, and knowing where a path leads is the caller's job.
//!
//! ```
//! use karyon::cli::{args, build};
//!
//! let words: Vec<String> = "chr1:1-60 --coverage depth.txt --label depth"
//!     .split_whitespace()
//!     .map(String::from)
//!     .collect();
//!
//! let args::Request::Draw(invocation) = args::parse(&words)? else {
//!     unreachable!("that command line draws something")
//! };
//!
//! // No disk anywhere: the one source it names is answered from memory.
//! let svg = build(&invocation, |_source| Ok("chr1\t1\t60\t7\n".to_string()))?;
//! assert!(svg.starts_with("<svg"));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub mod args;
pub mod stack;

pub use crate::cli::args::{parse, Invocation, Request};
pub use crate::cli::stack::{build, BuildError};
