#![warn(clippy::pedantic)]
#![warn(clippy::nursery)]
#![allow(clippy::cast_lossless)]
#![allow(clippy::inline_always)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::return_self_not_must_use)]
#![allow(clippy::similar_names)]
#![allow(clippy::suspicious_operation_groupings)]
#![allow(clippy::wildcard_imports)]
#![allow(clippy::unreadable_literal)]
#![allow(clippy::should_panic_without_expect)]

pub mod common;
pub mod consts;
pub mod ffi;
pub mod low_level;
pub mod qcd;
pub mod types;
pub mod utils;
pub mod ym;

pub use consts::I;
pub use types::{Integral, Num, C, NCTYPE, R};
