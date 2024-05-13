//! Utility functions.

use crate::{Num, C, I, R};
use std::f64::consts::PI;

/// Finds (at most) one zero of a real or complex function using the
/// [secant method](https://en.wikipedia.org/wiki/Secant_method).
///
/// `f` is the (real or complex) function whose zero is to be searched,
/// `guesses` contains two distinct initial guesses, `tol` is the tolerance
/// requested for the absolute value of the function at its candidate zero,
/// `max_iter` is the maximum number of iterations. After `max_iter` iterations
/// the algorithm is forced to fail and the function returns `None`.
///
/// When the algorithm doesn't fail, `find_root` returns `Some((z0, f0))`, where
/// `z0` is the candidate zero and `f0 = f(z0)`, with `abs(f0)` guaranteed to be
/// smaller than `tol`.
///
/// # Panics
///
/// This function panics if `guesses.0 == guesses.1` or if `tol <= 0.`.
pub fn find_root<T: Num, F: Fn(T) -> T>(
    f: &F,
    guesses: (T, T),
    tol: R,
    max_iter: usize,
) -> Option<(T, T)> {
    assert!(guesses.0 != guesses.1);
    assert!(tol > 0.);

    let (mut z0, mut z1) = guesses;
    let (mut f0, mut f1) = (f(z0), f(z1));
    for _ in 0..max_iter {
        if f1.re().is_nan() || f1.im().is_nan() {
            return None;
        } else if f1.abs() < tol {
            return Some((z1, f1));
        }
        let z2 = (z0 * f1 - z1 * f0) / (f1 - f0);
        (f0, f1) = (f1, f(z2));
        (z0, z1) = (z1, z2);
    }
    None
}

/// Computes the residue of a function `f` at point `z` numerically with
/// an epsilon `eps`.
pub fn compute_residue_with_eps<T: Num, F: Fn(T) -> T>(f: &F, z: T, eps: R) -> T {
    (f(z + eps).inv() - f(z - eps).inv()).inv() * 2. * eps
}

/// Computes the residue of a function `f` at point `z` numerically with
/// an epsilon of `1E-12`.
pub fn compute_residue<T: Num, F: Fn(T) -> T>(f: &F, z: T) -> T {
    compute_residue_with_eps(f, z, 1E-12)
}

/// Computes the spectral function of `f` at point `z` numerically with an
/// epsilon `eps`. The spectral function is defined as
/// `i[f(z+i*eps)-f(z-i*eps)]/2pi`.
pub fn spectral_function_with_eps<F: Fn(C) -> C>(f: F, z: C, eps: R) -> C {
    I * (f(z + I * eps) - f(z - I * eps)) / (2. * PI)
}

/// Computes the spectral function of `f` at point `z` numerically with an
/// epsilon of `1E-12`. The spectral function is defined as
/// `i[f(z+i*eps)-f(z-i*eps)]/2pi`.
pub fn spectral_function<F: Fn(C) -> C>(f: F, z: C) -> C {
    spectral_function_with_eps(f, z, 1E-12)
}
