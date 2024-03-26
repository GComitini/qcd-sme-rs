use crate::{Num, R};

pub fn find_root<T: Num, F: Fn(T) -> T>(
    f: F,
    guesses: (T, T),
    tol: R,
    max_iter: usize,
) -> Option<(T, T)> {
    let (mut z0, mut z1) = guesses;
    let (mut f0, mut f1) = (f(z0), f(z1));
    for _ in 0..max_iter {
        if f1.abs() < tol {
            return Some((z1, f1));
        }
        let z2 = (z0 * f1 - z1 * f0) / (f1 - f0);
        (f0, f1) = (f1, f(z2));
        (z0, z1) = (z1, z2);
    }
    None
}
