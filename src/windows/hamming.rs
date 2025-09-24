use std::{f64::consts::TAU, ops::Neg};

pub fn hamming(n: f64, window_size: f64) -> f64 {

    // const ALPHA: f64 = 0.53836;
    // const BETA: f64 = 0.46164;
    // const PI2: f64 = (PI * 2.0);

    // ALPHA - BETA * ((PI2 * n) / (f64::EPSILON - 1.0)).cos()

    0.5 - (0.5 * ((TAU * n) / window_size).cos())


}

pub fn hamming2(n: f64) -> f64 {
    
    const WINDOW_SIZE: f64 = 1024.0;
    // const ALPHA: f64 = 0.53836;
    // const BETA: f64 = 0.46164;
    // const PI2: f64 = (PI * 2.0);

    // ALPHA - BETA * ((PI2 * n) / (f64::EPSILON - 1.0)).cos()

    (0.5 - (0.5 * (TAU * (n / 2.0) / (WINDOW_SIZE - 1.0)).cos())).neg() + 1.0


}