use std::f64::consts::PI;

pub fn gaussian_kernel(size: usize, sigma: f64) -> Vec<f64> {
    let mut kernel = Vec::with_capacity(size);
    let center = (size as f64) / 2.0;

    for i in 0..size {
        let x = i as f64 - center;
        let weight = (1.0 / (sigma * (2.0 * PI).sqrt())) * (-x.powi(2) / (2.0 * sigma.powi(2))).exp();
        kernel.push(weight);
    }

    kernel
}

pub fn gaussian(x: f64) -> f64 {
    // Set fixed values for mu and sigma
    let mu = 0.0; // Example: Mean = 0
    let sigma = 1.0; // Example: Standard deviation = 1

    let coefficient = 1.0 / (sigma * (2.0_f64).sqrt());
    let exponent = -((x - mu).powi(2)) / (2.0 * sigma.powi(2));
    coefficient * exponent.exp()
}