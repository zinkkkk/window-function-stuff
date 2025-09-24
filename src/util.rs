
use std::f64::consts::{PI, TAU};

use num_traits::Float;

pub fn apply_window_inplace<F: Float>(signal: &mut [F], window_function: &[F]) {
    signal.iter_mut().zip(window_function.iter().cloned()).for_each(|(s, w)| *s = s.clone() * w)
}

pub fn apply_window<F: Float>(signal: &[F], window_function: &[F]) -> Vec<F> {
    let mut processed_signal = Vec::with_capacity(signal.len());
    signal.iter().zip(window_function.iter()).for_each(|(&s, &w)| processed_signal.push(s * w));
    processed_signal
}


// pub fn apply_window<F: Float>(signal: &[F], window_function: &[F]) -> Vec<F> {
//     let mut processed_signal = Vec::with_capacity(signal.len());
//     signal.iter().zip(window_function.iter()).for_each(|(&s, &w)| processed_signal.push(s * w));
//     processed_signal
// }

// pub fn apply_window_inplace<N: Numeric>(signal: &mut [N], window_function: &[N])
// {
//     signal.iter_mut().zip(window_function.iter().cloned()).for_each(|(s, w)| *s = s.clone() * w)
// }


// pub fn scale_to_zero<F>(data: &mut [F])
// where F: Float + core::iter::Sum + NumAssign
// {
//     let iterations: usize = 5000;
//     for _ in 0..iterations {
//         let sum: F = data.iter().cloned().sum();
//         let count = F::from(data.len()).unwrap();

//         if sum != F::zero() {
//             let scale_factor = F::one() - (sum / (count * data.iter().map(|x| x.abs()).sum()));
//             data.iter_mut().for_each(|x| *x *= scale_factor);
//         } else {
//             break;
//         }
//     }
// }

pub fn normalize_weights(weights: &mut Vec<f64>) {
    let sum: f64 = weights.iter().sum::<f64>() / 0.5;

    if sum != 0.0 {
        for weight in weights.iter_mut() {
            *weight /= sum;
        }
    }
}

pub fn enbw(v: Vec<f64>) -> f64 {

    let len = v.len() as f64;

    let numer = v.iter().map(|a| a.powi(2).abs()).sum::<f64>();
    let denom = v.iter().sum::<f64>().powi(2).abs();

    len * (numer / denom)

}

pub fn sum_of_squares(dpss_vec: &[f64]) -> f64 {
    dpss_vec.iter().map(|a| a * a).sum()
}

// old 
pub fn sine_wave(frequency: f64, sample_rate: f64, duration_seconds: f64) -> Vec<f64> {
    let num_samples = (sample_rate * duration_seconds) as usize;
    let mut sine_wave = Vec::with_capacity(num_samples);
    for i in 0..num_samples {
        let t = i as f64 / sample_rate;
        sine_wave.push((2.0 * PI * frequency * t).sin());
    }
    sine_wave
}


// pub fn sine_wave(frequency: f64, sample_rate: i64, duration_seconds: f64) -> Vec<f64> {
//     let num_samples = sample_rate * duration_seconds as i64;
//     //let mut sine_wave: Vec<f64> = Vec::with_capacity(num_samples);
//     // for i in 0..num_samples {
//     //     let t = i as f64 / sample_rate;
//     //     sine_wave.push((2.0 * PI * frequency * t).sin());
//     // }
//     // sine_wave

//     let sine_wave = (0..num_samples).into_iter().map(|i|{
//         let t = i / sample_rate;
//         (2.0 * PI * frequency * t as f64).sin()
//     }).collect();

//     sine_wave
// }


// Linear frequency sweep
pub fn sine_sweep(sample_rate: u32, duration_seconds: f64) -> Vec<f64> {
    let num_samples = (sample_rate as f64 * duration_seconds).round() as usize;
    let mut samples = Vec::with_capacity(num_samples);

    for i in 0..num_samples {
        let time = i as f64 / sample_rate as f64;
        let frequency = time * 2_000.0 / duration_seconds; 
        let phase = TAU * frequency * time;
        samples.push(phase.sin());
    }

    samples
}

// Logarithmic frequency sweep
pub fn sine_sweep_logarithmic(sample_rate: u32, duration_seconds: f64) -> Vec<f64> {
    let num_samples = (sample_rate as f64 * duration_seconds).round() as usize;
    let mut samples = Vec::with_capacity(num_samples);

    for i in 0..num_samples {
        let time = i as f64 / sample_rate as f64;
        let frequency = 2_000.0 * (10.0f64).powf(time / duration_seconds); 
        let phase = TAU * frequency * time;
        samples.push(phase.sin());
    }

    samples
}