use plotly::{common::Mode, Layout, Plot, Scatter};
use iter_num_tools::lin_space;
use rustfft::FftPlannerAvx;
use rustfft::num_complex::Complex;

/// Computes the full symmetrical frequency response of a window.
///
/// # Arguments
/// * `window_weights`: The weights of the window.
/// * `fft_size`: The number of points for the FFT.
///
/// # Returns
/// A tuple `(frequencies, magnitudes_db)`:
/// * `frequencies`: Vector of normalized frequencies (-0.5 to ~0.5).
/// * `magnitudes_db`: Vector of corresponding magnitudes in dB.
pub fn compute_window_frequency_response_symmetric(
    window_weights: &[f64],
    fft_size: usize,
) -> (Vec<f64>, Vec<f64>) {
    if window_weights.is_empty() {
        return (Vec::new(), Vec::new());
    }
    if fft_size < window_weights.len() {
        panic!("FFT size must be greater than or equal to window length.");
    }

    let mut padded_window = vec![0.0; fft_size];
    padded_window[..window_weights.len()].copy_from_slice(window_weights);

    // --- FFT Calculation (using rustfft is recommended for performance) ---
    // Placeholder: Direct DFT - REPLACE WITH RUSTFFT FOR EFFICIENCY
    // let mut complex_fft_result = vec![(0.0, 0.0); fft_size]; // (real, imag)
    // for k in 0..fft_size {
    //     let mut real_sum = 0.0;
    //     let mut imag_sum = 0.0;
    //     for n in 0..fft_size {
    //         let angle = -2.0 * PI * (k as f64) * (n as f64) / (fft_size as f64);
    //         real_sum += padded_window[n] * angle.cos();
    //         imag_sum += padded_window[n] * angle.sin();
    //     }
    //     complex_fft_result[k] = (real_sum, imag_sum);
    // }

    let mut planner = FftPlannerAvx::new().unwrap();
    let fft = planner.plan_fft_forward(fft_size);
    let mut buffer: Vec<Complex<f64>> = padded_window.into_iter().map(|x| Complex::new(x, 0.0)).collect();
    fft.process(&mut buffer);
    let complex_fft_result_rustfft: Vec<(f64, f64)> = buffer.into_iter().map(|c| (c.re, c.im)).collect();
    let complex_fft_result = complex_fft_result_rustfft;

    //let magnitudes_raw: Vec<f64> = buffer.iter().map(|a| a.abs()).collect();

    // Calculate magnitudes
    let magnitudes_raw: Vec<f64> = complex_fft_result
        .iter()
        .map(|(re, im)| (re.powi(2) + im.powi(2)).sqrt())
        .collect();

    // phase
    // let magnitudes_raw: Vec<f64> = complex_fft_result
    //     .iter()
    //     .map(|(re, im)| im.atan2(*re))
    //     .collect();

    // Normalize magnitudes (so peak is 0 dB) and convert to dB
    let max_mag = magnitudes_raw.iter().cloned().fold(0.0/0.0, f64::max);
    let mut magnitudes_db_full = Vec::with_capacity(fft_size);
    if max_mag > 1e-9 {
        for mag in &magnitudes_raw {
            let val_db = 10.0 * (mag / max_mag).log10();
            magnitudes_db_full.push(val_db.max(-144.0));
        }
    } else {
        magnitudes_db_full.resize(fft_size, -144.0);
    }

    // Create "fftshifted" version for plotting frequencies from -0.5 to 0.5
    let mut magnitudes_db_shifted = vec![0.0; fft_size];
    let mut frequencies_shifted = vec![0.0; fft_size];

    let midpoint = fft_size / 2;
    for i in 0..fft_size {
        if i < midpoint {
            // Second half of FFT output maps to negative frequencies
            magnitudes_db_shifted[i] = magnitudes_db_full[midpoint + i];
            frequencies_shifted[i] = (midpoint + i) as f64 / fft_size as f64 - 1.0;
        } else {
            // First half of FFT output maps to positive frequencies
            magnitudes_db_shifted[i] = magnitudes_db_full[i - midpoint];
            frequencies_shifted[i] = (i - midpoint) as f64 / fft_size as f64;
        }
    }
    // Ensure the very last point (Nyquist if fft_size is even) is correctly handled if needed,
    // but for plotting the symmetry, this shift is generally fine.
    // For an even fft_size, freq[fft_size/2] is Nyquist, and its symmetric counterpart is itself.
    // The shift above effectively puts DC at `midpoint` in the shifted array.

    (frequencies_shifted, magnitudes_db_shifted)
}


pub fn draw_f64(s2: &Vec<f64>) {

    // let min = s2.iter().cloned().reduce(f64::min).unwrap();
    // let max = s2.iter().cloned().reduce(f64::max).unwrap();

    let x: Vec<f64> = lin_space(-1.0..=1.0, s2.len()).collect();

    let trace = Scatter::new(x.clone(), s2.clone())
    .mode(Mode::Lines);

    let mut plot = Plot::new();
    let layout = Layout::new().height(1000).width(1000).auto_size(false);
    plot.set_layout(layout);

    plot.add_trace(trace);
    plot.show();
}

pub fn draw_f64_twoax(s: Vec<f64>, s2: Vec<f64>) {

    let trace = Scatter::new(s, s2)
    .mode(Mode::Lines);

    let mut plot = Plot::new();
    let layout = Layout::new().height(1000).width(1000).auto_size(false);
    plot.set_layout(layout);

    plot.add_trace(trace);
    plot.show();
}   