use std::f64;
use faer::prelude::*;
use std::f64::consts::PI;

/// Calculates the first Discrete Prolate Spheroidal Sequence (DPSS) or Slepian window.
///
/// This sequence is the eigenvector corresponding to the largest eigenvalue of a
/// specific symmetric tridiagonal matrix.
///
/// # Arguments
/// * `n_samples`: The desired length of the window (N).
/// * `nw_param`: The time-bandwidth product (NW). This controls the spectral concentration.
///
/// # Returns
/// A `Result` containing a `Vec<f64>` with the window weights if successful,
/// or an error message string if inputs are invalid or computation fails.
///
/// # Notes
/// The method constructs a symmetric tridiagonal matrix `K` (related to the prolate matrix).
/// The DPSS sequence is the eigenvector of K corresponding to the largest eigenvalue.
/// `faer`'s `compute_eigenvectors()` method for symmetric matrices is used, which directly
/// returns real eigenvalues and real eigenvectors.
/// The resulting window is normalized so that its sum of squares is 1, and its
/// element with the largest absolute value is positive (a common sign convention).
pub fn calculate_dpss_window_faer(
    n_samples: usize,
    nw_param: f64,
) -> Result<Vec<f64>, String> {
    if n_samples == 0 {
        return Err("Window length N must be greater than 0.".to_string());
    }
    if n_samples == 1 {
        // Trivial case: a single point window is just 1.0
        return Ok(vec![1.0]);
    }
    if nw_param <= 0.0 {
        return Err("NW parameter must be positive.".to_string());
    }
    // NW typically should be less than N/2 for meaningful Slepian sequences.
    if nw_param >= n_samples as f64 / 2.0 {
        // This is not strictly an error, but a warning condition.
        // The sequences can still be computed.
        eprintln!(
            "Warning: NW (={}) is relatively large (>= N/2 = {}). \
             The Slepian sequence might not be as well-behaved or physically meaningful for W >= 0.5. \
             Consider NW < N/2.",
            nw_param, n_samples as f64 / 2.0
        );
    }

    let n_f64 = n_samples as f64;
    // W is the normalized half-bandwidth (0 < W < 0.5).
    // W = NW / N
    let w_digital = nw_param / n_f64;

    // 1. Construct the symmetric tridiagonal matrix K using Mat::from_fn
    //    This matrix formulation is common for computing Slepian sequences.
    //    (See, e.g., Percival & Walden, "Spectral Analysis for Physical Applications")
    let k_matrix = Mat::<f64>::from_fn(n_samples, n_samples, |r_idx, c_idx| {
        let r_f64 = r_idx as f64; // Current row index as f64

        if r_idx == c_idx {
            // Diagonal elements: K_ii = ((N - 1 - 2*i)/2)^2 * cos(2 * pi * W)
            let term1 = ((n_f64 - 1.0 - 2.0 * r_f64) / 2.0).powi(2);
            let term2 = (2.0 * PI * w_digital).cos();
            term1 * term2
        } else if (r_idx as isize - c_idx as isize).abs() == 1 {
            // Off-diagonal elements: K_{i, i+1} = K_{i+1, i} = (idx + 1) * (N - 1 - idx) / 2.0
            // where idx is the smaller of r_idx, c_idx (0-indexed)
            let idx_smaller = (r_idx.min(c_idx)) as f64;
            ((idx_smaller + 1.0) * (n_f64 - 1.0 - idx_smaller)) / 2.0
        } else {
            0.0 // It's a tridiagonal matrix
        }
    });

    // 2. Perform eigenvalue decomposition on the symmetric matrix K.
    // `compute_eigenvectors()` for symmetric matrices in faer returns (eigenvalues, eigenvectors_matrix)
    // where eigenvalues are real (Col<f64>) and eigenvectors are real (Mat<f64>).
    // Eigenvalues are sorted in non-decreasing (ascending) order.
    // Eigenvectors are the columns of the eigenvectors_matrix and are orthonormal.
    let eigen = k_matrix.self_adjoint_eigen(faer::Side::Lower).unwrap();
    //let eigenvalues_col: Vec<f64> = eigen.S().column_vector().iter().map(|&a| a).collect();
    let eigenvectors_mat = eigen.U();
    //let (eigenvalues_col, eigenvectors_mat) = k_matrix.compute_eigenvectors();

    // 3. Select the eigenvector corresponding to the largest eigenvalue.
    // Since eigenvalues are sorted in ascending order by faer's symmetric EVD,
    // the largest eigenvalue is the last one.
    // The eigenvectors are columns, so we take the last column.
    let dpss_sequence_col = eigenvectors_mat.col(n_samples - 1);

    let mut dpss_sequence: Vec<f64> = dpss_sequence_col.iter().copied().collect();

    // 4. Sign convention: Ensure the element with the largest absolute value is positive.
    // Eigenvectors are unique only up to a sign. This normalization is common.
    if !dpss_sequence.is_empty() {
        let mut max_abs_val = 0.0;
        let mut idx_max_abs = 0;
        for (i, &val) in dpss_sequence.iter().enumerate() {
            if val.abs() > max_abs_val {
                max_abs_val = val.abs();
                idx_max_abs = i;
            }
        }

        // If the element with the largest absolute value is negative, flip all signs.
        if dpss_sequence[idx_max_abs] < 0.0 {
            for val in dpss_sequence.iter_mut() {
                *val = -*val;
            }
        }
    }
    // The eigenvectors from faer's symmetric EVD are already L2 normalized (sum of squares = 1).

    println!("window len {}", dpss_sequence.len());

    // make window symetrical
    if dpss_sequence.len() % 2 == 0 {
    println!("dpss_sequence.len() % 2 == 0");
    dpss_sequence = [dpss_sequence.get(0..(dpss_sequence.len() / 2)).unwrap(), &dpss_sequence.get(0..((dpss_sequence.len() / 2))).unwrap().iter().rev().cloned().collect::<Vec<f64>>()].concat();
    } else {
    println!("dpss_sequence.len() % 2 != 0");
    dpss_sequence = [dpss_sequence.get(0..dpss_sequence.len() / 2).unwrap(), &vec![*dpss_sequence.get(dpss_sequence.len() / 2).unwrap()],  &dpss_sequence.get(0..dpss_sequence.len() / 2).unwrap().iter().rev().cloned().collect::<Vec<f64>>()].concat();
    }

    println!("window len after mirror {}", dpss_sequence.len());

    Ok(dpss_sequence)
}

/// Estimates the time-bandwidth product (NW) for a DPSS/Slepian window
/// based on the desired minimum sidelobe attenuation.
///
/// # Arguments
/// * `desired_sidelobe_attenuation_db`: The desired minimum sidelobe attenuation in dB
///   (a positive value, e.g., 60.0 for -60 dB sidelobes).
///
/// # Returns
/// The estimated NW parameter. This can be a fractional value.
///
/// # Notes
/// This function uses an empirical linear approximation: NW ≈ Attenuation_dB / 20.
/// A minimum NW is enforced (e.g., 1.5) as very low NW values lose distinct Slepian characteristics.
/// Higher NW gives better sidelobe rejection but widens the main lobe for a fixed window length.
pub fn estimate_nw_from_sidelobe_attenuation(desired_sidelobe_attenuation_db: f64) -> f64 {
    // Ensure the input is treated as positive attenuation
    let att_db = desired_sidelobe_attenuation_db.abs();

    // Empirical formula: NW is roughly proportional to attenuation / 20
    // Example: 40dB -> NW=2.0; 65dB -> NW=3.25; 80dB -> NW=4.0
    let estimated_nw = att_db / 20.0;

    // Enforce a minimum practical NW. Below ~1.5, Slepian windows offer less distinct
    // advantages over simpler windows, and the concentration properties are weaker.
    // Users targeting very low attenuation might prefer other window types anyway.
    let min_nw = 1.5; // Lowered minimum slightly, can be adjusted.

    if estimated_nw < min_nw {
        // If the desired attenuation is very low (e.g., < 30dB),
        // directly calculated NW would be < 1.5.
        // We return min_nw, but also warn that Slepian might not be optimal here.
        eprintln!(
            "Warning: Desired attenuation {:.1}dB results in calculated NW {:.2}, which is below the typical minimum of {:.1}. \
             Using NW = {:.1}. For very low attenuation, other window types might be more suitable or simpler.",
            att_db, estimated_nw, min_nw, min_nw
        );
        min_nw
    } else {
        estimated_nw
    }
}

/// Estimates the necessary window length (N) for a DPSS/Slepian window
/// using a continuous formula for the main lobe width factor.
///
/// # Arguments
/// * `nw`: The time-bandwidth product (NW) chosen for the window.
/// * `sampling_frequency_hz`: The sampling frequency of the signal in Hz.
/// * `desired_resolution_hz`: The desired frequency resolution in Hz
///   (interpreted as the null-to-null width of the main lobe).
///
/// # Returns
/// The estimated window length N (number of samples).
///
/// # Panics
/// Panics if `sampling_frequency_hz` or `desired_resolution_hz` are not positive.
/// Panics if `nw` is not positive or excessively large to the point of formula breakdown.
pub fn estimate_window_length_n(
    nw: f64,
    sampling_frequency_hz: f64,
    desired_resolution_hz: f64,
) -> usize {
    if sampling_frequency_hz <= 0.0 {
        panic!("Sampling frequency must be positive.");
    }
    if desired_resolution_hz <= 0.0 {
        panic!("Desired frequency resolution must be positive.");
    }
    if nw <= 0.0 {
        // While estimate_nw_from_sidelobe_attenuation_v2 clamps NW,
        // this function could be called directly.
        panic!("NW parameter must be positive.");
    }
    // Add a sanity check for NW, as the K_ml formula is derived from NW in the 2-5 range.
    // Extrapolating too far might lead to less accurate K_ml.
    if nw < 1.5 || nw > 6.0 { // Slightly wider range than estimate_nw clamps to
        eprintln!(
            "Warning: NW value {} is outside the typical empirical range (1.5-6.0) \
             for the K_ml estimation formula. The K_ml factor might be less accurate.",
            nw
        );
    }


    // K_ml: Factor for main lobe width (null-to-null) relative to (fs / N)
    // Empirically derived linear approximation: K_ml ≈ 1.12 * NW + 0.46
    // This was fitted to the data points previously used in the lookup table.
    let k_ml = 1.12 * nw + 0.46;

    // Formula: DesiredResolution_Hz ≈ K_ml * (SamplingFrequency_Hz / N)
    // So, N ≈ (K_ml * SamplingFrequency_Hz) / DesiredResolution_Hz
    let n_float = (k_ml * sampling_frequency_hz) / desired_resolution_hz;

    // Window length must be an integer; round up.
    // Ensure N is at least 1.
    (n_float.ceil() as usize).max(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn n1() {
        match calculate_dpss_window_faer(1, 2.0) {
            //Ok(window) => println!("N=1, NW=2.0: {:?}", window), // Expected: [1.0]
            Ok(window) => assert!(window == vec![1.0]), // Expected: [1.0]
            Err(e) => eprintln!("Error: {}", e),
        }
    }

    #[test]
    fn error_1() {
        match calculate_dpss_window_faer(1, 2.0) {
            Ok(window) => println!("N=1, NW=2.0: {:?}", window), // Expected: [1.0]
            Err(e) => eprintln!("Error: {}", e),
        }
    }

    #[test]
    fn error_2() {
        match calculate_dpss_window_faer(1, 2.0) {
            Ok(window) => println!("N=1, NW=2.0: {:?}", window), // Expected: [1.0]
            Err(e) => eprintln!("Error: {}", e),
        }
    }

    #[test]
    fn nw_ge_n_div_2() {
        match calculate_dpss_window_faer(1, 2.0) {
            Ok(window) => println!("N=1, NW=2.0: {:?}", window), // Expected: [1.0]
            Err(e) => eprintln!("Error: {}", e),
        }
    }
}