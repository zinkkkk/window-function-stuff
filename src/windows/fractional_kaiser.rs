//use num_traits::Pow;

// this is a mess i cant remember what bits work and what bits are half failed/broken good luck...
// im pretty sure one/some of it works

// fn fib(n: f64) -> f64 {

//     let nn: u64 = n.round() as u64;

//     if nn == 0 {
//         return 0.0;
//     } else if nn == 1{
//         return 1.0;
//     }
//     fibu(nn) as f64
// }

// fn fibu(n: u64) -> u64 {


//     if n == 0 {
//         return 0;
//     } else if n == 1{
//         return 1;
//     }
//     fibu(n-1) + fibu(n-2)
// }

// fn fib89(n: f64) -> f64 {

//     (1.0 / n) / 89.0

// }

// pub fn factorial(n: u64) -> f64 {
//     if n == 0 {
//         1.0
//     } else {
//         (1..=n).product::<u64>() as f64
//     }
// }

// Factorial function that returns f64 to avoid overflow and for direct use
fn factorial(num: u128) -> f64 {
    if num == 0 {
        1.0 // 0! = 1
    } else {
        (1..=num).fold(1.0, |acc, v| acc * v as f64)
    }
}

// Factorial function that returns f64 to avoid overflow and for direct use
fn factorial_f64(num: f64) -> f64 {
    if num == 0.0 {
        1.0 // 0! = 1
    } else {
        (1..=num as u128).fold(1.0, |acc, v| acc * v as f64)
    }
}

// hmmm?
// fn gamma_recursive(z: f64) -> f64 {
//     if z == 1.0 {
//         1.0
//     } else if z == 0.5 {
//         PI.sqrt()
//     } else if z > 1.0 {
//         (z - 1.0) * gamma_recursive(z - 1.0)
//     } else if z > 0.0 {

//         // Series expansion for 0 < z < 1 (Euler's Reflection Formula)
//         // Note: This is computationally expensive and may have convergence issues.

//         // hmm?
//         // let mut sum = 1.0 / z;
//         // let mut term = 1.0;
//         // for n in 1..100 { // Adjust the number of iterations as needed
//         //     term *= -z / (n as f64);
//         //     sum += term / (z + n as f64);
//         // }

//         // Euler's Reflection formula Pi / (sin(pi*z) * gamma(z) = gamma(1-z)
//         PI / (z * PI).sin() / gamma_recursive(1.0 - z)

//     } else if z < 0.0 && z == (z as i64) as f64 {
//       panic!("Gamma function is undefined for negative integers!");
//     } else if z < 0.0 {
//       PI / ((z * PI).sin() * gamma_recursive(1.0-z))
//     } else {
//         f64::NAN  // Handle invalid input
//     }
// }

// fn gammaa(z: f64) -> f64 {
//     gamma_recursive(z)
// }

// fn gamma2(mut z: f64) -> f64 {
//     if z == 1.0 {
//         1.0
//     } else if z == 0.5 {
//         PI.sqrt()
//     } else if z > 1.0 {
//         let mut acc = 1.0;
//         while z > 1.0 {
//             z -= 1.0;
//             acc *= z;
//         }
//         acc * if z == 0.5 {
//            PI.sqrt()
//         } else {
//             gamma_small(z)
//         }

//     } else if z > 0.0 {
//         gamma_small(z)
//     } else if z < 0.0 && z == (z as i64) as f64 {
//       panic!("Gamma function is undefined for negative integers!");
//     } else if z < 0.0 {
//       PI / ((z * PI).sin() * gamma2(1.0-z))
//     } else {
//         f64::NAN  // Handle invalid input
//     }
// }

fn gamma_small(z: f64) -> f64 {
      // Series expansion for 0 < z < 1 (Euler's Reflection Formula)
        let mut sum = 1.0 / z;
        let mut term = 1.0;
        for n in 1..64 { // Adjust the number of iterations as needed
            term *= -z / (n as f64);
            sum += term / (z + n as f64);
        }
        sum
}

pub fn besselja(x: f64) -> f64 {
    
    const ALPHA: u128 = 0;
    const N: u128 = 32;

    let mut result = 0.0;

    for m in 0..N {

    let num = -1.0_f64.powf(m as f64);
    let den = factorial(m) as f64 * ((m + ALPHA + 1) as f64).gamma();

    let div = num as f64 / den;

    let x2pow2ma = (x / 2.0).powf(((2 * m) + ALPHA) as f64);

    let mul = div * x2pow2ma;

    result += mul;

    }

    result
}

pub fn modified_bessel(z: f64) -> f64 {
    let v = 1.0;
    const N: usize = 51;
    let mut result = z;

    for i in 0..N {
    let k = i as f64;
    let eg = gamma_small(k + 1.0) * gamma_small(k + v + 1.0);
    let eg = eg.recip();
    let eg = eg * (z / 2.0).powf((2.0 * k) + v);
    result *= eg;
    }
    //println!("result {}", result);
    result
}

pub fn besselj0(x: f64) -> f64 {

    //const N_TERMS: u32 = 85; // Increased slightly for potentially better behavior at x=25
    let mut result = 0.0;
    let x_over_2 = x / 2.0;

    for m in 0..85 {
        let sign = if m % 2 == 0 { 1.0 } else { -1.0 };
        let fact_m = factorial_f64(m as f64);
        if fact_m.is_infinite() { break; }
        let denominator = fact_m * fact_m;
        if denominator == 0.0 {
            if x_over_2 == 0.0 && m > 0 { continue; }
            else if x_over_2 == 0.0 && m == 0 {} // Handled
            else { break; }
        }
        let power_term = x_over_2.powi((2 * m) as i32);
        let term = sign * power_term / denominator;
        result += term;
    }
    result
}

pub fn ia(x: f64) -> f64 {

    // const Y: f64 = 0.95;
    // const V: f64 = 0.0;

    const A: f64 = 1.0;
    //const Y: f64 = 0.7;

    // // Special cases
    // if x == 0.0 {
    //     if V == 0.0 {
    //         return 1.0;
    //     } else {
    //         return 0.0;
    //     }
    // }

    // Use series representation
    //let mut sum = 1.0;
    //let mut term = 1.0;

    let mut sum = 0.0;
    let mut term;

    let x_div = x / 2.0;

    // Calculate using power series
    for m in 0..20 {

        let m_f = m as f64;
        
        //let m_fact = factorial(m);

        let m_fact = (m_f + 1.0).gamma();
        let gam = (m_f + A + 1.0).gamma();

        term = (m_fact * gam).recip();
        let pow_t = x_div.powf((2.0 * m_f) + A);

        term = term * pow_t;

        //term = term * x * x / ((4.0 * Y) * (k_f + 1.0) * ((k_f * Y) + V + 1.0));
        sum += term;

        // if term.abs() < 1e-15 * sum.abs() {
        //     break;
        // }
    }

    // Apply leading factor
    //let factor = (x / 2.0).powf(Y + 1.0) / (1.0_f64 + Y).gamma();

    sum
}

// pub fn i0_f(x: f64) -> f64 {

//     const Y: f64 = 0.0;
//     const V: f64 = 0.95;



//     // Special cases
//     // if x == 0.0 {
//     //     if V == 0.0 {
//     //         return 1.0;
//     //     } else {
//     //         return 0.0;
//     //     }
//     // }

//     // Use series representation
//     let mut sum = 0.0;
//     let mut term = 0.0;

//     // Calculate using power series
//     for k in 0..20 {
//         let k_f = k as f64;
//         term = term * x * x / ((4.0 * Y) * (k_f + 1.0) * ((k_f * Y) + V + 1.0));
//         sum += term;

//         if term.abs() < 1e-15 * sum.abs() {
//             break;
//         }
//     }

//     // Apply leading factor
//     let factor = (x / 2.0).powf(Y + 1.0) / (1.0_f64 + Y).gamma();

//     factor * sum
// }

// pub fn fractional_modified_bessel(z: f64) -> f64 {


//     const Y: f64 = 0.9;
//     const V: f64 = 0.0;

//     const N_TERMS: usize = 20;

//     let z_over_2 = z / 2.0;

//     let mut result = 0.0;

//     for k in 0..N_TERMS {
//     //let sign = if k % 2 == 0 { -1.0 } else { 0.0 };
//     //let div = (factorial_f64(k as u32)).recip();

//     //let k_f = k as f64;

//     let power_term = z_over_2.powi((k as i32) * 2);

//     let gam_1 = ((k + 1) as f64).gamma();
//     let gam_2 = ((k as f64 * Y) + 1.0).gamma();

//     let denom = (gam_1 * gam_2);
//     let val = denom.recip();

//     let term = val * power_term;

//     result += term;
//     }
//     //println!("result {}", result);
//     result
// }

pub fn ia_f(x: f64) -> f64 {

    const Y: f64 = 0.7;
    const V: u128 = 1;

    let mut sum = 0.0;
    let mut term;

    let x_div = x / 2.0;

    for m in 0..30 {

        let m_f = m as f64;
        

        let numer = factorial(2 * m + V);

        let gam = ((Y * (2.0 * m_f + V as f64)) + 1.0).gamma();
        let m_fact1 = factorial(m);
        let m_fact2 = factorial(m + V);

        let denom = m_fact1 * m_fact2 * gam;

        term = numer / denom;

        let pow_t = x_div.powf((2.0 * m_f) + (V as f64));

        term = term * pow_t;

        sum += term;

    }

    sum
}

// might be wrong but i think this does work?
pub fn fractional_kaiser_window(n: f64) -> f64 {

    const WINDOW_SIZE: f64 = 1000.0;
    const ALPHA: f64 = 2.0;

    let numerator = ia_f(ALPHA * (1.0 - ((2.0 * n) / (WINDOW_SIZE - 1.0)).powi(2)).sqrt());
    let denominator = ia_f(ALPHA);
    numerator / denominator
}