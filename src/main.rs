use iter_num_tools::lin_space;
use window_fn::*;

fn main() {

    let mut wave = util::sine_wave(10.0, 100.0, 3.0);
    //let mut wave: Vec<_> = wave.iter().map(|&a| a as f64).collect();
    println!("wave len {}", wave.len());
    //let win = dpss::calculate_dpss_window_faer(wave.len(), 3.0).unwrap();

    let points = lin_space(-5.0..=5.0, 1000).collect::<Vec<f64>>();
    let win: Vec<f64> = points.iter().map(|a| fractional_kaiser::ia_f(*a)).collect();
    //let win = fractional_kaiser::besselj0(x);

    // let mut win = Vec::with_capacity(wave.len());
    // let pts: Vec<f64> = lin_space(-100.0..=100.0, wave.len()).collect();
    // pts.iter().for_each(|&a| win.push(windows::fractional_kaiser::fractional_kaiser_window(a)));
    //let win: Vec<f64> = win.iter().map(|&a| a as f128).collect();
    println!("win {:?}", win);
    //plotting::compute_window_frequency_response_symmetric(&win, win.len());
    //plotting::draw_f64(&wave);
    //plotting::draw_f64(&win);
    plotting::draw_f64_twoax(points, win.clone());

    util::apply_window_inplace(wave.as_mut_slice(), &win);
    //let b = wave;
    plotting::draw_f64(&wave);

}