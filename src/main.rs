// Constants and conversion ratios
static KBOLTZ: f64  = 1.380658e-16;
static AMU: f64     = 1.660539e-24;
static EVTOERG: f64 = 1.602000e-12;
static MASSH: f64   = 1.00784 * AMU;
static MASSE: f64   = 0.00054858 * AMU;
static HERATIO: f64 = MASSE/MASSH;

fn calc_h_coll(t: f64, a: f64, n_h: f64)->f64 {
    let kt     = KBOLTZ * t;
    let eh_kev = 133.0 * a;
    let eh     = 1.0e3 * EVTOERG * eh_kev;
    let h_n    = 1.0 - (1.0 + eh / (2.0 * kt)) * (-eh/kt).exp();
    // Determine collisional energy loss for colliding H atoms
    let h_coll = 1.26e-19 * a.powf(2.0) * t.powf(1.5) * n_h * h_n;
    h_coll
}

fn calc_h_e(x_e: f64)->f64 {
    let nzmax   = 401;  // Number of integraton bins
    let logzmin = -2.0; // Minimum log(z) value to consider
    let logzmax = 2.0;  // Maximum log(z) value to consider
    // Get log spacing
    let dlogz   = (logzmax - logzmin) / ((nzmax - 1) as f64);
    // Initialise arrays
    let mut f     = vec![0.0;nzmax];
    let mut z     = vec![0.0;nzmax];
    let mut expmz = vec![0.0;nzmax];
    let mut dz    = vec![0.0;nzmax - 1];
    // Create z bins
    for n in 0..nzmax {
        let fln  = n as f64;
        z[n]     = f64::powf(10.0,logzmin + fln*dlogz);
        expmz[n] = (-z[n]).exp();
    }
    for n in 0..nzmax-1 {
        dz[n] = z[n+1] - z[n];
    }
    // Evaluate the function at each z position
    // Integrate using the trapezium rule
    let x_ep1p5 = x_e.powf(1.5);
    for n in 0..nzmax {
        let zpxe = z[n] + x_e;
        let f1   = zpxe.powf(1.5) - x_ep1p5;
        f[n]     = zpxe*f1.powf(2.0/3.0) * expmz[n];
    }
    let mut intf = 0.0;
    for n in 0..nzmax-1 {
        intf += 0.5 * (f[n+1] + f[n]) * dz[n];
    }
    let ix_star = 0.5 * (-x_e).exp() * intf;
    let h_e     = 1.0 - ix_star;
    // Return h_e
    h_e
}

fn calc_h_el(t:f64, a:f64, n_e:f64)->f64 {
    let kt      = KBOLTZ * t;
    let e_e_kev = 23.0 * a.powf(2.0/3.0);
    let e_e     = e_e_kev * 1.0e3 * EVTOERG;
    let x_e     = e_e / kt;
    let h_e     = calc_h_e(x_e);
    let h_el    = (1.26e-19 * a.powf(2.0) * t.powf(1.5) * n_e * h_e) / (HERATIO).sqrt();
    // Return H_el
    h_el
}

fn main() {
    // Gas parameters
    let rho = 1e-20;  // Gas density (g/cm^3)
    let a   = 5e-3;   // Grain radius (micron)
    // Calculate number density, currently only solar abundance but can fix
    let n_h = rho * (10.0/14.0) / MASSH;
    let n_e = 1.2 * n_h;
    // Build temperature range array
    let logtmin = 0.0;
    let logtmax = 9.0;
    let ntmax   = 201;
    let dlogt   = (logtmax - logtmin) / ((ntmax - 1) as f64);
    // Initialise temperature range array
    let mut t_arr = vec![0.0;ntmax];
    // Calculate T
    for n in 0..ntmax {
        let nf    = n as f64;
        let log_t = logtmin + (nf * dlogt);
        t_arr[n]  = f64::powf(10.0,log_t);
    }
    // Calculate lambda
    let mut lambda_arr = Vec::new();
    // Process temperature bins
    for t in t_arr {
        let h_coll = calc_h_coll(t,a,n_h);
        let h_el   = calc_h_el(t,a,n_e);
        let lambda = (h_coll + h_el)/ n_h;
        // Add to array containing values for lambda
        lambda_arr.push(lambda);
    }
}