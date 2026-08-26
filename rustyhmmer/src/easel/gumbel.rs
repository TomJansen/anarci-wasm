



use super::constants::ESL_SMALLX1;













#[inline]
pub fn esl_gumbel_pdf(x: f64, mu: f64, lambda: f64) -> f64 {
    let y = lambda * (x - mu);
    lambda * (-y - (-y).exp()).exp()
}













#[inline]
pub fn esl_gumbel_logpdf(x: f64, mu: f64, lambda: f64) -> f64 {
    let y = lambda * (x - mu);
    lambda.ln() - y - (-y).exp()
}













#[inline]
pub fn esl_gumbel_cdf(x: f64, mu: f64, lambda: f64) -> f64 {
    let y = lambda * (x - mu);
    (-(-y).exp()).exp()
}













#[inline]
pub fn esl_gumbel_logcdf(x: f64, mu: f64, lambda: f64) -> f64 {
    let y = lambda * (x - mu);
    -(-y).exp()
}















#[inline]
pub fn esl_gumbel_surv(x: f64, mu: f64, lambda: f64) -> f64 {
    let y = lambda * (x - mu);
    let ey = -(-y).exp();

    
    if ey.abs() < ESL_SMALLX1 {
        -ey
    } else {
        1.0 - ey.exp()
    }
}

















#[inline]
pub fn esl_gumbel_logsurv(x: f64, mu: f64, lambda: f64) -> f64 {
    let y = lambda * (x - mu);
    let ey = -(-y).exp();

    
    
    
    
    if ey.abs() < ESL_SMALLX1 {
        -y
    } else if ey.exp().abs() < ESL_SMALLX1 {
        -ey.exp()
    } else {
        (1.0 - ey.exp()).ln()
    }
}












#[inline]
pub fn esl_gumbel_invcdf(p: f64, mu: f64, lambda: f64) -> f64 {
    mu - ((-p.ln()).ln() / lambda)
}












#[inline]
pub fn esl_gumbel_invsurv(p: f64, mu: f64, lambda: f64) -> f64 {
    
    
    
    let log_part = if p < ESL_SMALLX1 {
        (p.powf(p) - 1.0) / p
    } else {
        (-(1.0 - p).ln()).ln()
    };

    mu - (log_part / lambda)
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPSILON: f64 = 1e-10;

    fn approx_eq(a: f64, b: f64) -> bool {
        if a.is_infinite() && b.is_infinite() {
            return a.signum() == b.signum();
        }
        if a.abs() < EPSILON && b.abs() < EPSILON {
            return true;
        }
        (a - b).abs() / a.abs().max(b.abs()).max(1.0) < EPSILON
    }

    #[test]
    fn test_gumbel_pdf_at_mode() {
        
        
        let pdf = esl_gumbel_pdf(0.0, 0.0, 1.0);
        let expected = 1.0 / std::f64::consts::E;
        assert!(
            approx_eq(pdf, expected),
            "PDF at mode: {} vs {}",
            pdf,
            expected
        );
    }

    #[test]
    fn test_gumbel_cdf_properties() {
        
        let cdf = esl_gumbel_cdf(0.0, 0.0, 1.0);
        let expected = (-1.0_f64).exp();
        assert!(
            approx_eq(cdf, expected),
            "CDF at mu: {} vs {}",
            cdf,
            expected
        );
    }

    #[test]
    fn test_gumbel_surv_cdf_relationship() {
        
        let x = 1.5;
        let mu = 0.0;
        let lambda = 1.0;
        let surv = esl_gumbel_surv(x, mu, lambda);
        let cdf = esl_gumbel_cdf(x, mu, lambda);
        assert!(
            approx_eq(surv + cdf, 1.0),
            "surv + cdf should equal 1: {} + {} = {}",
            surv,
            cdf,
            surv + cdf
        );
    }

    #[test]
    fn test_gumbel_logpdf_pdf_relationship() {
        
        let x = 0.5;
        let mu = 0.0;
        let lambda = 1.0;
        let pdf = esl_gumbel_pdf(x, mu, lambda);
        let logpdf = esl_gumbel_logpdf(x, mu, lambda);
        assert!(
            approx_eq(logpdf, pdf.ln()),
            "logpdf should equal ln(pdf): {} vs {}",
            logpdf,
            pdf.ln()
        );
    }

    #[test]
    fn test_gumbel_invcdf_cdf_inverse() {
        
        let x = 1.5;
        let mu = -20.0;
        let lambda = 0.4;
        let p = esl_gumbel_cdf(x, mu, lambda);
        let recovered = esl_gumbel_invcdf(p, mu, lambda);
        assert!(
            approx_eq(x, recovered),
            "invcdf(cdf(x)) should equal x: {} vs {}",
            x,
            recovered
        );
    }
}
