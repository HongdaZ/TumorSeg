pmfLabel <- function( alpha_star, beta_star, gamma_star,
                      det_sigma, inv_sigma, sigma_mu,
                      q_sigma, pw, hoi, neighbor_label, prob, p ) {
    .Call( "pmfLabel", alpha_star, beta_star, gamma_star,
           det_sigma, inv_sigma, sigma_mu,
           q_sigma, pw, hoi, neighbor_label, prob )
}
