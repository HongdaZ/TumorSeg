pmfLabel <- function( alpha_star, beta_star, gamma_star, neighbor_label,
                      det_sigma, inv_sigma, sigma_mu, q_sigma, pw, hoi ) {
    .Call( "pmfLabel", alpha_star, beta_star, gamma_star, neighbor_label,
           det_sigma, inv_sigma, sigma_mu, q_sigma, pw, hoi )
}
