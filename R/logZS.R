logZS <- function( new_patient, alpha_star, beta_star, gamma_star, p ) {
    .Call( "logZS", new_patient, pw, hoi, alpha_star, beta_star,
           gamma_star, p )
}
