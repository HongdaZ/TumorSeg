predictLabel <- function( old, new, burnin = 1000, n = 10000 ) {
    s <- length( old ) # # of patients in training set
    alpha <- matrix( 1, nrow = 4, ncol = s )
    beta <- matrix( 1, nrow = 3, ncol = s )
    gamma <- rep( 1, s )
    alpha_star <- rep( 1, 4 )
    beta_star <- rep( 1, 3 )
    gamma_star <- 1
    lambda_alpha <- rep( 1, 4 )
    lambda_beta <- rep( 1, 3 )
    lambda_gamma <- 1
    mu <- vector( "list", 4 )
    sigma <- vector( "list", 4 )
    initial <- sumPatient( old )
    for( i in 1 : 4 ) {
        mu[[ i ]] <- initial$sum_y[[ i ]] / initial$n_type[ i ]
        sigma[[ i ]] <- 1 / ( initial$n_type[ i ] - 1 ) *
            ( initial$sum_cross_y[[ i ]] - initial$n_type[ i ] *
                  mu[[ i ]] %*%
                  t( mu[[ i ]] ) )
    }
    # Gamma( a, b )
    a <- b <- .001
    m <- 4
    tau_0 <- .001

    for( i in 1 : ( burnin + n ) ) {
        for( j in 1 : s ) {
            # update alpha's
            for( k in 1 : 4 ) {
                prop_alpha <- alpha[ k, j ] + rnorm( 1, mean = 0, sd = .1 )

            }
        }
    }

}
