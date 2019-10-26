predictLabel <- function( old, new, burnin = 1000, n = 10000, core = 12 ) {
    s <- length( old ) # # of patients in training set
    alpha <- matrix( 1, nrow = 4, ncol = s )
    beta <- matrix( 1, nrow = 3, ncol = s )
    gamma <- rep( 1, s )
    alpha_star <- rep( 1, 4 )
    alpha_star_trace <- matrix( nrow = 4, ncol = burnin + n + 1 )
    alpha_star_trace[ , 1 ] <- alpha_star
    beta_star <- rep( 1, 3 )
    beta_star_trace <- matrix( nrow = 3, ncol = burnin + n + 1 )
    beta_star_trace[ , 1 ] <- beta_star
    gamma_star <- 1
    sd <- .1
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
                old_alpha <- alpha[ , j ]
                new_alpha <- old_alpha
                new_alpha[ k ] <- old_alpha[ k ] +
                    rnorm( 1, mean = 0, sd = sd )
                proportion <- exp( ( old_alpha[ k ] - new_alpha[ k ] ) *
                                       ( lambda_alpha[ k ] +
                                             old[[ j ]]$n_type[ k ] ) -
                                       ( logZ( old[[ j ]], old_alpha,
                                               beta[ , j ], gamma[ j ], core ) -
                                             logZ( old[[ j ]], new_alpha,
                                                   beta[ , j ],
                                                   gamma[ j ], core ) ) )
                rho <- min( 1, proportion )
                if( rbinom( 1, 1, rho ) == 1 ) {
                    alpha[ k, j ] <- new_alpha[ k ]
                }
            }
            # update beta's
            for( k in 1 : 3 ) {
                old_beta <- beta[ , j ]
                new_beta <- old_beta
                new_beta[ k ] <- old_beta[ k ] +
                    rnorm( 1, mean = 0, sd = sd )
                proportion <- exp( ( old_beta[ k ] - new_beta[ k ] ) *
                                       ( lambda_beta[ k ] +
                                             old[[ j ]]$beta_sum[ k ] ) -
                                       ( logZ( old[[ j ]], alpha[ , j ],
                                               old_beta, gamma[ j ], core ) -
                                             logZ( old[[ j ]], alpha[ , j ],
                                                   new_beta, gamma[ j ],
                                                   core ) ) )
                rho <- min( 1, proportion )
                if( rbinom( 1, 1, rho ) == 1 ) {
                    beta[ k, j ] <- new_beta[ k ]
                }
            }
            # update gamma
            old_gamma <- gamma[ j ]
            new_gamma <- old_gamma + rnorm( 1, 0, sd )
            proportion <- exp( ( old_gamma - new_gamma ) *
                                   ( lambda_gamma  + old[[ j ]]$gamma_sum ) -
                                   ( logZ( old[[ j ]], alpha[ , j ],
                                           beta[ , j ], old_gamma, core ) -
                                         logZ( old[[ j ]], alpha[ , j ],
                                               beta[ , j ], new_gamma, core ) ) )
            rho <- min( 1, proportion )
            if( rbinom( 1, 1, rho ) == 1 ) {
                gamma[ j ] <- new_gamma
            }
        }
        # get beta_sum and gamma_sum
        logZS( new, alpha_star, beta_star, gamma_star, new$beta_sum,
               new$gamma_sum, core )
        new_n_type <- vector( "numeric", 4 )
        for( k in 1 : 8 ) {
            new_n_type <- new_n_type + new[[ k ]]$n_type
        }
        # update star's
        # update alpha_star's
        for( k in 1 : 4 ) {
            old_alpha <- alpha_star
            new_alpha <- old_alpha
            new_alpha[ k ] <- old_alpha[ k ] +
                rnorm( 1, mean = 0, sd = sd )
            proportion <- exp( ( old_alpha[ k ] - new_alpha[ k ] ) *
                                   ( lambda_alpha[ k ] +
                                         old[[ j ]]$n_type[ k ] ) -
                                   ( logZ( old[[ j ]], old_alpha,
                                           beta[ , j ], gamma[ j ], core ) -
                                         logZ( old[[ j ]], new_alpha,
                                               beta[ , j ],
                                               gamma[ j ], core ) ) )
            rho <- min( 1, proportion )
            if( rbinom( 1, 1, rho ) == 1 ) {
                alpha[ k, j ] <- new_alpha[ k ]
            }
        }
        # update beta's
        for( k in 1 : 3 ) {
            old_beta <- beta[ , j ]
            new_beta <- old_beta
            new_beta[ k ] <- old_beta[ k ] +
                rnorm( 1, mean = 0, sd = sd )
            proportion <- exp( ( old_beta[ k ] - new_beta[ k ] ) *
                                   ( lambda_beta[ k ] +
                                         old[[ j ]]$beta_sum[ k ] ) -
                                   ( logZ( old[[ j ]], alpha[ , j ],
                                           old_beta, gamma[ j ], core ) -
                                         logZ( old[[ j ]], alpha[ , j ],
                                               new_beta, gamma[ j ],
                                               core ) ) )
            rho <- min( 1, proportion )
            if( rbinom( 1, 1, rho ) == 1 ) {
                beta[ k, j ] <- new_beta[ k ]
            }
        }
        # update gamma
        old_gamma <- gamma[ j ]
        new_gamma <- old_gamma + rnorm( 1, 0, sd )
        proportion <- exp( ( old_gamma - new_gamma ) *
                               ( lambda_gamma  + old[[ j ]]$gamma_sum ) -
                               ( logZ( old[[ j ]], alpha[ , j ],
                                       beta[ , j ], old_gamma, core ) -
                                     logZ( old[[ j ]], alpha[ , j ],
                                           beta[ , j ], new_gamma, core ) ) )
        rho <- min( 1, proportion )
        if( rbinom( 1, 1, rho ) == 1 ) {
            gamma[ j ] <- new_gamma
        }

    }

}
