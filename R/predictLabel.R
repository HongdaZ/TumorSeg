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
    gamma_star_trace <- vector( "numeric", burnin + n + 1 )
    gamma_star_trace[ 1 ] <- gamma_star
    sd <- .1
    lambda_alpha <- rep( 1, 4 )
    lambda_beta <- rep( 1, 3 )
    lambda_gamma <- 1
    mu <- vector( "list", 4 )
    sigma <- vector( "list", 4 )
    sum_old <- sumPatient( old )
    initial <- sum_old
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
                if( new_alpha[ k ] > 0 ) {
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
            }
            # update beta's
            for( k in 1 : 3 ) {
                old_beta <- beta[ , j ]
                new_beta <- old_beta
                new_beta[ k ] <- old_beta[ k ] +
                    rnorm( 1, mean = 0, sd = sd )
                if( new_beta[ k ] > 0 ) {
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
            }
            # update gamma
            old_gamma <- gamma[ j ]
            new_gamma <- old_gamma + rnorm( 1, 0, sd )
            if( new_gamma > 0 ) {
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
        # get beta_sum and gamma_sum
        logZS( new, alpha_star, beta_star, gamma_star, new$beta_sum,
               new$gamma_sum, core )
        sum_new <- sumPatient( new[ 1 : 8 ] )
        # update star's
        # update alpha_star's
        for( k in 1 : 4 ) {
            old_alpha <- alpha_star
            new_alpha <- old_alpha
            new_alpha[ k ] <- old_alpha[ k ] +
                rnorm( 1, mean = 0, sd = sd )
            if( new_alpha[ k ] > 0 ) {
                proportion <- exp( ( old_alpha[ k ] - new_alpha[ k ] ) *
                                       ( lambda_alpha[ k ] +
                                             sum_new$n_type[ k ] ) -
                                       ( logZS( new, old_alpha,
                                                beta_star, gamma_star, new$beta_sum,
                                                new$gamma_sum, core ) -
                                             logZS( new, new_alpha,
                                                    beta_star,gamma_star,
                                                    new$beta_sum,
                                                    new$gamma_sum, core ) ) )
                rho <- min( 1, proportion )
                if( rbinom( 1, 1, rho ) == 1 ) {
                    alpha_star[ k ] <- new_alpha[ k ]
                    alpha_star_trace[ k, i + 1 ] <- new_alpha[ k ]
                } else {
                    alpha_star_trace[ k, i + 1 ] <- old_alpha[ k ]
                }
            } else {
                alpha_star_trace[ k, i + 1 ] <- old_alpha[ k ]
            }
        }
        # update beta_star's
        for( k in 1 : 3 ) {
            old_beta <- beta_star
            new_beta <- old_beta
            new_beta[ k ] <- old_beta[ k ] +
                rnorm( 1, mean = 0, sd = sd )
            if( new_beta[ k ] > 0 ) {
                proportion <- exp( ( old_beta[ k ] - new_beta[ k ] ) *
                                       ( lambda_beta[ k ] +
                                             new$beta_sum[ k ] ) -
                                       ( logZS( new, alpha_star,
                                                old_beta, gamma_star, new$beta_sum,
                                                new$gamma_sum, core ) -
                                             logZS( new, alpha_star,
                                                    new_beta, gamma_star,
                                                    new$beta_sum,
                                                    new$gamma_sum,
                                                    core ) ) )
                rho <- min( 1, proportion )
                if( rbinom( 1, 1, rho ) == 1 ) {
                    beta_star[ k ] <- new_beta[ k ]
                    beta_star_trace[ k, i + 1 ] <- new_beta[ k ]
                } else {
                    beta_star_trace[ k, i + 1 ] <- old_beta[ k ]
                }
            } else {
                beta_star_trace[ k, i + 1 ] <- old_beta[ k ]
            }
        }
        # update gamma_star
        old_gamma <- gamma_star
        new_gamma <- old_gamma + rnorm( 1, 0, sd )
        if( new_gamma > 0 ) {
            proportion <- exp( ( old_gamma - new_gamma ) *
                                   ( lambda_gamma  + new$gamma_sum ) -
                                   ( logZS( new, alpha_star,
                                            beta_star, old_gamma, new$beta_sum,
                                            new$gamma_sum, core ) -
                                         logZS( new, alpha_star,
                                                beta_star, new_gamma, new$beta_sum,
                                                new$gamma_sum, core ) ) )
            rho <- min( 1, proportion )
            if( rbinom( 1, 1, rho ) == 1 ) {
                gamma_star <- new_gamma
                gamma_star_trace[ i + 1 ] <- new_gamma
            } else {
                gamma_star_trace[ i + 1 ] <- old_gamma
            }
        } else {
            gamma_star_trace[ i + 1 ] <- old_gamma
        }

        # update lamdba_alpha
        lambda_alpha <- rgamma( 4, shape =  s + a + 1,
                                rate = rowSums( alpha ) + b + alpha_star )
        # update lambda_beta
        lambda_beta <- rgamma( 3, shape =  s + a + 1,
                               rate = rowSums( beta ) + b + beta_star )
        # update lambda_gamma
        lambda_gamma <- rgamma( 1, shape =  s + a + 1,
                                rate = sum( gamma ) + b + gamma_star )
        # update mu_k and sigma_k
        for( k in 1 : 4 ) {
            mean <- ( sum_old$sum_y[[ k ]] + sum_new$sum_y[[ k ]] ) /
                ( 1 + sum_old$n_type[ k ] + sum_new$n_type[ k ] )
            cov <- sigma[[ k ]] / ( 1 + sum_old$n_type[ k ] +
                                        sum_new$n_type[ k ] )
            mu[[ k ]] <- mvrnorm( 1, mean, cov )
        }
        for( k in 1 : 4 ) {
            v <- m + 1 + sum_old$n_type[ k ] + sum_new$n_type[ k ]
            scl<- tau_0 * diag( 1, nrow = 4, ncol = 4 ) +
                sum_old$sum_cross_y[[ k ]] +
                sum_new$sum_cross_y[[ k ]] -
                2 * mu[[ k ]] %*%
                t( sum_old$sum_y[[ k ]] + sum_new$sum_y[[ k ]] ) +
                ( 1 + sum_old$n_type[ k ] + sum_new$n_type[ k ] ) *
                mu[[ k ]] %*% t( mu[[ k ]] )
            sigma[[ k ]] <- riwish( v, scl )
        }
        det_sigma <- vector( "numeric", 4 )
        inv_sigma <-matrix( 0, 16, 4 )
        sigma_mu <- matrix( 0, 4, 4 )
        q_sigma <- vector( "numeric", 4 )
        for( l in 1 : 4 ) {
            det_sigma[ l ] <- det( sigma[[ l ]] )
            inv_sigma[ , l ] <- solve( sigma[[ l ]] )
            sigma_mu[ , l ] <- solve( sigma[[ l ]] ) %*% mu[[ l ]]
            q_sigma[ l ] <- t( mu[[ l ]] ) %*% solve( sigma[[ l ]] ) %*%
                mu[[ l ]]
        }
        # update x_star
        for( k in 1 : 8 ) {
            pmfLabel( alpha_star, beta_star, gamma_star,
                      det_sigma, inv_sigma, sigma_mu,
                      q_sigma, new[[ k ]]$neighbor_label,
                      new[[ k ]]$sub_modality_mat,
                      new[[ k ]]$prob, core )
            updatePred( new, k )
            len_type <- 4
            sum_cross_y <- vector( "list", len_type )
            sum_y <- vector( "list", len_type )
            n_type <- tabulate( new[[ k ]]$pred_seg, nbins = len_type )
            for( l in 1 : len_type ) {
                index <- new[[ k ]]$pred_seg == l
                y <- matrix( new[[ k ]]$sub_modality_mat[ , index ],
                             nrow = len_type )
                sum_cross_y[[ l ]] <- tcrossprod( y )
                sum_y[[ l ]] <- rowSums( y )
            }
            new[[ k ]]$sum_cross_y <- sum_cross_y
            new[[ k ]]$sum_y <- sum_y
            new[[ k ]]$n_type <- n_type
            if( k < 8 ) {
                next_seq <- k + 1
                fillNbrlabel( new$cube, new[[ next_seq ]]$neighbor_index,
                              new[[ next_seq ]]$neighbor_label, core )
            }
        }
        for( k in 1 : 8 ) {
            fillNbrlabel( new$cube, new[[ k ]]$neighbor_index,
                          new[[ k ]]$neighbor_label, core )
        }
    }
    list( new = new, alpha_star_trace = alpha_star_trace,
          beta_star_trace = beta_star_trace,
          gamma_star_trace = gamma_star_trace )
}
