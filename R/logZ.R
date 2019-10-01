logZ <- function( patient, alpha, beta, gamma, s, p ) {
    .Call( "logZ", patient, alpha, beta, gamma, s, p )
}
