logZ <- function( patient, alpha, beta, gamma, p ) {
    .Call( "logZ", patient, alpha, beta, gamma, p )
}
