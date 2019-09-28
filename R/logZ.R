logZ <- function( patient, alpha, beta, gamma ) {
    unequal2 <- patient$neighbor_info %/% 10 ^ 8
    r <- patient$neighbor_info %% 10 ^ 8
    unequal1 <- r %/% 10 ^ 7
    r <- r %% 10 ^ 7
    unequal3 <- r %/% 10 ^ 6
    r <- r %% 10 ^ 6
    type_count <- vector( "list", 8 )

}
