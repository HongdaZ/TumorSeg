sumPatient <- function( patient ) {
    n <- length( patient )
    sum_cross_y <- list( matrix( 0, nrow = 4, ncol = 4 ),
                         matrix( 0, nrow = 4, ncol = 4 ),
                         matrix( 0, nrow = 4, ncol = 4 ),
                         matrix( 0, nrow = 4, ncol = 4 ) )
    sum_y <- list( vector( "numeric", 4 ),
                   vector( "numeric", 4 ),
                   vector( "numeric", 4 ),
                   vector( "numeric", 4 ) )
    n_type <- vector( "integer", 4 )

    for( i in 1 : n ) {
        for( j in 1 : 4 ) {
            sum_cross_y[[ j ]] <- sum_cross_y[[ j ]] +
                patient[[ i ]]$sum_cross_y[[ j ]]
            sum_y[[ j ]] <-sum_y[[ j ]] +
                patient[[ i ]]$sum_y[[ j ]]
        }
        n_type <- n_type +
            patient[[ i ]]$n_type
    }
    list( sum_cross_y = sum_cross_y,
          sum_y = sum_y,
          n_type = n_type )
}
