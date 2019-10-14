split3D <- function( dim, start ) {
    l1 <- ceiling( dim[ 1 ] / 2  )
    rows_1 <- ( 1 : l1 ) * 2 - 1
    l1 <- floor( dim[ 1 ] / 2  )
    rows_2 <- ( 1 : l1 ) * 2
    l2 <- ceiling( dim[ 2 ] / 2  )
    cols_1 <- ( 1 : l2 ) * 2 - 1
    l2 <- floor( dim[ 2 ] / 2  )
    cols_2 <- ( 1 : l2 ) * 2
    l3 <- ceiling( dim[ 3 ] / 2 )
    slices_1 <- ( 1 : l3 ) * 2 - 1
    l3 <- floor( dim[ 3 ] / 2 )
    slices_2 <- ( 1 : l3 ) * 2
    seq <- vector( "list", 8 )
    seq[[ 1 ]] <- t( expand.grid( i = rows_1, j = cols_1, k = slices_1 ) )
    seq[[ 2 ]] <- t( expand.grid( i = rows_1, j = cols_2, k = slices_1 ) )
    seq[[ 3 ]] <- t( expand.grid( i = rows_2, j = cols_1, k = slices_1 ) )
    seq[[ 4 ]] <- t( expand.grid( i = rows_2, j = cols_2, k = slices_1 ) )
    seq[[ 5 ]] <- t( expand.grid( i = rows_1, j = cols_1, k = slices_2 ) )
    seq[[ 6 ]] <- t( expand.grid( i = rows_1, j = cols_2, k = slices_2 ) )
    seq[[ 7 ]] <- t( expand.grid( i = rows_2, j = cols_1, k = slices_2 ) )
    seq[[ 8 ]] <- t( expand.grid( i = rows_2, j = cols_2, k = slices_2 ) )
    seq <- lapply( seq, function( x ) x - 1 + start )
    return( seq )
}
