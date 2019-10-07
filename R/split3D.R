split3D <- function( dim ) {
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
    seq[[ 1 ]] <- expand.grid( rows_1, cols_1, slices_1 )
    seq[[ 2 ]] <- expand.grid( rows_1, cols_2, slices_1 )
    seq[[ 3 ]] <- expand.grid( rows_2, cols_1, slices_1 )
    seq[[ 4 ]] <- expand.grid( rows_2, cols_2, slices_1 )
    seq[[ 5 ]] <- expand.grid( rows_1, cols_1, slices_2 )
    seq[[ 6 ]] <- expand.grid( rows_1, cols_2, slices_2 )
    seq[[ 7 ]] <- expand.grid( rows_2, cols_1, slices_2 )
    seq[[ 8 ]] <- expand.grid( rows_2, cols_2, slices_2 )
    return( seq )
}
