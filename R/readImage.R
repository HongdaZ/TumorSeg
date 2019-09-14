readImage <- function( patient_file, type ) {
    seg_array <- readNifti( patient_file[[ 2 ]] )
    img_dim <- dim( seg_array )
    pairwise_clique <- pairwise
    higher_index <- higher_order_index
    G <- prod( img_dim )
    # "flair", "t1", "t1ce", "t2"
    modality_mat <- matrix( nrow =  G, ncol = 4 )
    modality_mat[ , 1 ] <- readNifti( patient_file[[ 1 ]] )
    for( i in 2 : 4 ) {
        modality_mat[ , i ] <- readNifti( patient_file[[ i + 1 ]] )
    }
    nonzero <- rep( TRUE, G )
    for ( i in 1 : 4 ) {
        nonzero <- nonzero & ( modality_mat[ , i ] != 0 )
    }
    modality_mat <- modality_mat[ nonzero, ]
    modality_mat <- scale( modality_mat )
    seg_valid <- seg_array[ nonzero ]
    seg_array[ !nonzero ] <- NA
    # sum of cross products of y's and sum of y's
    sum_cross_y <- vector( "list", 4 )
    sum_y <- vector( "list", 4 )
    n_type <- tabulate( seg_valid + 1, nbins = 5 )[ -4 ]
    for( i in 1 : 4 ) {
        index <- seg_valid == type[ i ]
        y <- modality_mat[ index, ]
        sum_cross_y[[ i ]] <- crossprod( y )
        sum_y[[ i ]] <- colSums( y )
    }
    # neighborhood informations
    vec_index <- which( nonzero == TRUE )
    # reshape modality_mat to be of 4 rows
    dim( modality_mat ) <- c( 4, length( vec_index ) )

    matrix_index <- vec2array( vec_index, dim = img_dim )
    last <- length( vec_index )
    up_mat <- matrix( rep( img_dim, dim( pairwise_clique )[ 2 ] ), nrow = 3 )
    valid_pos <- vector( "integer", dim( up_mat )[ 2 ] )

    five_base <- matrix( 1, nrow = 5, ncol = 8 )
    ten_base <- vector( "numeric", 26 )
    for( i in 1 : 6 ) {
        ten_base[ i ] <- 10 ^ 7
    }
    for( i in 7 : 18 ) {
        ten_base[ i ] <- 10 ^ 8
    }
    for( i in 19 : 26 ) {
        ten_base[ i ] <- 10 ^ 6
    }
    for( i in 1 : 7 ) {
        five_base[ , i ] <- 5 ^ ( 8 - i )
    }
    type_count <- matrix( nrow = 5, ncol = 8 )
    shift <- matrix( 1, nrow = 7, ncol = 8 )
    for( i in 1 : 7 ) {
        shift[ , i + 1 ] <- i * 5 + 1
    }
    # pb <- tkProgressBar( "readImage", "Analyzing data",
    #                      0, 100, 0, width = 600 )
    # step <- last / 100
    # n_step <- 1

    for( i in 1 : last ) {
        # pairwise neighbor
        pairwise_array <- pairwise_clique +
            matrix( rep( matrix_index[ i, ],
                         dim( pairwise_clique )[ 2 ] ), nrow = 3 )
        # find out of bound voxel
        valid_index <- pairwise_array > 0 & pairwise_array <= up_mat
        pairwise_array[ !valid_index ] <- NA
        valid_pos <- pairwise_array[ 1, ] +
            img_dim[ 1 ] * ( pairwise_array[ 2, ] - 1 ) +
            img_dim[ 1 ] * img_dim[ 2 ] * ( pairwise_array[ 3, ] - 1 )
        pairwise_seg <- seg_array[ valid_pos ]
        # 8th order neighbor
        higher_seg <- pairwise_seg[ higher_index ]

        type_count <- tabulate( higher_seg + shift, 40 )
        dim( type_count ) <- c( 5, 8 )
        type_count <- type_count != 0

        for( j in 1 : 4 ) {
            # pairwise
            unequal_label <- pairwise_seg != type[ j ]
            res <- sum( unequal_label * ten_base, na.rm = T )
            # number of distinct cell types in a higher-order clique
            type_count2 <- type_count
            type_count2[ type[ j ] + 1, ] <- T
            res <- res + sum( type_count2 * five_base )

            modality_mat[ j, i ] <- res
        }
        # if( i > ( n_step * step ) ) {
        #     info <- sprintf( "Analyzing data: %d%%",
        #                      n_step )
        #     setTkProgressBar( pb, info, "readImage", info )
        #     n_step <- n_step + 1
        # }

    }
    # close( pb )

    list( sum_cross_y = sum_cross_y,
          sum_y = sum_y,
          n_type = n_type,
          neighbor_info = as.integer( modality_mat ) )
}
