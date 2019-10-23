readImage <- function( patient_file ) {
    type <- c( 0, 1, 2, 4 )
    seg_array <- readNifti( patient_file[[ 2 ]] )
    # change label to 1 ~ # of type
    changeLabel( seg_array, type )
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
    len_type <- length( type )
    sum_cross_y <- vector( "list", len_type )
    sum_y <- vector( "list", len_type )
    n_type <- tabulate( seg_valid, nbins = len_type )
    beta_sum <- vector( "numeric", 3 )
    gamma_sum <- vector( "numeric", 1 )
    for( i in 1 : len_type ) {
        index <- seg_valid == i
        y <- modality_mat[ index, ]
        sum_cross_y[[ i ]] <- crossprod( y )
        sum_y[[ i ]] <- colSums( y )
    }
    # neighborhood informations
    vec_index <- which( nonzero == TRUE )
    modality_mat <- matrix( nrow = len_type, ncol = length( vec_index ) )
    matrix_index <- vec2array( vec_index, dim = img_dim )
    last <- length( vec_index )
    up_mat <- matrix( rep( img_dim, dim( pairwise_clique )[ 2 ] ), nrow = 3 )
    valid_pos <- vector( "integer", dim( up_mat )[ 2 ] )

    five_base <- matrix( 1, nrow = len_type, ncol = 8 )
    ten_base <- vector( "numeric", 26 )
    ten_base[ 1 : 6 ] <- 10 ^ 7
    ten_base[ 7 : 18 ] <- 10 ^ 8
    ten_base[ 19 : 26 ] <- 10 ^ 6
    for( i in 1 : 7 ) {
        five_base[ , i ] <- 5 ^ ( 8 - i )
    }
    type_count <- matrix( nrow = len_type, ncol = 8 )
    shift <- matrix( 0, nrow = 7, ncol = 8 )
    for( i in 1 : 7 ) {
        shift[ , i + 1 ] <- i * len_type
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

        type_count <- tabulate( higher_seg + shift, 8 * len_type )
        dim( type_count ) <- c( len_type, 8 )
        type_count <- type_count != 0

        for( j in 1 : len_type ) {
            # pairwise
            unequal_label <- pairwise_seg != j
            res <- sum( unequal_label * ten_base, na.rm = T )
            # number of distinct cell types in a higher-order clique
            type_count2 <- type_count
            type_count2[ j, ] <- T
            res <- res + sum( type_count2 * five_base )
            modality_mat[ j, i ] <- res
            # add information for alpha, beta
            if( seg_valid[ i ] == j ) {
                beta_sum[ 1 ] <- beta_sum[ 1 ] + sum( unequal_label[ 1 : 6 ],
                                                      na.rm = T )
                beta_sum[ 2 ] <- beta_sum[ 2 ] + sum( unequal_label[ 7 : 18 ],
                                                      na.rm = T )
                beta_sum[ 3 ] <- beta_sum[ 3 ] +
                    sum( unequal_label[ 19 : 26 ], na.rm = T )
                gamma_sum <- gamma_sum + sum( log( colSums( type_count2 ) ) )
            }
        }
        # if( i > ( n_step * step ) ) {
        #     info <- sprintf( "Analyzing data: %d%%",
        #                      n_step )
        #     setTkProgressBar( pb, info, "readImage", info )
        #     n_step <- n_step + 1
        # }

    }
    # close( pb )

    storage.mode( modality_mat ) <- "integer"
    list( sum_cross_y = sum_cross_y,
          sum_y = sum_y,
          n_type = n_type,
          neighbor_info = modality_mat,
          beta_sum = beta_sum,
          gamma_sum = gamma_sum )
}
