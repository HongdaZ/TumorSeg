readNew <- function( patient_file, start, size ) {
    type <- c( 0, 1, 2, 4 )
    true_seg_array <- readNifti( patient_file[[ 2 ]] )
    # change label to 1 ~ # of type
    changeLabel( true_seg_array, type )
    img_dim <- dim( true_seg_array )
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
    modality_mat[ nonzero, ] <- scale( modality_mat[ nonzero, ] )
    true_seg_array[ !nonzero ] <- 0
    vec_index <- 1 : G
    groups <- split3D( size, start )
    groups_local <- split3D( size, c( 1, 1, 1 ) )
    n_group <- length( groups )
    new_patient <- vector( "list", n_group + 1 )
    up_mat <- matrix( rep( size, dim( neighbor )[ 2 ] ), nrow = 3 )
    for( i in 1 : n_group ) {
        # index in whole image
        sub_vec_index <- apply( groups[[ i ]], 2, array2vec, dim = img_dim  )
        sub_nonzero <- nonzero[ sub_vec_index ]

        # Remove invalid voxels
        sub_vec_index <- sub_vec_index[ sub_nonzero ]
        # return NULL if majority of the voxels are invalid
        if( length( sub_vec_index ) == 0 ) {
            cat( "Too little brain cells in this area!\n" )
            return()
        }
        sub_mat_local <- groups_local[[ i ]][ , sub_nonzero ]
        sub_vec_local <- apply( sub_mat_local, 2, array2vec, dim = size )
        sub_modality_mat <- t( modality_mat[ sub_vec_index, ] )
        sub_true_seg <- true_seg_array[ sub_vec_index ]
        # normal cell as default predicted value for new image
        pred_seg <- rep( 1, length( sub_vec_index ) )

        # sum of cross products of y's and sum of y's
        len_type <- 4
        sum_cross_y <- vector( "list", len_type )
        sum_y <- vector( "list", len_type )
        n_type <- tabulate( pred_seg, nbins = len_type )
        for( j in 1 : len_type ) {
            index <- pred_seg == j
            y <- sub_modality_mat[ , index ]
            sum_cross_y[[ j ]] <- tcrossprod( y )
            sum_y[[ j ]] <- rowSums( y )
        }
        neighbor_index <- matrix( 0, nrow = 26,
                                  ncol = dim( sub_mat_local )[ 2 ] )
        for( j in 1 : dim( sub_mat_local )[ 2 ] ) {
            neighbor_mat_local <- sub_mat_local[ , j ] + neighbor
            valid_index <- neighbor_mat_local > 0 &
                neighbor_mat_local <= up_mat
            neighbor_mat_local[ !valid_index ] <- NA
            valid_pos <- neighbor_mat_local[ 1, ] +
                size[ 1 ] * ( neighbor_mat_local[ 2, ] - 1 ) +
                size[ 1 ] * size[ 2 ] * ( neighbor_mat_local[ 3, ] - 1 )
            valid_pos[ is.na( valid_pos ) ] <- 0
            valid_pos[ is.na( valid_pos ) ] <- 0
            neighbor_index[ , j ] <- valid_pos
        }
        neighbor_label <- matrix( 0, nrow = dim( neighbor_index )[ 1 ],
                                  ncol = dim( neighbor_index )[ 2 ] )
        count <- matrix( 0, nrow = 4,
                         ncol = dim( neighbor_index )[ 2 ] )
        prob <- matrix( 0, nrow = 4,
                         ncol = dim( neighbor_index )[ 2 ] )
        new_patient[[ i ]] <- list( pred_seg = pred_seg,
                                    sub_true_seg = sub_true_seg,
                                    sub_vec_index = sub_vec_index,
                                    sub_vec_local = sub_vec_local,
                                    sub_modality_mat = sub_modality_mat,
                                    neighbor_index = neighbor_index,
                                    neighbor_label = neighbor_label,
                                    sum_cross_y = sum_cross_y,
                                    sum_y = sum_y,
                                    n_type = n_type )


    }
    cube <- array( 1, dim = size )
    new_patient[[ n_group + 1 ]] <- cube
    new_patient
}
