updatePred <- function( new_patient, i ) {
    len_type <- 4
    # Sample the labels
    new_pred <- apply( new_patient[[ i ]]$prob, 2, sample, x = 1 : len_type,
                       size = 1,
                       replace = FALSE  )
    storage.mode( new_pred ) <- "numeric"
    .Call( "updatePred", new_patient[[ i ]]$pred_seg, new_pred )
    .Call( "updateCube", new_patient$cube, new_patient[[ i ]]$sub_vec_local,
           new_pred )
    .Call( "updateCount", new_patient[[ i ]]$count, new_pred )
}
