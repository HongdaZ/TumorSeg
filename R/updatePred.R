updatePred <- function( new_patient, i ) {
    # Sample the labels
    new_pred <- apply( new_patient[[ i ]]$prob, 2, sample, x = 1 : 4, size = 1,
                       replace = FALSE  )
    storage.mode( new_pred ) <- "numeric"
    .Call( "updatePred", new_patient[[ i ]]$pred_seg, new_pred )
    .Call( "updateCube", new_patient[[ 9 ]], new_patient[[ i ]]$sub_vec_index,
           new_pred )
    .Call( "updateCount", new_patient[[ i ]]$count, new_pred )
}
