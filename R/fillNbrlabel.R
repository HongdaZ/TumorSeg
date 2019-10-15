fillNbrlabel <- function( cube, neighbor_index, neighbor_label, p = 12 ) {
    .Call( "fillNbrlabel", cube, neighbor_index, neighbor_label, p = 12 )
}
