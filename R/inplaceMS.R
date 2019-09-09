inplaceMS <-function( A, B ) {
    storage.mode( A ) <- storage.mode( B ) <- "double"
    .Call( "inplaceMS", A, B )
}
