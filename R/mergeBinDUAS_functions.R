.mergeBinDUAS <- function( du, as, targets, contrast = NULL ) { 
  
  targets <- .condenseTargetsConditions( targets )
  bas <- joint( as )
  bdu <- binsDU( du )
  bas$event <- NULL
  
  if ( is.null( contrast )) contrast <- .getDefaultContrasts( targets$condition )
  
  psir <- t( bas[ rownames(bdu) , getConditions( targets ) ] )
  psir <- psir * contrast
  psir <- t( psir[ contrast != 0, ] )
  psir <- data.frame( delta = rowSums( psir ) )
  
  cbind( bdu, bas[rownames(bdu),], psir )
  
}