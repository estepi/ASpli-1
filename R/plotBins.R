.plotBins <- function( 
    bin, 
    counts, 
    as, 
    factorsAndValues, 
    targets,
    main = NULL,
    binCountsColor  = '#2F7955',
    geneCountsColor = '#79552F',
    psirColor       = '#465579',
    j1CountsColor   = '#A04935',
    j2CountsColor   = '#752020',
    j3CountsColor   = '#A07C35',
    useBarplots     = NULL ) {
  
  # -------------------------------------------------------------------------- #
  # Manipulate factors
  # If just one factor, add another fictional factor with only one value to
  # make computation easier.
  justOneFactor <- FALSE
  if ( length ( factorsAndValues) == 1 ) {
    factorsAndValues <- append( factorsAndValues, list( single = 1 ) )
    targets <- cbind( targets, single = 1 )
    justOneFactor <- TRUE
    
  }
  
  factorsAndValues <- factorsAndValues[rev(names(factorsAndValues))]
  
  mainFactorIndex <- length( factorsAndValues )
  
  nPoints <- length( factorsAndValues[[mainFactorIndex]] )
  

    gridPanels <- expand.grid( factorsAndValues[-mainFactorIndex] )
    
    condInTargets <-  unique( apply( targets[ , names( gridPanels ) , drop= FALSE], 1,
            function(x) paste0(x,collapse = "_") ) )
    
    condInGrid <- apply( gridPanels, 1,
        function(x) paste0(x,collapse = "_") )
    
    gridPanels <- gridPanels[ condInGrid %in% condInTargets , , drop=FALSE ]
  
  # -------------------------------------------------------------------------- #
  
  
  # -------------------------------------------------------------------------- #
  # Extract data from bin
  binCounts <- .extractCountColumns( countsb( counts )[ bin,] , targets )
  geneLocus <- countsb( counts )[ bin,'locus']
  
  geneCounts <- .extractCountColumns( countsg( counts )[ geneLocus,] , targets )
  
  psir <- joint( as )[ bin, as.character( unique( targets$condition ) ) ]
  psir [ is.na(psir) ] <- 0
  
  J123 <- joint( as )[ bin,  colnames( joint( as ) ) %in%  rownames(targets)]
  ac <- 0
  J1   <- J123[ ( ac + 1 ) : ( ac + nrow(targets) ) ]
  ac <- ac + nrow(targets)
  J2   <- J123[ ( ac + 1 ) : ( ac + nrow(targets) ) ]
  ac <- ac + nrow(targets)
  J3   <- J123[ ( ac + 1 ) : ( ac + nrow(targets) ) ]
  # -------------------------------------------------------------------------- #
  
  
  prepareData <- function( samplesProfileData, gridPanels, targets, factorsAndValues ) {
    mainFactorIndex <- length( factorsAndValues )
    targets$samplesNames <- rownames( targets)
    
    data <- lapply ( 1:nrow( gridPanels ) , function ( i ) {
          samplesData <- merge( targets, gridPanels[i,], by=colnames(gridPanels[i,]) )
          plotData <- sapply ( factorsAndValues[[mainFactorIndex]], function(x) {
                samples <- samplesData[ 
                    samplesData[ , names(factorsAndValues[mainFactorIndex]) ] == x  , 
                    'samplesNames' ]
                rowMeans ( samplesProfileData[ rownames( targets ) %in% samples ] )
              } )
          plotData[ is.na( plotData )] <- 0
          return( plotData )
        } )
    return(data)
  }
  
  makePlot <- function ( data, useBarPlot = FALSE, main, col) {
    
    mainFactorIndex <- length( factorsAndValues )
    
    if ( ! justOneFactor ) {
      
      xTicksLabels <- rep( apply( gridPanels, 1, 
              function( x ) paste0(x,collapse=  ".") ), sapply( data, length  ) )
      
      xTicksLabels <- apply( gridPanels, 1, function( x ) paste0(x,collapse=  ".") )
      
      xTicksLabels <- expand.grid( factorsAndValues[[ mainFactorIndex]], xTicksLabels )
      
      xTicksLabels <- paste0( xTicksLabels[,2] ,".",  xTicksLabels[,1] )
    } else {
      xTicksLabels <- factorsAndValues[[mainFactorIndex]]
    }
    
    maxValue <- max ( sapply ( data, max ))
    minValue <- min ( sapply ( data, min ))
    
    dataCons <- Reduce( function( a, b) { a <- append( a, b )}, data )
    
    par( bty = 'n')
    
    if ( ! useBarPlot ) {
      
      plot( x = 1:length( dataCons ) , 
          y = dataCons,
          pch = 20,
          type = 'p',
          col = 'gray',
          ylim = c( minValue, maxValue),
          xlab='',
          xaxt = "n",
          ylab = "Counts",
          yaxt = 'n')
      
    } else {
      
      dataCons <- matrix( dataCons, ncol = nPoints, byrow = TRUE ) 
      barplot ( height = t( dataCons ),
          col = col,
          ylim = c( 0, maxValue),
          xlab='',
          yaxt = 'n',
          width=0.9,
          space = c(0.111111, 0.5),
          border = '#DDDDDDFF',
          beside = T,
          las = 2 )
      legend( 
          x = "topleft", 
          legend = main,
          adj=c( 0,0),
          bty="n",
          cex=1 )
      
      ticks <- matrix( 1, 
          ncol = ncol(dataCons),
          nrow = nrow(dataCons))
      ticks <- cbind( ticks, 0.35 )
      
      ticks <- as.vector( t( ticks ) )
      
      ac <- 0
      tickSum <- c()
      for ( i in 1 : (length( ticks ) -1) ) {
        tickSum <- append( tickSum , ac + ticks[i])
        ac <- tickSum[i]
      }
      
      tlblmat <- matrix ( xTicksLabels, byrow=TRUE, ncol = ncol(dataCons) ) 
      tlblmat <- cbind( tlblmat, '' )
      xTicksLabels <- as.vector( t(tlblmat) )[1:( length( tickSum ) )]
      
      axis( 
          side = 1, 
          at = tickSum - 0.11, 
          labels= xTicksLabels, las=2 , 
          cex.axis=0.7, tck=0,
          mgp=c(3,0.3,0),
          col= rgb(0.9,0.9,0.9, 0) )
      
    }
    axis( 2, las=1, cex.axis=0.7 )
    
    if (!useBarPlot) {
      secFactorTable <- table( gridPanels[ , ncol(gridPanels)] )
      st = 0
      for ( i in 1:nrow( secFactorTable )) {
        from <-  (st + 1)
        to <- st + secFactorTable[i] * nPoints
        at <- c(from:to)
        axis( side = 1, 
            at = at, 
            labels= xTicksLabels[from:to], las=2 , 
            cex.axis=0.7, tck=0,
            mgp=c(3,0.3,0),
            col= rgb(0.9,0.9,0.9) )
        st = st + to
      }
      
      st = 0
      for( i in 1:nrow( gridPanels )) { 
        lines( x = st + 1:length( factorsAndValues[[ mainFactorIndex ]] ) , 
            y = data[[i]],
            col = col)
        st <- st + length( factorsAndValues[[mainFactorIndex]] )
      }
    }
  }

  par( mfrow= c( 6 , 1 ) )
  
  par( oma = c( 0, 0, 2.4, 0),
       mar = c( 2.1, 3.1, 1.1, 1.1) )
  
  # Add Bin count data
  data <- prepareData( binCounts, gridPanels, targets, factorsAndValues )
  makePlot( data, TRUE, main = "Bin counts", col = binCountsColor )
  
  # Add Gene count data
  data <- prepareData( geneCounts, gridPanels, targets, factorsAndValues )
  makePlot( data, TRUE, main = "Gene counts", col = geneCountsColor)
  
  # Add PSI/PIR count data
  makePlot( psir, TRUE, main = "PSI / PIR" , col = psirColor)
  
  # Add J1 Exc count data
  data <- prepareData( J1, gridPanels, targets, factorsAndValues )
  makePlot( data, TRUE, main = "Incl. junc 1", col = j1CountsColor  )

  # Add J2 Exc count data
  data <- prepareData( J2, gridPanels, targets, factorsAndValues )
  makePlot( data, TRUE, main = "Incl. junc 2", col =  j2CountsColor )
  
  # Add J3 Exc count data
  data <- prepareData( J3, gridPanels, targets, factorsAndValues )
  makePlot( data, TRUE, main = "Exclu. junc", col = j3CountsColor)
  
  title( main = bin, outer=TRUE)
  
}


















