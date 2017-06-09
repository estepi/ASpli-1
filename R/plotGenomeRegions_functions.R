.getBamMatrix <- function( conditionMatrix, mergedFiles, gene ) { 

  apply( conditionMatrix, c(1,2), 
      function( x ) {
        a <- mergedFiles[ mergedFiles$genes == gene & mergedFiles$cond == x , 'bams' ]
        if (any( is.na(a) ) ) { return(NA) }
        return(a)
      }
  )
} 

.mergeBamsByCondition <- function ( targets , regions , tempFolder = './tmp' ) {
  
  destBams <- expand.grid( cond = getConditions( targets ), genes = names( regions) )
  
  destBams$bams <- file.path( tempFolder, paste0( 'cond.',destBams[,1], '.gene.',destBams[,2],".bam") ) 
  
  file.exists( tempFolder) || dir.create( tempFolder , recursive = TRUE) 
  
  mapply( function( cond, gene, bam ) { 
        sourceBams <- as.character( targets$bam[ targets$condition == cond ] )
        mergeBam( files = sourceBams , destination = bam , region = regions[ gene ], overwrite = TRUE )
        indexBam( bam )
      } ,
      destBams$cond, destBams$genes, destBams$bams
  )
  return( destBams )
}


.definePlottingRegions <- function( counts, x, xIsBin ) {
  
  binCounts <-  countsb( counts )
  geneCounts <- countsg( counts )

  if ( xIsBin ) {
    selectedGenes <- unique( countsb( counts )[ rownames(  countsb( counts ) ) %in% x ,'locus'] )
  } else {
    selectedGenes <- x
  }
  if ( length( genes ) > 0 ) {
    
    geneCoordinates <- geneCounts [ selectedGenes, 'gene_coordinates', drop = FALSE ]
    
    geneCoordinates <- as.character( geneCoordinates[ 
            !duplicated( geneCoordinates), 1 ] )
    geneCoordinates <- strsplit( geneCoordinates, '[-:]')
    geneCoordinates <- do.call( rbind , geneCoordinates )
    geneCoordinates <- data.frame(geneCoordinates, stringsAsFactors = FALSE)
    
    geneCoordinates[,2] <- as.integer( geneCoordinates[,2] )
    geneCoordinates[,3] <- as.integer( geneCoordinates[,3] )
    
    regions <- GRanges( geneCoordinates[,1], 
        IRanges( geneCoordinates[,2] , geneCoordinates[,3]) , strand='+')
    
    names( regions ) <- selectedGenes   
    
    return( regions )
  
  } else {
    
    return( GRanges( NULL, IRanges( NULL, NULL ) ) )
    
  }
  
  
}

.arrangeLayout <- function( targets, layout, colors, plotTitles ) {
  
  # Check if auto layout is required
  if ( is.character( layout ) && tolower( layout ) == 'auto' ) {
    
    .colorsFromConditionMatrix <- function( conditionMatrix , colors ) {
      if ( length( colors ) == 1 && tolower( colors ) == 'auto' ) {
        colors <- colors()[1:(length(conditionMatrix)) + 15]
        matrix( colors, ncol = ncol(conditionMatrix) )
      } else {
        matrix( colors, 
            ncol = ncol( conditionMatrix ),
            nrow = nrow( conditionMatrix ),
            byrow = TRUE)
      }
    }
    
    expFactors <- colnames( targets )[ 
        ! colnames( targets ) %in% c( 'bam', 'condition' ) ]
    
    nExpFactors <- length( expFactors )
    
    # ------------------------------------------------------------------------ #
    # Case 1: just one experimental factor
    # result matrix is 1 x n
    if ( nExpFactors == 1 ) {
      conditionMatrix <- matrix( getConditions(targets), 
          ncol = length( getConditions(targets) ) )
      plotTitles <- conditionMatrix
      colors <- .colorsFromConditionMatrix( conditionMatrix, colors )
    } 
    # ------------------------------------------------------------------------ #
    
    # ------------------------------------------------------------------------ #
    # Case 2: more than one experimental factor and n * m * ... conditions. 
    # result matrix is n * ( m * ... ), where n is the number of values of the 
    # first factor and m * ... are the number of combination of remaining 
    # factors.
    if ( nExpFactors > 1 ) {

      nComb <- prod( apply( targets[, expFactors], 2 , 
              function( x ) length( unique( x ) ) ) )
      
      if ( length( getConditions( targets ) == nComb ) ) {
        conditionMatrix <- matrix( getConditions( targets ), 
            nrow = length( unique( targets[expFactors][,1]) ),
            byrow = TRUE )
        plotTitles <- conditionMatrix
        colors <- .colorsFromConditionMatrix(conditionMatrix, colors)
        
      } else {
        conditionMatrix <- matrix( getConditions( targets ), 
            ncol = length( getConditions(targets) ) )
        plotTitles <- conditionMatrix     
        colors <- .colorsFromConditionMatrix(conditionMatrix, colors)
      }
    } 
    # ------------------------------------------------------------------------ #
    
  } else {
    
    # ------------------------------------------------------------------------ #
    # Case 3: No auto arrange required. 
    # Matching matrix dimensions of layout, colors and plotTitles matrixes is 
    # checked.
    if ( ! is.matrix( layout ) ) { 
      stop(simpleError("Layout must be a matrix or 'auto'" ))
    }
    if ( is.null( colors ) ) {
      colors <- .colorsFromConditionMatrix( layout, colors )
    }
    if ( ncol( plotTitles ) != ncol(layout) || 
        nrow( plotTitles ) != nrow(layout) ||
        ncol( colors ) != ncol(layout) || 
        nrow( colors ) != nrow(layout) ) {
      stop( simpleError(
              "The number of columns and rows of colors and plotTitles must match to layout" ))
    }
    conditionMatrix <- layout
    # ------------------------------------------------------------------------ #
  }
  result <- list()
  result$conditionMatrix <- conditionMatrix
  result$plotTitles <- plotTitles
  result$colors <- colors
  return( result )
}

.plotGenomeRegions <- function( x, genomeTxDb, counts, targets, xIsBin = TRUE, 
    layout = 'auto', colors = 'auto', plotTitles = 'auto', sashimi = FALSE, 
    zoomOnBins= FALSE, deviceOpt = NULL, highLightBin = TRUE, outfolder = NULL, 
    outfileType = 'png', mainFontSize = 24, annotationHeight = 0.2,
    annotationCol = 'black', annotationFill = 'gray', annotationColTitle = 'black' ) {
  
  targets <- .condenseTargetsConditions( targets )
  
  # -------------------------------------------------------------------------- #
  # Autoarrange 
  layoutPars <- .arrangeLayout( targets, layout, colors, plotTitles )
  conditionMatrix <- layoutPars$conditionMatrix
  colors <- layoutPars$colors
  print( colors )
  plotTitles <- layoutPars$plotTitles
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # Keep existing genes and bins in the data set
  if ( xIsBin ) {
    x <- x[ x %in% rownames( countsb( counts ) ) ]
  } else {
    x <- x[ x %in% rownames( countsg( counts ) ) ]
  }
  if ( length( x ) == 0 ) {
    stop( simpleError( "Bin and/or genes names are incorrect." ) )
  } 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Collect and merge bam files
  regions <- .definePlottingRegions( counts, x, xIsBin)
  mergedFiles <- .mergeBamsByCondition( targets, regions )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Plot
  for ( xi in x ) {

    currentGene <- if ( xIsBin ) countsb( counts )[ xi, c('locus') ] else xi
    highLightBin <- highLightBin & xIsBin
    zoomOnBins <- zoomOnBins & xIsBin
    bamfiles <- .getBamMatrix( conditionMatrix , mergedFiles , currentGene )
    genLims <- c( start( regions[currentGene] ), end( regions[currentGene] ) )
    binLims <- if ( xIsBin ) as.integer( countsb( counts )[ 
                  xi, c( 'start', 'end' ) ] ) else NULL 
    
    if ( ! is.null( outfolder ) ) {
      outfile <- file.path( outfolder, 
          .makeValidFileName( paste0( xi,'.gr.', outfileType ) ) )
    } else { 
      outfile <- NULL
    }
    
    .makeGenomeRegionPlot( 
        main = xi, 
        genLims = genLims ,
        binLims = binLims ,
        chromosome = seqlevels(regions[currentGene]),
        outfile = outfile,
        genome = genomeTxDb,
        bamFiles = bamfiles,
        plotNames = conditionMatrix,
        colors = colors,
        sashimi = sashimi,
        zoomOnBins = zoomOnBins,
        outfileType = outfileType,
        deviceOpt = deviceOpt,
        highLightBin = highLightBin,
        annotationHeight = annotationHeight,
        annotationCol = annotationCol,
        annotationFill = annotationFill,
        annotationColTitle = annotationColTitle
    )
  }
  # -------------------------------------------------------------------------- #
}

.makeGenomeRegionPlot <- function ( 
    main , 
    genLims , 
    binLims , 
    chromosome ,  
    genome , 
    outfile , 
    bamFiles = NULL ,    # a matrix of H x V 
    plotNames = NULL ,   # a matrix of H x V
    colors = NULL ,      # a matrix of H x V  
    mainFontSize = 24,
    sashimi = TRUE , 
    zoomOnBins = FALSE, 
    individualWidth = 1024,
    individualHeight = 800,
    highLightBin = TRUE,
    outfileType = 'png',
    deviceOpt = NULL,
    annotationHeight = 0.2,
    annotationCol = 'black',
    annotationFill = 'gray',
    annotationColTitle = 'black') {
  
  # -------------------------------------------------------------------------- #
  # Extrae la cantidad de plots horizontales y verticales
  hplots <- ncol ( bamFiles )
  vplots <- nrow ( bamFiles )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Zoom on Bins
  if (zoomOnBins) {
    if ( genLims[2] - genLims[1] - ( binLims[2] - binLims[1] ) > 2000 ) {
      genLims[1] <- binLims[1] - as.integer( abs(binLims[2]-binLims[1])*0.1 ) 
      genLims[2] <- binLims[2] + as.integer( abs(binLims[2]-binLims[1])*0.1 )
    } 
  } else {
    genLims[1] <- genLims[1]- as.integer( abs(genLims[2]-genLims[1])*0.1 ) 
    genLims[2] <- genLims[2]+ as.integer( abs(genLims[2]-genLims[1])*0.1 ) 
  }
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # Functions to define tracks
  alnTrack <- function( rangeFile, name, colors, size  ) {
    AlignmentsTrack( 
        range = rangeFile,
        isPaired = FALSE, 
        name=name,
        col= colors,
        fill.coverage = colors,
        col.coverage = colors,
        background.title = 'transparent',
        type = c( 'coverage', 'sashimi' )[c( TRUE, sashimi ) ],
        sashimiHeight = 0.5,
        coverageHeight = 0.5,      
        lwd.sashimiMax = 3,
        lwd.sashimiMin = 0.5,        
        showTitle = TRUE,
        lwd = 0,
        col.axis = 'black',
        cex.axis = 0.6,
        cex.title = 0.6,
        col.title = 'black',
        size = size,
        margin= 20
      )
  }
  
  binHighLight <- function ( tracks ) {
    tracks <- tracks[ ! is.na (tracks )]
    HighlightTrack(
        trackList = tracks, 
        start = binLims[1], 
        end = binLims[2],
        chromosome = chromosome,
#        col="transparent",
        col=rgb(0.5,0.5,0.5,0.3),
        fill=rgb(0.5,0.5,0.5,0.25),
        lwd=1,
        frame = TRUE,
        inBackground=FALSE
    )
  }
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Define tracks
  alntracks <- lapply(1:hplots, function ( y ) { lapply(1:vplots, function( z ) { NA } ) } )
  for ( h in 1:hplots ) {
    for ( v in 1:vplots) {
      if ( ! is.na( bamFiles[ v, h ] ) ) {
        cTrack <- alnTrack ( 
            rangeFile = bamFiles [ v, h ],
            name      = plotNames[ v, h ],
            colors    = colors   [ v, h ],
            size      = 10000 * ( 1 - annotationHeight )  / sum( ! is.na( bamFiles[,h]))) 
        alntracks[[ h ]][[ v ]] <- cTrack 
      } 
    }
  }
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  geneAnnotationTrack <- GeneRegionTrack(genome, 
      chromosome = chromosome,
      start = genLims[1], 
      end = genLims[2],
      options(ucscChromosomeNames=FALSE),
      name = "", 
      transcriptAnnotation = "symbol",
      size = annotationHeight * 10000,  
      fill = annotationFill ,
      col = annotationCol,
      shape = "arrow",
      col.title = annotationColTitle
  )
  # -------------------------------------------------------------------------- #
  

  
  if ( highLightBin ) {
    columnTracks <- lapply( alntracks, function (x) {
        binHighLight ( x )
      } )
  } else {
    columnTracks <- alntracks
  }

  
  # -------------------------------------------------------------------------- #
  # creates the graphic device
  outputIsAValidFile <- ( ! is.null( outfile ) ) && 
      ( file.access( dirname( outfile ) , 2 ) == 0 ) 
  
  if ( outputIsAValidFile ) {
    
    devicesNames <- c( 'png', 'bmp', 'jpeg', 'tiff', 'pdf')
    devices      <- list( png, bmp, jpeg, tiff, pdf)
    
    deviceIndex <- match( tolower( outfileType ), devicesNames )
    
    if ( is.na(deviceIndex) ) { 
      dev.new()
      message( paste('Format',outfileType,'is not recognized.',
              'Select png, bmp, jpeg, tiff or pdf') )
      
    } else {
      device <- devices[[ deviceIndex ]]
      undefDeviceOpt <-  is.null( deviceOpt )
      
      if ( (! is.na( deviceIndex )) && deviceIndex %in%  c(1:4) ) {
        fileNameOpt <- list( filename = outfile)
        defaultDeviceOpt <- list( res = 200, height = individualHeight * vplots, 
            width = individualWidth  * hplots, 
            pointsize = 12, units = 'px')
      }
      
      if ( (! is.na( deviceIndex )) && deviceIndex ==  5 ) {
        fileNameOpt <- list( file = outfile)
        defaultDeviceOpt <- list( height = 8 * vplots, width = 8 * hplots, 
            pointsize = 12)        
      }
      
      deviceOpt <- if ( undefDeviceOpt ) 
            append( defaultDeviceOpt, fileNameOpt ) 
          else 
            append( deviceOpt, fileNameOpt )
      do.call( device, deviceOpt )    
      
    }
    
  } else {
    dev.new()
    if ( ! is.null ( outfile ) ) {
      message( paste( "File:",outfile,"cannot be created." ) )
    }
  }
  # -------------------------------------------------------------------------- #
  
  grid.newpage()
  
  pushViewport( 
      viewport( layout = grid.layout( 
              nrow = 2, 
              ncol = hplots , 
              heights = unit( c ( 1, 5 ) , "null") ) ) )

  grid.text(main, gp=gpar( fontsize = mainFontSize ) , 
      vp = viewport(layout.pos.row = 1, layout.pos.col = 1:hplots) )
  
  first <- TRUE
  for ( i in 1:hplots) {
    pushViewport(viewport(layout.pos.col=i, layout.pos.row=2))
    plotTracks( c( columnTracks [[ i ]] ,geneAnnotationTrack ), 
        chromosome        = chromosome, 
        from              = genLims[1], 
        to                = genLims[2], 
        min.height        = 0, 
        minCoverageHeight = 0, 
        min.width         = 2, 
        min.distance      = 5,
        background.title  = 'transparent',
        margin            = 10,
        innerMargin       = 5,
        add               = TRUE
    )
    first <- FALSE
    popViewport(1)
  }
  
  # -------------------------------------------------------------------------- #
  # Close graphic device 
  if ( outputIsAValidFile ) {
    output <- capture.output( dev.off() )
  } 
  # -------------------------------------------------------------------------- #
}
