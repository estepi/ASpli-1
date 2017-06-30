.extractCountColumns <- function ( aDataframe, targets ) {
  result <- aDataframe[ , match( row.names(targets), colnames( aDataframe ) ) ]
  colnames( result ) <- as.character( row.names(targets) )
  return( result )
}

.extractDataColumns <- function ( aDataframe, targets ) {
  result <- aDataframe[ , - match( row.names(targets) , colnames( aDataframe ) ) ]
  return( result )
}

.condenseTargetsConditions <- function ( targets ) {
  if( ! "condition" %in% colnames( targets ) ) {
    
    targets <- data.frame( 
        targets, 
        condition = apply( targets[ , -1 , drop = FALSE] ,1 ,paste,collapse="_"))
    
  }
  return( targets )
}

.DUreport <- function( 
    counts, 
    targets, 
    minGenReads  = 10,
    minBinReads  = 5,
    minRds = 0.05,
    offset = FALSE,
    offsetAggregateMode = c( "geneMode", "binMode" )[2],
    offsetUseFitGeneX = TRUE,
    contrast = NULL,
    forceGLM = FALSE,
    ignoreExternal = TRUE,
    ignoreIo = TRUE, 
    ignoreI = FALSE
    # ---------------------------------------------------------------------- #
    # Comment to disable priorcounts usage in bin normalization 
    # , priorCounts = 0 
    # ---------------------------------------------------------------------- #
  ) {
  
  # Create result object                   
  du <- new( Class="ASpliDU" )
  
  # Generate conditions combining experimental factors
  targets <- .condenseTargetsConditions( targets ) 
  
  # ------------------------------------------------------------------------ #
  # Filter genes and calculates differential usage of genes
  du <- .DUreportGenes( du, counts, targets, minGenReads, minRds, contrast, 
      forceGLM )
  message("Genes DE completed")
  # ------------------------------------------------------------------------ #
  
  # ------------------------------------------------------------------------ #
  # Filter bins and calculates differential usage of bins 
  du <- .DUReportBins( du, 
                       counts, 
                       targets, 
                       minGenReads, 
                       minBinReads, 
                       minRds,
                       offsetAggregateMode, 
                       offsetUseFitGeneX, 
                       offset, 
                       contrast, 
                       forceGLM, 
                       ignoreExternal, 
                       ignoreIo, 
                       ignoreI ) 
  message("Bins DE completed")
  # ------------------------------------------------------------------------ #
  
  return(du)
}

.DUreportBinSplice <- function (  
    counts, 
    targets, 
    minGenReads  = 10,
    minRds = 0.05,
    contrast = NULL,
    forceGLM = FALSE ) {
 
  # Create result object                   
  du <- new( Class="ASpliDU" )
  
  # Generate conditions combining experimental factors
  targets <- .condenseTargetsConditions( targets ) 
  
  # ------------------------------------------------------------------------ #
  # Filter genes and calculates differential usage of genes
  du <- .DUreportGenes( du, counts, targets, minGenReads, minRds, contrast, 
      forceGLM )
  message("Genes DE completed")
  # ------------------------------------------------------------------------ #
  
  # ------------------------------------------------------------------------ #
  # Filter bins and calculates differential usage of bins 
  du <- .DUReportBinsSpliceDGE(counts, targets, contrast )
  message("Bins DE completed")
  # ------------------------------------------------------------------------ #
  
}

.DUreportGenes <- function (
    du, 
    counts, 
    targets, 
    minGenReads  = 10,
    minRds = 0.05,
    contrast = NULL,
    forceGLM = FALSE ) {
  
  # Filter genes and calculates differential usage of genes
  dfG0 <- .filterByReads( df0 = countsg( counts ), targets = targets,
      min = minGenReads, type = "any" )
  
  dfGen <- .filterByRdGen( df0 = dfG0, targets = targets,
      min = minRds, type = "any" ) 
  
  genesde <- .genesDE( df=dfGen, targets = targets,
      contrast = contrast, forceGLM = forceGLM )
  
  du@genes <- genesde
  
  return( du )
  
}

.DUReportBins <- function( 
    du, 
    counts, 
    targets, 
    minGenReads, 
    minBinReads, 
    minRds,
    offsetAggregateMode, 
    offsetUseFitGeneX, 
    offset, 
    contrast, 
    forceGLM, 
    ignoreExternal, 
    ignoreIo, 
    ignoreI ) {
  
  # Filter bins
  dfG0 <- .filterByReads( df0 = countsg( counts ), targets = targets,
      min = minGenReads, type = "all" )
  
  dfGen <- .filterByRdGen( df0 = dfG0, targets = targets, 
      min = minRds, type = "all" )
  
  dfBin <- countsb(counts)[countsb(counts)[,"locus"]%in%row.names(dfGen),]
  
  if( ignoreIo ) dfBin <- dfBin[dfBin[,"feature"]!="Io",]
  
  df1 <- .filterByReads( df0=dfBin, targets=targets, 
      min=minBinReads, type="any" )
  
  df2 <- .filterByRdBinRATIO( dfBin=df1, dfGen=dfGen,
      targets=targets, min=minRds, type="any" )
  
  # Set offset matrix is required
  if( offset ) {
    mOffset <- .getOffsetMatrix(
        dfBin,
        dfGen,
        targets,
        offsetAggregateMode = offsetAggregateMode,
        offsetUseFitGeneX   = offsetUseFitGeneX )
    mOffset <- mOffset[ rownames( df2 ), ]
  } else {
    mOffset <- NULL
  }
  
  # Calculate bins DU
  binsdu <- .binsDU(
      df = df2,
      dfGen = dfGen,
      targets = targets,
      mOffset = mOffset,
      contrast = contrast,
      forceGLM = forceGLM,
      ignoreExternal = ignoreExternal,
      ignoreIo = ignoreIo, 
      ignoreI = ignoreI
      # -------------------------------------------------------------------- # 
      # Comment to disable priorcounts usage 
      , priorCounts = 0 
#    , priorCounts = priorCounts 
  # -------------------------------------------------------------------- # 
  ) 
  
  du@bins <- binsdu
 
  return( du )
}

.DUReportBinsSpliceDGE <- function( counts, targets, contrast ) {
  
  countData <- countsb( counts )
  countData <- countData[ countData[,'feature'] != "Io", ]
  
  group <- targets$condition
  
  if( is.null( contrast ) ) constrast <- .getDefaultContrasts(group)
  
  y <- DGEList( counts = .extractCountColumns( countData, targets ),
      group = targets$condition,
      genes = data.frame( 
          locus = countData$locus, 
          bin   = rownames( countData ) ) )
  
  y <- DGEList( counts = .extractCountColumns( countData, targets ),
      group = targets$condition,
      genes = .extractDataColumns(countData, targets) )       
  
  keep <- rowSums( cpm( y ) > 1) >= 2
  y <- y[ keep, , keep.lib.sizes = FALSE ]
  y <- calcNormFactors( y )
  
  design <- model.matrix( ~targets$condition )
  y      <- estimateDisp( y, design )
  fit    <- glmFit( y, design, contrast )
  captured <- capture.output(
      ds <- diffSpliceDGE( 
          fit, contrast = contrast, geneid = "locus", exonid = "exonid" ) )
  tsp    <- topSpliceDGE( ds, test = "exon", FDR = 1, number = INF )
  
  du@bins <- tsp
  
  return( du )
}


.filterByReads <- function( df0, targets, min, type, pair = NULL ) {

  # subset a working data frame
  cropped <- .extractCountColumns( df0, targets )
 
  # Modify working dataframe with pair being compared
  # TODO: pair no se requiere mas, modificar por constrastes?
  if ( is.null( pair ) ) {
    pair <- as.character( unique( targets$condition ) )
  }
  cropped <- cropped[ , targets$condition %in% pair ]
  
  # Calculates row means for conditions being compared
  list <- matrix( unlist( 
    lapply( pair, function( x ) {
    rowMeans( cropped[ , targets$condition == x ] ) >= min  } ) ),  
    nrow=nrow( cropped ), 
    byrow = FALSE )

  # Filter original gene counts 
  if (type == "all"){
    df=df0[ rowSums(list) == length( pair ), ]
  } else {
    df=df0[ rowSums(list) > 0,]
  } 
  
  return ( df )
}


.filterByRdGen <- function( df0, targets, min, type ) {
  
  dens <- .extractCountColumns( df0, targets ) / df0$effective_length
  
  list <- matrix( unlist(
    lapply( unique( targets$condition ), 
            function( x ) {
              rowMeans( dens[ , targets$condition == x ] ) >= min } )),  
    nrow = nrow( dens ), 
    byrow = FALSE )

  #keeps those genes which ave rd > min in any condition

  if (type=="any")  { 
    df <- df0[ rowSums( list ) > 0             , ]
  } else  {
    df <- df0[ rowSums( list ) == ncol( list ) , ]
  }
  
  return (df)
}


.getDefaultContrasts <- function ( conditions ) {
	contrast <- rep( 0, length( unique( conditions ) ) )
	contrast[1:2] <- c(-1,1)
  return( contrast )
}

.genesDE <- function( df, targets, contrast = NULL, forceGLM = FALSE ) { 
  
  if( is.null( contrast ) ) contrast <- .getDefaultContrasts(targets$condition)
  
  cols <- match( rownames( targets ), colnames( df ) )
  
  group <- targets$condition
  
  er <- DGEList( counts = df[ , cols ], samples=targets, group=group)
  er <- calcNormFactors( er )
  
  justTwoConditions <- sum( contrast != 0 ) == 2
  
  if( justTwoConditions & ! forceGLM ){
    capture.output( er   <- estimateDisp( er ) )
    pair <- which( contrast != 0 )
    et   <- exactTest(er, pair=pair)
  } else {
    design <- model.matrix( ~0 + group, data = er$samples )
    captured <- capture.output(
      er     <- estimateDisp( er, design = design )
    )
    glf    <- glmFit( er, design = design)  
    et     <- glmLRT( glf, contrast = contrast)
  } 
  
  fdr.gen <- p.adjust( et$table$PValue, method="BH" )
  
  cols <- match(rownames( targets ), colnames( df ) )
  geneData <- .extractDataColumns( df, targets )
  genesFull <- data.frame( geneData ,
      logFC = as.numeric( et$table$logFC ), 
      pvalue = as.numeric( et$table$PValue ), 
      gen.fdr = as.numeric(fdr.gen), 
      stringsAsFactors = FALSE)
  rownames( genesFull ) <- rownames( df )
  return( genesFull )
}


.filterByRdBinRATIO <- function( 
    dfGen,
    dfBin,
    targets, 
    min, 
    type ) {
  

  # -------------------------------------------------------------------------- #
  # Get densities for genes and bins
  genes.rd <- .extractCountColumns( dfGen, targets ) / dfGen$effective_length 
  bins.rd  <- .extractCountColumns( dfBin, targets ) / dfBin$length
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Calculates avg. bin dentity by condition
  avRdBin <- matrix( unlist(
          lapply( unique( targets$condition ), 
              function( x ) { 
                rowMeans( bins.rd[ , targets$condition == x] ) 
              } ) ),  
      nrow = nrow( bins.rd ), 
      byrow = FALSE) 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Calculates avg. gene dentity by condition
  avRdGen <- matrix( unlist(
          lapply( unique( targets$condition ), 
              function(x) rowMeans(genes.rd[ , targets$condition == x]))),  
      nrow = nrow(genes.rd), 
      byrow = FALSE)#Ok
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Calculates bin to gene read density ratio
  gen.rdb <- avRdGen[ match( dfBin$locus, rownames( dfGen ) )  , ]
  bin.gen.rd <-  avRdBin / gen.rdb
  bin.gen.rd[ is.na( bin.gen.rd ) ] <- 0 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Apply filter
  if ( type == "any" )  { 
    dfBin <- dfBin[ rowSums( bin.gen.rd >= min ) > 0 ,]
  } else{ 
    dfBin <- dfBin[ rowSums( bin.gen.rd >= min ) == ncol( bin.gen.rd ) ,]
  }
  # -------------------------------------------------------------------------- #
  
return (dfBin)
}

# ---------------------------------------------------------------------------- #
.normalizeByGenFeature <- function( feature, gene, targets, priorCounts = 0 ) {

  # Search the column with the gene name in the features 
  colLocus <- match( c( "gene", "locus" ), colnames( feature ) )
  colLocus <- colLocus[ ! is.na( colLocus ) ][1]
  
  # Repeat gene rows by the number of bins in that gene 
  f.index <- match( feature[ , colLocus ], rownames(gene) ) 
  counts.genes.f <- gene[f.index,]

  # Extract count data
  counts.genes.f  <- .extractCountColumns( counts.genes.f, targets)
  counts.features <- .extractCountColumns( feature, targets) 
  
  # Calculate mean of gene counts 
  gen.mean.f <- rowMeans( counts.genes.f )   
  
  # Apply normalization
  normalizedCounts <- round( ( as.matrix(  counts.features ) + priorCounts ) / 
                       ( as.matrix( counts.genes.f ) + priorCounts ) * 
                       ( gen.mean.f + priorCounts ) )
  
  # Set invalid results to zero
  normalizedCounts[ is.na( normalizedCounts ) | is.infinite( normalizedCounts) ] <- 0 
  
  # Create result dataframe
  normalizedFull <- cbind( .extractDataColumns( feature , targets ) , normalizedCounts )

  return( normalizedFull )
}
# ---------------------------------------------------------------------------- #


.binsDU <- function (
      df,
      dfGen,
      targets,                   
      ignoreExternal= TRUE, 
      ignoreIo = TRUE, 
      ignoreI = FALSE,
      contrast = NULL,
      forceGLM = FALSE,
      mOffset  = NULL,
      priorCounts = 0 ) {

  # -------------------------------------------------------------------------- #
  # Filtrar bins de acuerdo de diferentes criterios
  df = df[ ! ignoreExternal | df$event != "external" ,] 
  df = df[ ! ignoreIo | df$feature != "Io" ,] 
  df = df[ ! ignoreI | df$feature != "I" ,] 
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Normalize bins by gen or filter offsets
  if( is.null( mOffset ) ){
    df <- .normalizeByGenFeature( feature=df, gene=dfGen, targets = targets, 
        priorCounts = priorCounts )
  } else {
    mOffset <- mOffset[ rownames( df ), ]
  }
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # perform edgeR extact or glm test 
  et <- .edgeRtest( df, dfGen, targets, mOffset, contrast, forceGLM )
  # -------------------------------------------------------------------------- #
  
  # -------------------------------------------------------------------------- #
  # Build result dataframe
  bin.fdr <- p.adjust( et$table$PValue, method="BH" )
  logFC   <- et$table$logFC
  pvalue  <- et$table$PValue
  
  cols <- match( rownames( targets ), colnames( df ) )

  splicing_full <- data.frame ( df[ , -cols ],
      logFC, 
      pvalue,
      bin.fdr,
      stringsAsFactors = FALSE )
  
  rownames( splicing_full ) <- rownames( df )
  
  return( splicing_full )
  # -------------------------------------------------------------------------- #
  
}

.junctionsDU_SUM <- function( df, 
                              dfGen,
                              targets, 
                              mOffset = NULL,
                              contrast = NULL,
                              forceGLM = FALSE,
                              priorCounts = 0 ) {

  # -------------------------------------------------------------------------- #
  # Inner function to compute junction ratio
  jratio <- function( junctions ) {
    junctions[ is.na( junctions ) ] <- 0 
    return( junctions[1] / ( junctions[1] + junctions[2] ) )
  }
  # -------------------------------------------------------------------------- #
  
  if( is.null( mOffset ) ) {
    df <- .normalizeByGenFeature( feature=df, gene=dfGen, targets, priorCounts )
  }
  
  et <- .edgeRtest( df, dfGen, targets, mOffset, contrast, forceGLM )
  
  fdr <- p.adjust( et$table$PValue, method="BH" )
  logFC <- et$table$logFC
  pvalue <- et$table$PValue
  
  group<- targets$condition
  
  jranges <- .createGRangesExpJunctions( rownames( df ) )

# TODO:  Es Necesario seguir preguntando por una version vieja de IRanges ?
# Ademas, las dos llamadas a la version vieja llaman diferente al parametro
# ignore.redundant ( en la otra llamada es ignoreRedundant )
#  if ( packageVersion("IRanges") < 2.6) {
#    j.start <- findOverlaps( jranges, ignoreSelf=TRUE, ignore.redundant=FALSE,
#        type="start")  
#  } else {
    j.start <- findOverlaps( jranges, drop.self=TRUE, drop.redundant=FALSE,
        type="start")
#  }
  
  jjstart <- as.data.frame( j.start )
  
  jjstart$queryHits <- names(jranges[jjstart$queryHits])
  jjstart$subjectHits <- names(jranges[jjstart$subjectHits])
  shareStart <- data.frame( aggregate( subjectHits ~ queryHits, 
          data = jjstart, paste, collapse=";"))
  
  
  start <- ncol( df ) - nrow(targets) + 1 
  end   <- ncol( df ) 

  dfCountsStart <- data.frame( names=jjstart$queryHits, 
      df[ jjstart$subjectHits, start:end],
      row.names=NULL) 
  
  dfSumStart <- data.frame(aggregate(. ~ names, data = dfCountsStart, sum))
  sumJ <- paste(colnames(dfSumStart), "jsum", sep=".")
  colnames(dfSumStart) <-  sumJ
  rownames(dfSumStart) <- dfSumStart$names.jsum
  
  dfSumStart$names.jsum <- NULL

  dffStart <- data.frame(matrix(NA, nrow =  nrow(df), ncol = ncol(dfSumStart)) )
  rownames(dffStart) <- rownames(df)
  colnames(dffStart) <- colnames(dfSumStart)  
  mSumStart <- match( row.names(dfSumStart), row.names(dffStart)) 
  
  dffStart[mSumStart,] <-  dfSumStart
  dffStart[is.na(dffStart)] <- 0
  
  mbin_start_hit <- match(shareStart$queryHits, row.names(dffStart))
  #aca reacomodo el bin_start_hit con el indexado de dffStart
  bin_start_hit <- as.character(rep("-", nrow(dffStart)) )
  bin_start_hit[mbin_start_hit] <- shareStart$subjectHits
  ################################################
  
  cols <- match( rownames( targets ),colnames(df))
  
  
  ratioStart <- data.frame( df[,cols],dffStart)
  colnames(ratioStart) <- rep(rownames( targets ),2)
  #aca hay que armar un df itnermedio con la suma por condicion:

 

  ff <- rep(group,2)
  colnames(ratioStart) <- paste(ff, rep(1:2,each=length(group)))
  dfSum <- t(apply(ratioStart, 1, function(x){tapply(as.numeric(x), 
                INDEX=colnames(ratioStart), 
                sum)}))
  colnames(dfSum) <- rep(unique(group), each=2) # old version
  
  jratioStartRes <- t(apply(dfSum, 1, function(x){tapply(as.numeric(x), 
                INDEX=colnames(dfSum), 
                jratio )}))

#  if (packageVersion("IRanges")<2.6) { 
#    j.end <- findOverlaps(jranges, 
#        ignoreSelf=TRUE,
#        ignoreRedundant=FALSE,
#        type="end")
#    
#  } else {
    j.end <- findOverlaps( jranges, 
        drop.self = TRUE,
        drop.redundant = FALSE,
        type = "end" )
#  }

  jjend <- as.data.frame(j.end)
  jjend$queryHits <- names(jranges[jjend$queryHits])
  jjend$subjectHits <- names(jranges[jjend$subjectHits])
  shareEnd <- data.frame( aggregate( subjectHits ~ queryHits, 
          data = jjend, paste, collapse=";") ) 
  dfCountsEnd <- data.frame( names=jjend$queryHits, 
      df[jjend$subjectHits,start:end],
      row.names=NULL) #recover counts
  dfSumEnd <- data.frame(aggregate(. ~ names, data = dfCountsEnd, sum))   
  sumJ <- paste(colnames(dfSumEnd), "jsum", sep=".")
  colnames(dfSumEnd) <- sumJ
  rownames(dfSumEnd) <- dfSumEnd$names.jsum
  dfSumEnd$names.jsum <- NULL
  dffEnd =data.frame(matrix(NA, nrow =  nrow(df), 
          ncol = ncol(dfSumEnd)) )
  rownames(dffEnd) <- rownames(df)
  colnames(dffEnd) <- colnames(dfSumEnd)
  ########################################################################
  mSumEnd <- match(row.names(dfSumEnd), row.names(dffEnd))
  dffEnd[mSumEnd,] <- dfSumEnd
  dffEnd[is.na(dffEnd)] <- 0
  ########################################################################
  mbin_end_hit <- match(shareEnd$queryHits, row.names(dffEnd))
  bin_end_hit <- rep("-", nrow(dffEnd))
  bin_end_hit[mbin_end_hit] <- shareEnd$subjectHits
  ########################################################################
  
  ratioEnd <- data.frame(df[,cols],dffEnd)
  ff <- rep(group,2)
  colnames(ratioEnd) <- paste(ff,rep(1:2, each = length(group)))
  dfSum <- t(apply(ratioEnd, 1, function(x){tapply(as.numeric(x), 
                INDEX=colnames(ratioEnd), sum  )}))

  colnames(dfSum) <- rep(unique(group), each = 2) # new version
  
  jratioEndRes <- t(apply(dfSum, 1, function(x){tapply(as.numeric(x), 
                INDEX=colnames(dfSum), jratio )}))

  
  et_merge <- data.frame(df,                       
      logFC,
      pvalue, 
      fdr,
      bin_start_hit,
      dffStart,
      jratioStartRes,
      bin_end_hit,
      dffEnd,
      jratioEndRes )
  return(et_merge)
  
}

# ---------------------------------------------------------------------------- #
# Calculates the matrix of offset values for normalization
# TODO: Que es el valor 10^-4 que aparece, se puede pasar como parametro?
.getOffsetMatrix <- function( df, dfGen, targets, 
    offsetAggregateMode = c( "geneMode","binMode" )[2], 
    offsetUseFitGeneX = TRUE) {
  
  locus  <- df[,"locus"]
  
  if( offsetAggregateMode=="geneMode" ) {
    if(offsetUseFitGeneX){
      countData = dfGen[,rownames(targets)]
      
      yg <- DGEList( counts = countData, 
                     group  = targets$condition, 
                     genes  = data.frame( locus = rownames(countData) ) )
      
      #filter lowcount genes
      keep <- rowSums(cpm(yg)>1) >= 2
      yg <- yg[keep, , keep.lib.sizes=FALSE]
      
      yg     <- calcNormFactors(yg)
      fc     <- targets$condition
      design <- model.matrix(~0+fc)
      yg     <- estimateDisp(yg,design)
      fitg   <- glmFit(yg,design)
      maux   <- fitg$fitted.values + 10^-4
    } else {
      maux   <- dfGen[,rownames(targets)] + 10^-4
    }
  } else {
    a <- by( data = df[,c("feature",rownames(targets))],
             INDICES = as.factor( locus ),
             FUN = function(x){ 
               apply( x [x[,1] == "E", 2:ncol(x) ], 2 , sum )
             })
    maux<- matrix(unlist(a),byrow=TRUE,ncol=nrow(targets))+10^-4
    rownames(maux)<-names(a)
    colnames(maux)<-names(a[[1]])
  }
  mOffset <- do.call( rbind, list( A = maux[ locus ,] ) )
  rownames( mOffset )<-rownames(df)
  return(mOffset)
}

.setDefaultOffsets <- function ( aDGEList , mOffset) {
	logNi <- apply( aDGEList$samples[,c("lib.size","norm.factors")],1,
			function(x){ log( prod( x ) ) } )
	
	# TODO: Check Alternative code 
	# logNi <- log( aDGEList$samples[,c("lib.size")] * aDGEList$samples[,"norm.factors"] )
	
	aDGEList$offset <- log(mOffset) + matrix( rep( logNi, nrow( mOffset )), 
			byrow = TRUE, ncol = length( logNi ) )
	return ( aDGEList )
}

.edgeRtest <- function( 
    df,
    dfGen,
    targets,
    mOffset = NULL,
    contrast = NULL,
    forceGLM = FALSE ) {
  
  cols <- match( rownames( targets ), colnames( df ) )
  
  group <- targets$condition
  
  if( is.null( contrast ) ) contrast <- .getDefaultContrasts(group)
  
  er <- DGEList( counts  = df[,cols],
                 samples = targets,
                 group   = group )
  
  er <- calcNormFactors(er)
  
  justTwoConditions <- sum( contrast != 0 ) == 2
  
  # TODO: Forzar GLM no tiene efecto si se pasa un offset. 
  if( ( !forceGLM ) & is.null( mOffset ) & justTwoConditions ){

    captured <- capture.output(
      er   <- estimateDisp( er )
    )
    pair <- which( contrast != 0 )
    testResult   <- exactTest( er, pair = pair )
    
#    message( "ExactTest... done" )
    
  } else {
    
    if( ! is.null( mOffset ) ) er <- .setDefaultOffsets( er, mOffset ) 
    
    design     <- model.matrix( ~0 + group, data = er$samples )
    er         <- estimateDisp( er, design = design )      
    glf        <- glmFit( er, design = design )  
    testResult <- glmLRT( glf, contrast = contrast )
    
#    message( "glmLRT... done" )

  } 
  return( testResult )
}


