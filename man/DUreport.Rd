\name{DUreport}

\alias{DUreport}

\title{
  Differential gene expression and differential bin usage estimation
}
\description{
  Estimate differential expression at gene level and differential usage at bin 
  level. When targets has only two conditions, and contrast is not set, the 
  estimation of differential expression and usage is done with an exact test,
  otherwise is estimated using a generalized linear model.    
}

\usage{
  DUreport( counts, 
            targets, 
            minGenReads = 10, 
            minBinReads = 5,
            minRds = 0.05, 
            offset = FALSE, 
            offsetAggregateMode =  c( "geneMode", "binMode" )[1],
            offsetUseFitGeneX = TRUE,
            contrast = NULL, 
            forceGLM = FALSE, 
            ignoreExternal = TRUE, 
            ignoreIo = TRUE, 
            ignoreI = FALSE,
            filterWithContrasted = FALSE,
            verbose = FALSE)
}

\arguments{
  \item{ counts }{ 
    An object of class ASpliCounts
  }
  \item{ targets }{
    A data.frame containing sample, bam and experimental factor columns.
  }
  \item{ minGenReads }{
    Genes with at least an average of \code{minGenReads} reads for any condition 
    are included into the differential expression test.
    Bins from genes with at least an average of \code{minGenReads} reads for all
    conditions are included into the differential bin usage test.
    Default value is 10 reads.
  }
  \item{ minBinReads }{ 
    Bins with at least an average of \code{minGenReads} reads for any condition 
    are included into the differential bin usage test.
    Default value is 5 reads.
  }
  \item{ minRds }{
    Genes with at least an average of \code{minRds} read density for any 
    condition are included into the differential expression test.
    Bins from genes with at least an average of \code{minRds} read density for 
    all conditions are included into the differential bin usage test.
    Bins with at least an average of \code{minRds} read density for any 
    condition are included into the differential bin usage test.
    Default value is 0.05.
    }
  \item{ignoreExternal}{
    Ignore Exon Bins at the beginning or end of the transcript. 
    Default value is TRUE. }
  \item{ignoreIo}{ Ignore original introns.  Default TRUE }
  \item{ignoreI}{Ignore intron bins, test is performed only for exons. 
    Default FALSE }
  \item{offset}{ Corrects bin expression using an offset matrix derived from 
    gene expression data. Default = FALSE}
  \item{offsetAggregateMode}{ Choose the method to aggregate gene counts to 
    create the offset matrix. When \code{offsetAggregateMode} is 'geneMode' and 
    option \code{offsetUseFitGeneX} is TRUE, a generalized linear model is used 
    to create the offset matrix. When \code{offsetAggregateMode} is 'geneMode' 
    and option \code{offsetUseFitGeneX} is FALSE, the offset matrix is generated
    by adding a prior count to the gene count matrix. When 
    \code{offsetAggregateMode} is 'binMode' a matrix from obtained from the sum
    of exonic bin counts, this only takes those bins that passes filters using 
    minGenReads, minBinReads and minRds.
    Options:=c( "geneMode", "binMode" )[ 1 ] },
  \item{offsetUseFitGeneX}{ Default= TRUE } 
  \item{contrast}{ Define the comparison between conditions to be tested.
    \code{contrast} should be a vector with length equal to the number of 
    experimental conditions defined by \code{targets}. The values of this vector
    are the coefficients that will be used to weight each condition, the order 
    of the values corresponds to the order given by \code{getConditions} 
    function. When \code{contrast} is NULL, defaults to a vector containing -1,
    as the first value, 1 as the second an zero for all the remaining values, 
    this corresponds to a pair comparison where the first condition is assumed 
    to be a control and the second condition is the treatment, all other 
    conditions are ignored. 
    Default = NULL}
  \item{forceGLM}{ Force the use of a generalized linear model to estimate 
    differential expression and usage. Default = FALSE }
  \item{filterWithContrasted}{
    A logical value specifying if bins, genes and junction will be filtered by
    read quantity and read density using data from those conditions that will 
    be used in the comparison, i.e. those which coefficients in contrast 
    argument are different from zero. The default value is FALSE, it is strongly
    recommended to do not change this value.}
  \item{verbose}{
    A logical value that indicates that detailed information about each step in 
    the analysis will be presented to the user. }
}

\value{ 
  An ASpliDU object with results at \code{genes}, \code{bins} level.}

\author{
  Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
}

\examples{
  # Create a transcript DB from gff/gtf annotation file.
  # Warnings in this examples can be ignored. 
  library(GenomicFeatures)
  genomeTxDb <- makeTxDbFromGFF( system.file('extdata','genes.mini.gtf', 
                                 package="ASpli") )
  
  # Create an ASpliFeatures object from TxDb
  features <- binGenome( genomeTxDb )
  
  # Define bam files, sample names and experimental factors for targets.
  bamFileNames <- c( "A_C_0.bam", "A_C_1.bam", "A_C_2.bam", 
                     "A_D_0.bam", "A_D_1.bam", "A_D_2.bam" )
  targets <- data.frame( 
               row.names = paste0('Sample_',c(1:6)),
               bam = system.file( 'extdata', bamFileNames, package="ASpli" ),
               factor1 = c( 'C','C','C','D','D','D') )
  
  # Load reads from bam files
  bams <- loadBAM( targets )
  
  # Read counts from bam files
  counts   <- readCounts( features, bams, targets, cores = 1, readLength = 100, 
                          maxISize = 50000 )
  
  # Calculate differential usage of genes and bins
  du       <- DUreport( counts, targets )

  # Export results  
  writeDU( du = du, output.dir = "only_du" )
}

\seealso{
  \code{\link{edgeR}}, \code{\link{junctionDUreport}}
  Accessors: \code{\link{genesDE}}, \code{\link{binsDU}}
  Export: \code{ \link{writeDU} }
}


