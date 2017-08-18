\name{normalizedBinCounts}

\alias{normalizedBinCounts}

\title{Extraction of bin counts normalized to genes}

\description{
 \code{normalizedBinCounts} method computes bin counts normalized by gene 
 counts as in \code{DUreport} method.
}

\usage{
  normalizedBinCounts( 
        counts, 
        targets, 
        minGenReads = 10,
        minBinReads = 5,
        minRds = 0.05,
        ignoreExternal = TRUE,
        ignoreIo = TRUE,
        ignoreI = FALSE,
        verbose = TRUE )
}

\arguments{
   \item{counts}{ 
    An object of class ASpliCounts
  }
  \item{targets}{
    A data.frame containing sample, bam and experimental factor columns.
  }
  \item{minGenReads}{
    Genes with at least an average of \code{minGenReads} reads for any condition 
    are included into the differential expression test.
    Bins from genes with at least an average of \code{minGenReads} reads for all
    conditions are included into the differential bin usage test.
    Default value is 10 reads.
  }
  \item{minBinReads}{ 
    Bins with at least an average of \code{minGenReads} reads for any condition 
    are included into the differential bin usage test.
    Default value is 5 reads.
  }
  \item{minRds}{
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
  \item{ignoreI}{Ignore intron bins. 
    Default FALSE }
  \item{verbose}{
    A logical value that indicates that detailed information about each step in 
    the analysis will be presented to the user. }
}

\author{
  Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz 
}

\value{
  A data frame containing feature, event, locus, locus_overlap, symbol, 
  gene coordinates, start of bin, end of bin, bin length, and normalized counts
  for each sample.
}

\seealso{
  \code{ \link{ASpliDU} },
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
  normBins  <- normalizedBinCounts( counts, targets )                   
  
}
