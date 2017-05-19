\name{DUreport}

\alias{DUreport}

\title{
  Differential gene expression and differential bin usage estimation
}
\description{
  Estimate differential expression at gene level and differential usage at bin 
  level.
}

\usage{
  DUreport( counts, 
            targets, 
            minGenReads = 10, 
            minBinReads = 5,
            minRds = 0.05, 
            offset = FALSE, 
            offsetAggregateMode =  c( "geneMode", "binMode" )[2],
            offsetUseFitGeneX = TRUE,
            contrast = NULL, 
            forceGLM = FALSE, 
            ignoreExternal = TRUE, 
            ignoreIo = TRUE, 
            ignoreI = FALSE)
}

\arguments{
  \item{ counts }{ 
    An object of class ASpliCounts
  }
  \item{ targets }{
    A dataframe containing sample, bam and experimental factor columns 
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
    Default value is TRUE. 
  }
  \item{ignoreIo}{ Ignore original introns.  Default TRUE }
  \item{ignoreI}{Ignore intron bins, test is performed only for exons. 
    Default FALSE }
  \item{offset}{ Corrects bin expression using an offset model. Default = FALSE }
  \item{offsetAggregateMode}{ Options:=c( “geneMode”, “binMode” )[ 2 ] },
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
  \item{forceGLM}{ Default = FALSE}
}


\value{ 
  An ASpliDU object with results at \code{genes}, \code{bins} }

\author{
  Estefania Mancini, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz
}

\examples{

library(RNAseqData.HNRNPC.bam.chr14)
chr14 <- system.file("extdata","chr14.sqlite", package="ASpli")
genome <- loadDb(chr14)
features <- binGenome(genome)
targets <- data.frame(bam=RNAseqData.HNRNPC.bam.chr14_BAMFILES,
                       condition=c(rep("CT",4),rep("KD",4)))
bam <- loadBAM(targets)
counts <- readCounts ( features, bam, targets, cores = 1, readLength=100L, maxISize=50000)
du <- DUreport(counts, targets, offset = FALSE, ignoreExternal=TRUE, ignoreIo = TRUE)
writeDU(du, output.dir="only_du")
}
\seealso{
\code{
\link{DEXSeq}, \link{edgeR}
Accesors: \link{genesDE}, \link{binsDU},\link{junctionsDU}
Export: \link{writeDU}
}
}


