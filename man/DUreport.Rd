\name{DUreport}
\alias{DUreport}
\title{
Differential gene expression and differential bin/junction usage estimation
}
\description{
Estimate differential expression at gene level and differential usage at bin and junction level.
}
\usage{
DUreport(counts, targets, minGenReads, minBinReads,
                 minRds, offset = FALSE, offsetAggregateMode,
                 offsetUseFitGeneX = TRUE,
                 contrast = NULL, 
                 forceGLM = FALSE, 
                 ignoreExternal =TRUE, 
                 ignoreIo = TRUE, 
                 ignoreI = FALSE)


}

\arguments{
\item{counts}{
An object of class ASpliCounts
}
\item{targets}{
A dataframe containing sample, bam and condition columns 
}
\item{minGenReads}{Default 10 reads}
\item{minBinReads}{Default 5 reads}
\item{minRds}{Default 0.05}
\item{ignoreExternal}{Ignore Exon Bins at the beggining or end of the transcript. Default TRUE}
\item{ignoreIo}{Ignore original introns. Default TRUE}
\item{ignoreI}{Ignore intron bins . Default FALSE}
\item{offset}{Default = FALSE}
\item{offsetAggregateMode}{Options:=c( “geneMode”, “binMode” )[ 2 ]},
\item{offsetUseFitGeneX}{Default= TRUE} 
\item{contrast}{ Default = NULL}
\item{forceGLM}{ Default = FALSE}

}


\value{ An ASpliDU object with results at \code{genes}, \code{bins} and \code{junctions} level }
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


