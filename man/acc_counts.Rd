\name{Counts accesors}
\alias{countsb}
\alias{countse1i}
\alias{countsg}
\alias{countsie2}
\alias{countsj}
\alias{rdsg}
\alias{rdsb}
\alias{countsb<-}
\alias{countse1i<-}
\alias{countsg<-}
\alias{countsie2<-}
\alias{countsj<-}
\alias{rdsg<-}
\alias{rdsb<-}

\title{
Accessors for ASpliCounts object
}
\description{
Accessors for ASpliCounts object
}
\usage{
countsb(x)
countse1i(x)
countsg(x)
countsie2(x)
countsj(x)
rdsg(x)
rdsb(x)
}
\arguments{
  \item{x}{   
  An ASpliCounts object 
  }
}
\value{
Returns dataframes with counts by sample and genomic metadata
}
\author{ Estefania Mancini, Marcelo Yanovsky, Ariel Chernomoretz  }
\examples{

# Get data tables from an ASpliCounts object

counts <- aspliCountsExample() 

cb1  <- countsb(counts)
ce1i <- countse1i(counts)
cg   <- countsg(counts)
cie2 <- countsie2(counts)
cj   <- countsj(counts)
rg   <- rdsg(counts)
rb   <- rdsb(counts)

# Set data tables to an ASpliCounts object

countsb(counts)   <- cb1
countse1i(counts) <- ce1i
countsg(counts)   <- cg
countsie2(counts) <- cie2
countsj(counts)   <- cj
rdsg(counts)      <- rg
rdsb(counts)      <- rb

}


