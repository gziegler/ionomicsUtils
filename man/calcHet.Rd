\name{calcHet}
\alias{calcHet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate fraction of hets for each SNP in a genotype table.
}
\description{
Takes a recoded (0,1,2) data.table and returns a vector of fraction of het SNPs for each row
Calculates fraction of called SNPs that were called as heterozygous (value of 1)
Doesn't include NA calls in calculation
}
\usage{
calcHet(genoTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genoTable}{
Data table where SNPs are rows and lines are columns, no metadata columns. Coded as 0,1,2.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns of vector of same length as input with translated values.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Greg Ziegler
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
alleleTable <- cbind(snpInfo,calcHet(genoTable))
hetThreshold <- 0.1 #change hets to NA for a snp if there are higher than this percentage of hets
length(alleleTable$FracHet[which(alleleTable$FracHet>=0.1)])
alleleTable[,discard := (FracHet>=hetThreshold)]

GenotypeData[which(alleleTable$discard==TRUE)][GenotypeData[which(alleleTable$discard==TRUE)]==1] <- NA
alleleTable[, discard := NULL]

}
