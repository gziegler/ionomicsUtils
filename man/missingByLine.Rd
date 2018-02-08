\name{missingByLine}
\alias{missingByLine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Takes a recoded (0,1,2) data.table and returns a single row data.table of length ncol with fraction of missing calls for each line
}
\description{
Takes a recoded (0,1,2) data.table and returns a single row data.table of length ncol with fraction of missing calls for each line
}
\usage{
missingByLine(genoTable)
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
Returns a data.table object with one row and ncol(genoTable) columns. The row contains the fraction of missing SNP calls for that genotype.
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
missingByLine(GenotypeData)
missingByLine[,which(missingByLine>0.2),with=F] #3 lines are missing more than 20% of the calls
snp.subset[,colnames(missingByLine[,which(missingByLine>0.2),with=F]) := NULL]

## The function is currently defined as
`missingByLine` <- function(genoTable) {
  return(genoTable[,lapply(.SD,function(x) length(which(is.na(x)))/length(x))])
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
