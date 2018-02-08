\name{recodeGenoTable}
\alias{recodeGenoTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Take a data.table genotype file and call `recode` to convert to 0,1,2
}
\description{
Expects rows are SNPs and columns are lines. Data.table shouldn't have any metadata columns.
recode changes the most common base pair to be major, second most common to be minor. Anything else is a het, except N, which is coded as NA.
}
\usage{
recodeGenoTable(genoTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genoTable}{
Data.table with values to be converted. Expects rows are SNPs and columns are lines. Data.table shouldn't have any metadata columns. 
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
recodeGenoTable(genoTable)

## The function is currently defined as
function(genoTable){
  require(data.table)
  genoTable[, (names(genoTable)) := as.list(recode(.SD)), by=1:nrow(genoTable)]
  return(genoTable)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
