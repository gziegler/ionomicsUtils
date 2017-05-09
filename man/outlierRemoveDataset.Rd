\name{outlierRemoveDataset}
\alias{outlierRemoveDataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calls is.outlier for all of the data columns (cols) of the provided data.frame (x)
}
\description{
Calls \code{"is.outlier"} for all of the data columns (cols) of the provided data.frame (x)
and returns the data frame with NA in place of outliers

}
\usage{
outlierRemoveDataset(x, mcut = 6.2, by = NA, cols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data.frame with sample data, metadata, etc.
}
  \item{mcut}{
Number of MADs a data point need to be from the median to be considered an outlier, default is 6.2
}
  \item{by}{
Column name to group data by for outlier removal (e.g. by line, by run, etc.), if not provided then by whole dataset
}
  \item{cols}{
Vector of column numbers or names in x to remove outliers from.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns data frame in the same format as input, but with outliers in each of the specified columns changed to NA.
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
set.seed(1)
x <- rnorm(100)
y <- rnorm(100)
z <- rnorm(100)
x <- c(-10, x, 10)
y <- c(-20, y, 20)
z <- c(-30, z, 30)
df <- data.frame(id=sample(LETTERS[1:5],length(x),replace=TRUE),x,y,z)
#By entire dataset
dfOR <- outlierRemoveDataset(df,6.2,by=NA,c("x","y","z"))
summary(dfOR)
#Look for outliers within groups
dfOR <- outlierRemoveDataset(df,6.2,by="id",c("x","y","z"))
summary(dfOR)

## The function is currently defined as
function (x, mcut = 6.2, by = NA, cols) 
{
    for (i in cols) {
        if (is.na(by)) {
            x[, i] <- is.outlier(x[, i], mcut)
        }
        else {
            for (j in unique(x[, by])) {
                if (is.na(j)) {
                  x[is.na(x[, by]), i] <- is.outlier(x[is.na(x[, 
                    by]), i], mcut)
                }
                else {
                  x[x[, by] == j & !(is.na(x[, by])), i] <- is.outlier(x[x[, 
                    by] == j & !(is.na(x[, by])), i], mcut)
                }
            }
        }
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
