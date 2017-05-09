\name{is.outlier}
\alias{is.outlier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove outliers from the provided vector and return vector with outliers changed to NA.
}
\description{
Remove outliers from vector using a MAD value as cutoff.
}
\usage{
is.outlier(x, mcut = 6.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vector of data points to remove outliers from
}
  \item{mcut}{
Number of MADs a data point needs to be from the median to be considered an outlier, default is 6.2
}
}
\details{
See: Davies, P.L. and Gather, U. (1993).
"The identification of multiple outliers" (with discussion)
J. Amer. Statist. Assoc., 88, 782-801.
}
\value{
A vector of same length as \code{"x"} with outliers changed to NA.
}
\references{
Davies, P.L. and Gather, U. (1993).
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
x <- c(-10, x, 10)
y <- is.outlier(x)
summary(x)
summary(y)
## The function is currently defined as
function (x, mcut = 6.2) 
{
    y <- na.omit(x)
    lims <- median(y) + c(-1, 1) * mcut * mad(y, constant = 1)
    for (j in 1:length(x)) {
        if (is.na(x[j]) | x[j] < lims[1] | x[j] > lims[2]) {
            x[j] <- NA
        }
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
