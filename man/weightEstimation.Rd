\name{weightEstimation}
\alias{weightEstimation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict weights of samples from ionomics data.
}
\description{
Takes a data frame of values with weights included, and a data frame of values that haven't been weighed. Finds list of stable elements (based on the assumption most samples should have a similar elemental concentraiton, so low RSDs), and uses those elements to estimate weights.
}
\usage{
weightEstimation(x, y, predictx, maxRSD = 25, maxEl = 10, df = 500/9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data frame of weight normalized elelmental concentrations.
}
  \item{y}{
Vector of weights for concentrations in x.
}
  \item{predictx}{
Data frame of values to predict weights for.
}
  \item{maxRSD}{
Maximum RSD threshold for elements to be included in estimation.
}
  \item{maxEl}{
Maximum number of elements to include in estimation.
}
  \item{df}{
Dilution factor for weight normalized concentrations.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, predictx, maxRSD = 25, maxEl = 10, df = 500/9) 
{
    normX <- x
    for (i in colnames(x)) {
        normX[, i] <- (x[, i]/y) * (df)
    }
    elRSD <- sapply(normX, sd, na.rm = TRUE)/sapply(normX, mean, 
        na.rm = TRUE) * 100
    elMean <- sapply(normX, mean, na.rm = TRUE)
    meetCutoff <- elRSD[which(elRSD < maxRSD)]
    if (length(meetCutoff) < maxEl) {
        numEl <- length(meetCutoff)
    }
    else {
        numEl <- maxEl
    }
    bestElements <- names(sort(meetCutoff)[1:numEl])
    weightEstimates <- predictx
    for (i in colnames(predictx)) {
        weightEstimates[, i] <- (predictx[, i]/elMean[i]) * (df)
    }
    weightEstimates <- as.data.frame(t(apply(weightEstimates[bestElements], 
        1, is.outlier, mcut = 3)))
    finalPrediction <- rowMeans(weightEstimates[bestElements], 
        na.rm = TRUE)
    return(list(weights = as.vector(finalPrediction), els = bestElements))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
