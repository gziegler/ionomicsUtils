\name{labelOutliers}
\alias{labelOutliers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a ggplot with outliers colored.
}
\description{
Not a fully fleshed out function, but somewhat recreating what plot.lm does in labeling outliers.
Labels outliers from a linear model in a plot (e.g. those with a std. resid. > 2) (note that plot.lm defaults to labelling 3 ids, regardless of extremity).
Takes a data frame with x, y, and label columns, performs an LM and extracts 'extreme' values and labels them on the plot.
}
\usage{
labelOutliers(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
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
df <- data.frame(id=sample(LETTERS[1:5],length(x),replace=TRUE),x,comp=x+y,z)
df$comp[4] <- 20
df$x[5] <- 10
labelOutliers(data.frame(x=df$x,y=df$comp,label=df$id))

function (df) 
{
    require(ggplot2)
    df <- na.omit(df)
    m = lm(df[, "x"] ~ df[, "y"])
    df.fortified = fortify(m)
    df$extreme = ifelse(abs(df.fortified$.stdresid) > 2, 1, 0)
    p <- ggplot(df, aes(x = x, y = y)) + geom_point() + geom_text(data = df[df$extreme == 
        1, ], aes(x = x, y = y, label = label), size = 3, hjust = -0.3) + 
        annotate("text", label = lm_eqn(df), parse = TRUE, x = Inf, 
            y = Inf, hjust = 1.1, vjust = 1.5)
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
