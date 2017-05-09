\name{lm_eqn}
\alias{lm_eqn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return a linear model equation with R^2 value
}
\description{
Return the linear model equation of y~x formatted for plotting.
Takes the first column of the data frame as the dependent variable and the second column as the independent variable.
}
\usage{
lm_eqn(df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
Two column data frame, first column is the dependent variable (y) and second column is the independent variable (x)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a latex formatted string with the linear model fit equation and r^2 value.
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
x <- c(-10, x, 10)
y <- c(-20, y, 20)
df <- data.frame(id=sample(LETTERS[1:5],length(x),replace=TRUE),x,y)
lm_eqn(df[,3:2])
###Put it on a plot###
library(ggplot2)
ggplot(df,aes(x=x,y=y)) + geom_point() + geom_smooth(method="lm") +annotate("text", label=lm_eqn(df[,3:2]), parse=TRUE, x=Inf, y=Inf, hjust=1.1, vjust=1.5)

## The function is currently defined as
function (df) 
{
    m = lm(df[, 1] ~ df[, 2])
    paste("italic(y)==", format(coef(m)[1], digits = 2), "+", 
        format(coef(m)[2], digits = 2), "\%.\%italic(x)*\",\"~~italic(r)^2==", 
        format(summary(m)$r.squared, digits = 3), sep = "")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
