\name{multiplot}
\alias{multiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multiple plot function
}
\description{
Put multiple plots on the same page.

}
\usage{
multiplot(..., plotlist = NULL, file, cols = 1, layout = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
List of ggplot objects to put on page
}
  \item{plotlist}{
Optional plotlist (list of ggplots) to put on page
}
  \item{file}{
Not implemented.
}
  \item{cols}{
Number of columns in layout.
}
  \item{layout}{
A matrix specifying the layout. If present, 'cols' is ignored.
}
}
\details{
If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
then plot 1 will go in the upper left, 2 will go in the upper right, and
3 will go all the way across the bottom.
}
\value{
No return value. Prints a plot to graphic device.
}
\references{
http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
}
\author{
Online
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(ggplot2)
set.seed(1)
x <- rnorm(100)
y <- rnorm(100)
x <- c(-10, x, 10)
y <- c(-20, y, 20)
df <- data.frame(id=sample(LETTERS[1:5],length(x),replace=TRUE),x,y)
p1 <- ggplot(df,aes(x=x,y=y)) + geom_point()
p2 <- ggplot(df,aes(x=x,y=y)) + geom_point() + geom_smooth(method="lm")
p3 <- ggplot(df,aes(x=id,y=y)) + geom_boxplot()
p4 <- ggplot(df,aes(x=id,y=y)) + geom_boxplot() + geom_jitter()
multiplot(p1, p2, p3, p4, layout=matrix(c(1,3,2,4),ncol=2))


## The function is currently defined as
function (..., plotlist = NULL, file, cols = 1, layout = NULL) 
{
    require(grid)
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), 
            ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots == 1) {
        print(plots[[1]])
    }
    else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), 
            ncol(layout))))
        for (i in 1:numPlots) {
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, 
                layout.pos.col = matchidx$col))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
