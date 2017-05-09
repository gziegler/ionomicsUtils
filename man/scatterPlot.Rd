\name{scatterPlot}
\alias{scatterPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make an ionomics scatter plot with points in order run on x-axis and values on y-axis.
}
\description{
Takes a vector of ionomics data, plots values colored by RSD on scatter plot, with control values as triangles. RSDs are broken into categories. <2,2-4,4-6,6-8,8-10, >10.
}
\usage{
scatterPlot(scatter, rsd = NA, shape = NA, main, xlab = "Sample No.", ylab = "Concentration", runBreaks = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scatter}{
Vector of points for y-axis
}
  \item{rsd}{
Vector of RSD values in same order as scatter column.
}
  \item{shape}{
Column specifying whether each data point is a control or a sample point.
}
  \item{main}{
String for main title of the plot
}
  \item{xlab}{
String for the plot x label.
}
  \item{ylab}{
String for the plot y label.
}
  \item{runBreaks}{
Vector of indices to drow vertical lines. Usually test tube or run breaks.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a ggplot object.
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
library(ggplot2)
set.seed(1)
x <- rnorm(100)
y <- rnorm(100)
x <- c(-10, x, 10)
y <- c(-20, y, 20)
df <- data.frame(id=sample(LETTERS[1:5],length(x),replace=TRUE),x,y)
df$type <- sapply(1:nrow(df),function(x){ifelse(x\%\%10,"Sample","Control")})
scatterPlot(scatter=df$y,rsd=abs(df$x*10),shape=df$type,main="Test data",ylab="CPS",runBreaks=NA)

## The function is currently defined as
function (scatter, rsd = NA, shape = NA, main, xlab = "Sample No.", 
    ylab = "Concentration", runBreaks = NA) 
{
    suppressPackageStartupMessages(require(ggplot2))
    suppressPackageStartupMessages(require(grid))
    data <- data.frame(scatter = scatter, rsd = rsd, shape = shape)
    smooth <- stats::lowess(na.omit(data[, 1]), f = 0.1)$y
    if (any(is.na(data$scatter))) {
        for (na.val in which(is.na(data$scatter))) {
            smooth <- append(smooth, NA, after = na.val - 1)
        }
    }
    data$smooth <- smooth
    if (!sum(is.na(data$rsd)) == length(data$rsd)) {
        data$colorCol <- cut(as.numeric(data$rsd), breaks = c(2 * 
            (0:5), Inf), labels = c("<2", "2-4", "4-6", "6-8", 
            "8-10", ">10"), include.lowest = TRUE)
    }
    data$x <- 1:nrow(data)
    p1 <- ggplot(data = data, aes(x = x, y = scatter))
    if (!sum(is.na(data$rsd)) == length(data$rsd) & !sum(is.na(data$shape)) == 
        length(data$shape)) {
        p1 <- p1 + geom_point(aes(colour = colorCol, shape = shape))
        p1 <- p1 + scale_colour_manual(values = colorRampPalette(c("#00007F", 
            "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
            "red", "#7F0000"))(7), name = "RSD")
        p1 <- p1 + scale_shape_manual(name = "Type", values = c(17, 
            16))
    }
    else if (!sum(is.na(data$rsd)) == length(data$rsd)) {
        p1 <- p1 + geom_point(aes(colour = colorCol))
        p1 <- p1 + scale_colour_manual(values = colorRampPalette(c("#00007F", 
            "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
            "red", "#7F0000"))(7), name = "RSD")
    }
    else if (!sum(is.na(data$shape)) == length(data$shape)) {
        p1 <- p1 + geom_point(aes(shape = shape))
        p1 <- p1 + scale_shape_manual(name = "Type", values = c(17, 
            16))
    }
    else {
        p1 <- p1 + geom_point()
    }
    if (nrow(data) > 50) {
        p1 <- p1 + geom_line(aes(y = smooth), size = 1.3, colour = "blue")
    }
    p1 <- p1 + labs(x = xlab, y = ylab, title = main)
    p1 <- p1 + theme(legend.background = element_rect(), legend.position = "bottom", 
        legend.box = "horizontal", legend.direction = "horizontal")
    p1 <- p1 + scale_x_continuous(limits = c(0, nrow(data)), 
        expand = c(0, 0))
    if (!is.na(runBreaks)) {
        p1 <- p1 + geom_vline(xintercept = runBreaks)
    }
    p1 <- p1 + geom_hline(yintercept = mean(data[, 1], na.rm = TRUE), 
        colour = "red", size = 1)
    return(p1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
