\name{visReps}
\alias{visReps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize where replicates occur in an ionomics run
}
\description{
Takes a vector of when each replicate occurred (e.g. if two samples A and B were run twice in this order A,B,A,B the vector would be 1,1,2,2
if the replicates were run A,A,B,B the vector would be 1,2,1,2)
}
\usage{
visReps(vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
Vector of replicate counts.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns ggplot object
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
df <- data.frame(id=sample(LETTERS,length(x),replace=TRUE),x,y,z)
numReps <- as.data.frame(table(df$id))
for(i in numReps$Var1){
  df$repNum[df$id == i & !is.na(df$id)] <- c(1:numReps$Freq[numReps$Var1 == i])
}  
visReps(df$repNum)


## The function is currently defined as
function (vec) 
{
    require(ggplot2)
    vec <- na.omit(vec)
    df <- data.frame(sampleNum = factor(1:length(vec), ordered = TRUE), 
        repNum = as.factor(vec), one = 1)
    p <- ggplot(df, aes(x = sampleNum, y = one, fill = repNum, 
        colour = repNum, width = 0.5)) + geom_bar(stat = "identity", 
        position = "dodge") + theme(axis.line = element_blank(), 
        axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank())
    print(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
