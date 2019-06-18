\name{calcHet}
\alias{calcHet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Astle-Balding Kinship using Synbreed Kinship Function
}
\description{
Takes a recoded (0,1,2,NA) matrix where SNPs are columns and a vector of minor allele frequencies for each column. Returns nxn kinship matrix
}
\usage{
realizedAB(W=genoMat,maf=snpInfo$MAF)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{
Matrix where SNPs are columns and lines are rows, no metadata columns. Coded as 0,1,2,NA.
}
  \item{maf}{
  Previously calculated minor allele frequencies for each row.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an nxn kinship matrix
}
\references{
See R package synbreed kinship calculations for original source.
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
genoMat <- as.matrix(genotype)
class(genoMat) <- "numeric"
genoMat <- t(genoMat)
colnames(genoMat) <- snpInfo$`rs#`
kinship <- realizedAB(genoMat,maf=snpInfo$MAF)

pdf(paste("../results/5.AstleBalding.synbreed.kinship.pdf",sep=""), width = 12, height = 12)
par(mar = c(25,25,25,25))
heatmap.2(kinship,  cexRow =.2, cexCol = 0.2, col=rev(heat.colors(256)), scale="none", symkey=FALSE, trace="none")
dev.off()

}

