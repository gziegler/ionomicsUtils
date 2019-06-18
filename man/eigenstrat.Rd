\name{calcHet}
\alias{calcHet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate PCs similar to eigenstrat of genotype matrix
}
\description{
Takes a recoded (0,1,2,NA) matrix where SNPs are rows and a vector of minor allele frequencies for each row. Returns object with PCs, mu, and sigma values. Can handle NAs, but simply removes any rows containing NA, so NA should be imputed first.
}
\usage{
eigenstrat(W=genoMat,maf=snpInfo$MAF)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{
Matrix where SNPs are rows and lines are columns, no metadata columns. Coded as 0,1,2,NA. NA rows will be removed, so most likely should be imputed first.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an eigenstrat object
}
\references{
This was modified from http://www.popgen.dk/software/index.php/Rscripts
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
#impute all missing to 0
genoMat <- as.matrix(genotype)
class(genoMat) <- "numeric"
genoMat[which(is.na(genoMat))] <- 0

####Genotype file should be in 0, 1, 2 coding
###rows are snps and columns are individuals

e<-eigenstrat(genoMat)
pdf("../results/5.Eigenstrat.structure.PC1vsPC2.and.VarianceExplained.pdf")
plot(e,xlab="PC1",ylab="PC2")
barplot(height=e$values[1:10]/sum(e$values))
dev.off()
structData <- e$vectors[,1:10]
rownames(structData) <- colnames(genotype)
varExplained <- e$values[1:10]/sum(e$values)

}
