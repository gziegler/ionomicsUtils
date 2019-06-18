\name{calcMAF}
\alias{calcMAF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Takes a recoded (0,1,2) data.table and returns a two column data.table of num missing and minor allele frequency.
}
\description{
Takes a recoded (0,1,2) data.table and returns a two column data.table of num missing and minor allele frequency ((minor alelle)/total)
FracMissing is number of NA values/number of lines
MAF is (number of minor allele + 0.5 the number of het allele)/(number of nonNA alleles)
}
\usage{
calcMAF(genoTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genoTable}{
Data table where SNPs are rows and lines are columns, no metadata columns. Coded as 0,1,2.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a data.table with two columns and rows in the same order as the input. One column contains the MAF and one contains the FracMissing for each row.
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
alleleTable <- cbind(alleleTable,calcMAF(genoTable))
#decide who to delete
alleleTable[,discard := (FracMissing>0.2 | MAF<0.1)]

#get rid of the rows in snp that are TRUE in alleleTable$discard
#memory efficient approach to removing row
# see http://stackoverflow.com/questions/10790204/how-to-delete-a-row-by-reference-in-data-table?lq=1
keepIdxs <- which(alleleTable$discard==FALSE)
cols <- names(GenotypeData)
firstCol <- cols[1]
snp.subset <- data.table(col1 = GenotypeData[[firstCol]][keepIdxs]) 
setnames(snp.subset,"col1",firstCol)
for(col in cols[2:length(cols)]){
  snp.subset[, (col) := GenotypeData[[col]][keepIdxs],]
  GenotypeData[, (col):= NULL,] #delete
}
rm(GenotypeData)
snpInfo <- alleleTable[discard==FALSE]

}
