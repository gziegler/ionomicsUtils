\name{recodeBiallele}
\alias{recodeBiallele}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Takes a biallelic SNP vector and converts to 0,1,2 coding
}
\description{
  Looks for AA,GG,CC,TT sets major to 0 (most frequently found basepair), minor to 2 (second most frequently found bp), het to 1 (e.g. AG, AC, AT, ...), anything else to NA

}
\usage{
  recodeBiallele(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Vector of values to be converted.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a vector of same length as input with translated values.
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
  recodeBiallele(c("AA","GG","GG","GG","GG","AA"))
  
  ## The function is currently defined as
function(x){
  x <- as.character(x)
  freqs <- names(sort(table(x[x %in% c("AA","GG","CC","TT")]),decreasing=TRUE))
  major <- freqs[1]
  if(length(freqs)>1){
    minor <- freqs[2]
  }else{
    minor <- "np"
  }
  
  x[which(x==major)] <-0
  x[which(x==minor)] <-2
  x[which(x %in% c("AG","AC","AT","GA","GC","GT","CA","CG","CT","TA","TG","TC"))] <- 1
  x[which(!(x==1|x==0|x==2))]   <- NA
  #x[which(!(x==2|x==0))]   <- NA
  #x[which(x=="N")] <- NA
  return(x)
}

}