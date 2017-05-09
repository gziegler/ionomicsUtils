\name{SymboltoElementName}
\alias{SymboltoElementName}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert element symbol to full element name
}
\description{
Takes an element symbol or symbol isotope transformation type (e.g. Na23_normConc) and converts to element name (e.g. Sodium)
If an element isn't found in the lookup table, the value is returned unchanged with a warning.
Translation currently only includes commonly used ionomic elements.
}
\usage{
SymboltoElementName(els)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{els}{
String or vector of values to be converted.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns of vector of same length as input with translated values.
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
SymboltoElementName(c("Na23_normConc","Na23_corrConc","Na23","S","weight","otherID"))

## The function is currently defined as
function (els) 
{
    els <- sub("_\\\w+", "", els)
    els <- sub("\\\d+", "", els)
    elTable <- data.frame(Element = c("Boron", "Sodium", "Magnesium", 
        "Aluminum", "Phosphorus", "Sulfur", "Potassium", "Calcium", 
        "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc", 
        "Arsenic", "Selenium", "Rubidium", "Strontium", "Molybdenum", 
        "Cadmium", "Indium", "Yttrium", "Lead"), Symbol = c("B", 
        "Na", "Mg", "Al", "P", "S", "K", "Ca", "Mn", "Fe", "Co", 
        "Ni", "Cu", "Zn", "As", "Se", "Rb", "Sr", "Mo", "Cd", 
        "In", "Y", "Pb"), stringsAsFactors = FALSE)
    out <- elTable[match(els, elTable$Symbol), "Element"]
    if (length(which(is.na(out))) > length(which(is.na(els)))) {
        warning("Some values from input not found in symbol lookup table. Original value returned")
    }
    out[which(is.na(out))] <- els[which(is.na(out))]
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
