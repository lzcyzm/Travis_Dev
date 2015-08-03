\name{getNeighborhood}
\alias{getNeighborhood}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getNeighborhood(comp, side = 5, binWidth = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comp}{
%%     ~~Describe \code{comp} here~~
}
  \item{side}{
%%     ~~Describe \code{side} here~~
}
  \item{binWidth}{
%%     ~~Describe \code{binWidth} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (comp, side = 5, binWidth = 5000) 
{
    tx_name <- names(comp)
    granges <- unlist(comp)
    tx <- granges[tx_name]
    strand <- as.character(as.data.frame(strand(tx))[[1]])
    chr <- as.character(as.data.frame(seqnames(tx))[[1]])
    if (side == 5) {
        checkpoints <- rep(1, length(tx))
        checkpoints_transcript <- GRanges(seqnames = chr, IRanges(start = checkpoints, 
            end = checkpoints, names = tx_name), strand = strand)
        checkPoints_genomic <- pmapFromTranscripts(checkpoints_transcript, 
            comp)
        checkRegion_genomic <- resize(x = checkPoints_genomic, 
            width = binWidth, fix = "end")
    }
    else if (side == 3) {
        checkpoints <- sum(width(comp))
        checkpoints_transcript <- GRanges(seqnames = chr, IRanges(start = checkpoints, 
            end = checkpoints, names = tx_name), strand = strand)
        checkPoints_genomic <- pmapFromTranscripts(checkpoints_transcript, 
            comp)
        checkRegion_genomic <- resize(x = checkPoints_genomic, 
            width = binWidth, fix = "start")
    }
    names(checkRegion_genomic) <- rep("", length(tx))
    sidelist <- split(checkRegion_genomic, tx_name, drop = TRUE)
    sidelist <- sidelist[tx_name]
    mapped_chr <- as.character(as.data.frame(seqnames(checkRegion_genomic))[[1]])
    mcols(sidelist) <- data.frame(mapped_chr)
    return(sidelist)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
