\name{TravisPlot}
\alias{TravisPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TravisPlot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TravisPlot(peak, TravisCoordsFromTxDb = NA, txdb = NA, genome = NA, noBins = 10, saveToPDFprefix = NA, returnCount = FALSE, includeNeighborDNA = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{peak}{
%%     ~~Describe \code{peak} here~~
}
  \item{TravisCoordsFromTxDb}{
%%     ~~Describe \code{TravisCoordsFromTxDb} here~~
}
  \item{txdb}{
%%     ~~Describe \code{txdb} here~~
}
  \item{genome}{
%%     ~~Describe \code{genome} here~~
}
  \item{noBins}{
%%     ~~Describe \code{noBins} here~~
}
  \item{saveToPDFprefix}{
%%     ~~Describe \code{saveToPDFprefix} here~~
}
  \item{returnCount}{
%%     ~~Describe \code{returnCount} here~~
}
  \item{includeNeighborDNA}{
%%     ~~Describe \code{includeNeighborDNA} here~~
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
function (peak, TravisCoordsFromTxDb = NA, txdb = NA, genome = NA, 
    noBins = 10, saveToPDFprefix = NA, returnCount = FALSE, includeNeighborDNA = FALSE) 
{
    suppressWarnings(if (is.na(TravisCoordsFromTxDb) & is.na(txdb) & 
        is.na(genome)) {
        stop("Most provide one of the three: TravisCoords, txdb or genome")
    })
    if (suppressWarnings(is.na(TravisCoordsFromTxDb))) {
        if (suppressWarnings(is.na(txdb))) {
            print("Downloading Transcriptome Information from UCSC ...")
            txdb <- suppressMessages(makeTxDbFromUCSC(genome = genome))
            print("Making Travis Coordinates ...")
            TravisCoords <- suppressMessages(makeTravisCoordsFromTxDb(txdb))
        }
        else {
            print("Making Travis Coordinates from provided TranscriptDb Object ...")
            TravisCoords <- makeTravisCoordsFromTxDb(txdb, noBins = noBins)
        }
    }
    else {
        print("Using provided Travis Coordinates")
        TravisCoords <- TravisCoordsFromTxDb
    }
    noGroup <- length(peak)
    group_names <- names(peak)
    m <- peak
    if (is.null(group_names)) {
        group_names <- paste("item", 1:noGroup)
    }
    for (i in 1:noGroup) {
        temp = .countTravisDensity(peak[[i]], TravisCoords)
        temp = cbind(temp, Feature = group_names[i])
        m[[i]] = temp
    }
    ct = .combineListOfDataFrame(m)
    ct[[4]] <- as.character(ct[[4]])
    count_result <- ct
    temp <- ct$category
    ct <- ct[ct$count > 0, ]
    ct[ct$category == "lncRNA", 2] <- ct[ct$category == "lncRNA", 
        2] + 6
    ct <- data.frame(ct, weight = ct$count)
    ct$weight <- as.numeric(ct$weight)
    ct1 <- ct[ct$category == "mRNA", ]
    ct2 <- ct[ct$category == "lncRNA", ]
    featureSet <- as.character(unique(ct$Feature))
    for (i in 1:length(featureSet)) {
        id <- (ct1$Feature == featureSet[i])
        ct1$weight[id] <- ct1$weight[id]/sum(ct1$weight[id])
        id <- (ct2$Feature == featureSet[i])
        ct2$weight[id] <- ct2$weight[id]/sum(ct2$weight[id])
    }
    if (includeNeighborDNA) {
        p1 <- ggplot(ct1, aes(x = pos, group = Feature, weight = 5 * 
            weight)) + scale_x_continuous(minor_breaks = seq(0, 
            5, 1)) + ggtitle("Distribution on mRNA") + theme(axis.ticks = element_blank(), 
            axis.text.x = element_blank()) + xlab("") + ylab("Frequency") + 
            annotate("rect", xmin = 1, xmax = 4, ymin = -0.9, 
                ymax = -0.5, alpha = 0.2, fill = "red") + annotate("rect", 
            xmin = 0, xmax = 1, ymin = -0.9, ymax = -0.1, alpha = 0.2, 
            fill = "green") + annotate("rect", xmin = 1, xmax = 2, 
            ymin = -0.5, ymax = -0.1, alpha = 0.2, fill = "yellow") + 
            annotate("rect", xmin = 2, xmax = 3, ymin = -0.5, 
                ymax = -0.1, alpha = 0.2, fill = "orange") + 
            annotate("rect", xmin = 3, xmax = 4, ymin = -0.5, 
                ymax = -0.1, alpha = 0.2, fill = "yellow") + 
            annotate("rect", xmin = 4, xmax = 5, ymin = -0.9, 
                ymax = -0.1, alpha = 0.2, fill = "green") + geom_density(adjust = 0.5, 
            aes(fill = factor(Feature)), alpha = 0.2) + annotate("text", 
            x = 2.5, y = -0.7, label = "mRNA") + annotate("text", 
            x = 1.5, y = -0.3, label = "5'UTR", size = 3) + annotate("text", 
            x = 2.5, y = -0.3, label = "CDS", size = 3) + annotate("text", 
            x = 0.5, y = -0.5, label = "Promoter") + annotate("text", 
            x = 4.5, y = -0.5, label = "Tail") + annotate("text", 
            x = 3.5, y = -0.3, label = "3'UTR", size = 3)
        p2 <- ggplot(ct2, aes(x = pos, group = Feature, weight = 3 * 
            weight)) + ggtitle("Distribution on lncRNA") + scale_x_continuous(minor_breaks = seq(6, 
            9, 1)) + theme(axis.ticks = element_blank(), axis.text.x = element_blank()) + 
            xlab("") + ylab("Frequency") + annotate("rect", xmin = 6, 
            xmax = 7, ymin = -0.9, ymax = -0.1, alpha = 0.2, 
            fill = "green") + annotate("rect", xmin = 7, xmax = 8, 
            ymin = -0.9, ymax = -0.1, alpha = 0.2, fill = "red") + 
            annotate("rect", xmin = 8, xmax = 9, ymin = -0.9, 
                ymax = -0.1, alpha = 0.2, fill = "green") + geom_density(adjust = 0.5, 
            aes(fill = factor(Feature)), alpha = 0.2) + annotate("text", 
            x = 7.5, y = -0.5, label = "lncRNA") + annotate("text", 
            x = 6.5, y = -0.5, label = "Promoter") + annotate("text", 
            x = 8.5, y = -0.5, label = "Tail")
    }
    else {
        p1 <- ggplot(ct1, aes(x = pos, group = Feature, weight = 5 * 
            weight)) + ggtitle("Distribution on mRNA") + theme(axis.ticks = element_blank(), 
            axis.text.x = element_blank()) + xlab("") + ylab("Frequency") + 
            annotate("rect", xmin = 1, xmax = 4, ymin = -0.9, 
                ymax = -0.5, alpha = 0.2, fill = "red") + annotate("rect", 
            xmin = 1, xmax = 2, ymin = -0.5, ymax = -0.1, alpha = 0.2, 
            fill = "yellow") + annotate("rect", xmin = 2, xmax = 3, 
            ymin = -0.5, ymax = -0.1, alpha = 0.2, fill = "orange") + 
            annotate("rect", xmin = 3, xmax = 4, ymin = -0.5, 
                ymax = -0.1, alpha = 0.2, fill = "yellow") + 
            geom_density(adjust = 0.5, aes(fill = factor(Feature)), 
                alpha = 0.2) + annotate("text", x = 2.5, y = -0.7, 
            label = "mRNA") + annotate("text", x = 1.5, y = -0.3, 
            label = "5'UTR", size = 3) + annotate("text", x = 2.5, 
            y = -0.3, label = "CDS", size = 3) + annotate("text", 
            x = 3.5, y = -0.3, label = "3'UTR", size = 3) + xlim(1, 
            4)
        p2 <- ggplot(ct2, aes(x = pos, group = Feature, weight = 3 * 
            weight)) + ggtitle("Distribution on lncRNA") + theme(axis.ticks = element_blank(), 
            axis.text.x = element_blank()) + xlab("") + ylab("Frequency") + 
            annotate("rect", xmin = 7, xmax = 8, ymin = -0.9, 
                ymax = -0.1, alpha = 0.2, fill = "red") + geom_density(adjust = 0.5, 
            aes(fill = factor(Feature)), alpha = 0.2) + annotate("text", 
            x = 7.5, y = -0.5, label = "lncRNA") + xlim(7, 8)
    }
    suppressWarnings(if (is.na(saveToPDFprefix)) {
        print("no figure saved ...")
    }
    else {
        f1 <- paste(saveToPDFprefix, "_Travis.pdf", sep = "")
        pdf(file = f1, width = 9, height = 9)
        .multiplot(p1, p2, cols = 1)
        dev.off()
        print(paste("Figures saved into", f1, "...", sep = " "))
    })
    multiplot(p1, p2, cols = 1)
    if (returnCount) {
        q <- list(original = ct, mRNA_normalized = ct1, ncRNA_normalized = ct2)
        return(q)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
