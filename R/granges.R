selectFeaturesGR <- function (gr, type="start", width=10, feature="gene", randomranges=NULL) {
  #' @title selectFeaturesGR
  #' @description Returns new GRanges object based on the selection criteria.
  #'
  #' @param type Could be "start", "end" or "random"
  #' @param gr \code{GRanges} object - it is assumed that it is full genomic annotation imported from GFF
  #' @param width Integer indicating length of selected ranges
  #' @param feature "gene" or "exon", to be selected from GFF annotation. "gene" is the default.
  #' @param randomranges optional argumument, required if \code{type="random"} is supplied; this is a
  #' vector of widths
  #' @import GenomicRanges
  #'
  #' @return GenomicRange object containing only selected features
  #' @export selectFeaturesGR
  #' #@examples selectFeaturesGR(GRfromGFF, type="end", width=25, feature="exon")
  #'
  #' #using random ranges
  #'
  #' #selectFeaturesGR(GRfromGFF, type="random", randomranges=c(89,34,1,4,25))

  stopifnot(verifySeqLevels(gr))
  out<-gr[which(elementMetadata(gr)$type == feature)]

  if (type == "start" | type == "end") {
    out <- resize(out, fix=type, width=width)

  }

  if (type == "random") {
    if (is.null(randomranges)){
      stop("Fill required parameter for random type: randomranges (a vector of widths)")
    }
    out <- randomGenomicRanges(randomranges, gr)
  }
  return(out)

}


randomGenomicRanges <- function(x, gr) {
  #' @title randomGenomicRanges
  #' @description Generate random Genomic Ranges given \code{GRanges} object. The resulting
  #' object will have unique starts, but ranges (depending on the vector of lengths)
  #' might have overlaps.
  #'
  #' @param x vector of widths
  #' @param gr \code{GRanges} object
  #' @import BSGenome
  #' @import GenomicRanges
  #'
  #' @return GenomicRange object containing randomly selected uniquue ranges of widths requested by vector \code{x}
  #' @export randomGenomicRanges
  #' #@examples randomGenomicRanges(c(4,89,325), GRobject)
  #'
  #'


  n <- length(x)
  names <- as.vector(seqnames(gr))
  available_ranges <- unlist(mapply(seq, from=start(gr), to=end(gr)))
  out <- GRanges(sample(names, n, replace=TRUE), IRanges(sample(available_ranges, n), width=x) )
  return(out)

}
