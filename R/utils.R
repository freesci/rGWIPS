
loadGWIPSdata <- function(forward, reverse){
  #' @title loadGWIPSdata
  #' @description Loads data from GWIPS database (in BigWig format) for forward and reverse strands into global GRanges objects. As a result, global objects "gwips_forw" and "gwips_rev" become available for further analyses.
  #'
  #' @param forward path to the BigWig file containing Riboseq data for forward strand
  #' @param reverse path to the BigWig file containing Riboseq data for forward strand
  #' @import rtracklayer
  #'
  #' @return NULL
  #' @export loadGWIPSdata


  stopifnot( file.exists( forward ) )
  stopifnot( file.exists( reverse ) )
  gwips_forw <<- rtracklayer::import(con=forward, format="bigWig")
  gwips_rev <<- rtracklayer::import(con=reverse, format="bigWig")

}


verifySeqLevels <- function(x) {
  #' @title verifySeqLevels
  #' @description Assumes that global variables gwips_forw and gwips_rev are available. This function checks if the seqLevels in all Genomic Ranges objects are the same. It's a sanity-check function, not to be used directly.
  #'
  #' @param x Genomic Ranges object from the genome annotation
  #' @import GenomicRanges
  #'
  #' @return 0 if seqLevels are different, 1 if everything is OK
  #' @export verifySeqLevels

  if (! ( (seqlevels(gwips_forw) == seqlevels(gwips_rev)) && (seqlevels(gwips_rev) == seqlevels(x)) ) ){
    print("Please correct seqlevels in your Granges objects.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}
