#!/usr/bin/Rscript
#source("http://bioconductor.org/biocLite.R")
#biocLite("DNAcopy")
#biocLite("optparse")

# load libraries
library("DNAcopy")

# function for the whole CBS.
CBS = function(bins, chromosomes, positions, filtering = 0, sampleid = "unknown", min_width = 5){
    # set random generator seed
    set.seed(5)
    # create CNA object
    CNA_f = CNA(cbind(bins), chromosomes, positions, data.type="logratio", sampleid=sampleid, presorted=TRUE)
    # smooth it
    smoothed.CNA_f = smooth.CNA(CNA_f)

    if (filtering) {
        # do a CBS with filtering of segments
        return(segment(smoothed.CNA_f, min.width=min_width, undo.splits="sdundo", undo.SD=2, verbose=0)$output)
    } else {
        # do a simple CBS
        return(segment(smoothed.CNA_f, min.width=min_width, verbose=0)$output)
    }
}
