#' Convert structural variant VCF to GRanges
#'
#' This function is just an opinionated approach to convert a VCF file
#' that encodes SVs into a GRanges object. The rationale for this
#' function is to be identify locations in the genome that disrupt
#' genes (or other regions of interest). An overlap of the returned
#' regions from this function with genes will result in a list of
#' putatively disrupted genes.
#'
#'
#' BND (break ends) are converted to single base pair
#' \code{GRanges}. DUPs and INVs are considered as having two breaks and are
#' converted to two individual single base pair regions. DELs are
#' considered as contiguous regions and are converted to one GRanges
#' object with \code{start} and \code{end} as specified in the VCF
#' file.
#'
#' @importFrom GenomicRanges GRanges
#'
#' @param vcf \code{\link[VariantAnnotation]{VCF}} object encoding
#'     variants in VCF format (see
#'     \link{\url{http://samtools.github.io/hts-specs/VCFv4.2.pdf}}
#' @param max_deletion_size integer(1) giving the upper limit on
#'     size of deletions to include. The idea is to treat focal
#'     deletions a little differently than large-scale deletions.
#'
#' @return \code{\link[GenomicRanges]{GRanges}} object with the
#'     structural variants coded as in Details
#' 
#' @export
structuralVariantVCFToGRanges = function(vcf, max_deletion_size = Inf) {
    if(!inherits(vcf,'VCF')) {
        stop('vcf parameter must inherit from VariantAnnotation::VCF class')
    }
    bndvars = vcf[info(vcf)$SVTYPE=="BND"]
    delvars = vcf[info(vcf)$SVTYPE=="DEL"]
    dupvars = vcf[info(vcf)$SVTYPE=="DUP"]
    invvars = vcf[info(vcf)$SVTYPE=="INV"]
    
    bnds = GRanges(seqnames(bndvars),IRanges(start(bndvars),width=1))
    mcols(bnds) = mcols(bndvars)
    bnds = bnds[width(bnds)<= max_deletion_size]
    dups = c(GRanges(seqnames(dupvars),
                     ranges = IRanges(start(dupvars),width=1)),
             GRanges(seqnames=seqnames(dupvars),
                     ranges=IRanges(start=info(dupvars)$END,width=1)))
    mcols(dups) = rbind(mcols(dupvars),mcols(dupvars))
    invs = c(GRanges(seqnames(invvars),
                     ranges = IRanges(start(invvars),width=1)),
             GRanges(seqnames=seqnames(invvars),
                     ranges=IRanges(start=info(invvars)$END,width=1)))
    mcols(invs) = rbind(mcols(invvars),mcols(invvars))
    dels = GRanges(seqnames=seqnames(delvars),
                   ranges=IRanges(start=start(delvars),end=info(delvars)$END))
    mcols(dels) = mcols(delvars)
    dels = dels[width(dels) <= max_deletion_size]

    c(bnds, dups, dels, invs)
}
