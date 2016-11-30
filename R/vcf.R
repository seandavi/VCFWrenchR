.anncols = function(anncol,headerstring) {
    anncols = strsplit(sub("Functional annotations: '",'',
                           headerstring),' \\| ')[[1]]
    dfannempty = data.frame(matrix(vector(), 0, length(anncols),
                                   dimnames=list(c(), anncols)),
                            stringsAsFactors=F)
    yy = data.frame(suppressWarnings(do.call(rbind,
                            c(dfannempty,lapply(lapply(anncol,`[`,1),
                                                function(x){strsplit(x,'\\|')[[1]]})))),
                    stringsAsFactors=FALSE)
    yy = data.frame(lapply(yy,type.convert))
    colnames(yy) = paste("ANN",anncols,sep="_")
    return(yy)
    }
    

#' Convert ExpandedVCF to data.frame
#'
#' This method converts an \linkS4class{ExpandedVCF} to a data.frame.
#'
#' @importClassesFrom VariantAnnotation ExpandedVCF
#'
#' @examples
#' require(VariantAnnotation)
#' vcffile = system.file(package='VariantAnnotation',path='extdata/chr22.vcf.gz')
#' v = readVcf(vcffile,'hg19')
#' ve = expand(v)
#' class(ve)
#' head(as.data.frame(v))
#' 
#' @export
setMethod('as.data.frame',
          signature = c('ExpandedVCF'),
          function(x,...) {
              df = as.data.frame(rowRanges(x))
              df = cbind(df,as.data.frame(info(x)))
              if('ANN' %in% names(info(x))) {
                  dfann = .anncols(df$ANN,info(header(x))['ANN',]$Description)
                  df = df[,colnames(df)!="ANN"]
                  df = cbind(df,dfann)
              }
              n  = names(geno(x))
              tmp = lapply(n,function(col) {
                  return(as.data.frame(geno(x)[[col]]))
              })
              ncols = sapply(tmp,ncol)
              tmp = do.call(cbind,tmp)
              colnames(tmp) = paste(rep(n,times=ncols),colnames(tmp),sep="_")
              df = cbind(df,tmp)
              return(df)
          }
          )

#' Convert CollapsedVCF to data.frame
#'
#' This method converts an \linkS4class{CollapsedVCF} to a data.frame.
#'
#' @importClassesFrom VariantAnnotation CollapsedVCF
#'
#' @examples
#' require(VariantAnnotation)
#' vcffile = system.file(package='VariantAnnotation',path='extdata/chr22.vcf.gz')
#' v = readVcf(vcffile,'hg19')
#' head(as.data.frame(v))
#' 
#' @export
setMethod('as.data.frame',
          signature = c('CollapsedVCF'),
          function(x,...) {
              message('Expanding VCF first, so number of rows may increase')
              return(as.data.frame(expand(x)))
          })



setGeneric('as.json',
           function(x,...) {
               standardGeneric('as.json')
           })

#' Convert VCF object to json
#'
#' This method converts an \linkS4class{VCF} to a JSON string.
#'
#' @importFrom jsonlite toJSON
#' @importClassesFrom VariantAnnotation VCF
#'
#' @examples
#' require(VariantAnnotation)
#' vcffile = system.file(package='VariantAnnotation',path='extdata/chr22.vcf.gz')
#' v = readVcf(vcffile,'hg19')
#' as.json(head(v,3),pretty=TRUE)
#' 
#' @export
setMethod('as.json',signature=c('VCF'),
          function(x,...) {
              return(jsonlite::toJSON(as.data.frame(x),...))
          })
