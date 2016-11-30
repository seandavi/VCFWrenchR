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
    


setMethod('as.data.frame',
          signature = c('ExpandedVCF'),
          function(x,...) {
              df = as.data.frame(rowRanges(x))
              df = cbind(df,as.data.frame(info(x)))
              dfann = .anncols(df$ANN,info(header(x))['ANN',]$Description)
              df = df[,colnames(df)!="ANN"]
              df = cbind(df,dfann)
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
setMethod('as.data.frame',
          signature = c('CollapsedVCF'),
          function(x,...) {
              message('Expanding VCF first, so number of rows may increase')
              return(as.data.frame(expand(x)))
          })




.longestRow = function(raggedList) {
    return(max(sapply(raggedList,length)))
}

.raggedListToMatrix = function(raggedList) {
    maxlen = max(sapply(raggedList, length))
    ret = do.call(rbind,lapply(raggedList,function(l) {
        length(l) = maxlen
        return(l)
    }))
    return(ret)
}

.dataFrametoExpandedDataFrame = function(df) {
    
    listcols = sapply(df,function(y) {inherits(y,'List')})
    maxlens = lapply(df[,listcols],function(col) {return(.longestRow(col))})
    z = DataFrame(do.call(
        cbind,
        lapply(df[,listcols],function(x1) {
            raggedListToDataFrame(x1)
        })
    ))
    colnames(z) = make.unique(rep(colnames(df)[listcols],maxlens))
    df = df[,!listcols]
    z = cbind(df,z)
    return(z)
}


#' Dump vcf file as csv or txt file
#'
#'
#'
#' @importFrom VariantAnnotation readVcf
#' @export
vcf2df = function(vcf,flatten=TRUE,collapse=":") {
    vr = as(vcf,'VRanges')
    df = mcols(vr)
    df = cbind(as.data.frame(ranges(vr)),df)
    df = cbind(data.frame(chrom = seqnames(vr)),df)
    effcols = strsplit(gsub("\\s+|\\[|\\)|\\]|\\'",'',sub(".*'Effect",'Effect',info(header(vcf))['EFF','Description'])),'\\||\\(')[[1]]
    anncols = strsplit(gsub("\\s+|\\[|\\)|\\]|\\'",'',sub(".*'Effect",'Effect',info(header(vcf))['ANN','Description'])),'\\||\\(')[[1]]
    # This mess of code is all to split the EFF column, then convert to a data.frame with
    # the appropriate types.
    # TAKES the first record ONLY--may not be correct for EFF field
    ## effects = data.frame(lapply(data.frame(
    ##     do.call(rbind,strsplit(sapply(vr$EFF,function(x) {
    ##         return(sub('\\)','',x[[1]]))}),'\\(|\\|'))),
    ##     function(z) {
    ##         return(type.convert(as.character(z)))}
    ##     ))
    ann = data.frame(lapply(data.frame(
        do.call(rbind.fill,strsplit(sapply(vr$ANN,function(x) {
            if(length(x)==0)
                return('')
            else
                return(sub('\\)','',x[[1]]))}),'\\(|\\|'))),
        function(z) {
            return(type.convert(as.character(z)))}
        ))
#    effcols = effcols[seq_len(ncol(effects))]
#    colnames(effects) = effcols
#    df = cbind(df,effects)
    anncols = anncols[seq_len(ncol(ann))]
    colnames(ann) = anncols
    df = cbind(df,ann)
    return(df)
}
