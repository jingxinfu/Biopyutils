#!/usr/bin/Rscript
# License           : GPL3
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 11/02/2020
# Last Modified Date: 20/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

main = function(){
    args <- parse_input() 
    exprsn <- fread(args$exprsn_path,data.table = F)
    rownames(exprsn) <- exprsn[,1]
    exprsn <- as.matrix(exprsn[,-1])
    geneSets <- customGeneSet(rownames(exprsn),args$geneset_path)
    result <- t(gsva(exprsn,geneSets,method=args$method,kcdf=args$kcdf,verbose=F))
    write.csv(result,args$output)
}

customGeneSet = function(candiates,gset_path){
    gset = read.table(gset_path,row.names=1,header=F,stringsAsFactor=F,sep='\t')
    filter_gset = lapply(rownames(gset),function(n){
        member = strsplit(gset[n,1],split=',')[[1]]
        member = member[member %in% candiates]
        if(length(member)>0){
            return(member)
        }else{
            return(NA)
        }
    })

    names(filter_gset) = rownames(gset)
    return(filter_gset[!is.na(filter_gset)])
}

parse_input = function(){
    require_pkgs <- c('data.table','GSVA','argparse')
    scripts_wd <- dirname(thisFile())
    source(file.path(scripts_wd,'utils.R'))
    load_package(require_pkgs)

    parser <- ArgumentParser(description='Execute R bioconductors GSVA')
    inputs <- parser$add_argument_group('Input Option')
    inputs$add_argument('--exprsn_path',help='Exrpsn Path')
    inputs$add_argument('--geneset_path',required=T,help='Gene set Path')

    outputs <- parser$add_argument_group('Output Option')
    outputs$add_argument('--output',help='Output Path')

    cmds <- parser$add_argument_group('GSVA command Option')
    kcdf_str <-'Character string denoting the kernel to use during the non-parametric estimation
    of the cumulative distribution function of expression levels across samples
    when method="gsva". By default, kcdf="Gaussian" which is suitable when
    input expression values are continuous, such as microarray fluorescent units in
    logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input
    expression values are integer counts, such as those derived from RNA-seq experiments,
    then this argument should be set to kcdf="Poisson". This argument
    supersedes arguments rnaseq and kernel, which are deprecated and will be
    removed in the next release.
    '
    cmds$add_argument('--kcdf',choices=c("Gaussian", "Poisson", "none"),default='Gaussian',help=gsub('\n',' ',kcdf_str))
    method_str <-'Method to employ in the estimation of gene-set enrichment scores per sample.
    By default this is set to gsva (Hanzelmann et al, 2013) and other options
    are ssgsea (Barbie et al, 2009), zscore (Lee et al, 2008) or plage (Tomfohr
    et al, 2005). The latter two standardize first expression profiles into z-scores
    over the samples and, in the case of zscore, it combines them together as their
    sum divided by the square-root of the size of the gene set, while in the case of
    plage they are used to calculate the singular value decomposition (SVD) over
    the genes in the gene set and use the coefficients of the first right-singular vector
    as pathway activity profile.'
    cmds$add_argument('--method',choices=c("gsva", "ssgsea", "zscore", "plage"),default='gsva',help=gsub('\n',' ',method_str))
    args <- parser$parse_args()
    return(args)
}

thisFile = function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  opt_input = cmdArgs[-match]
  if (length(match) > 0) {
    # Rscript
    return(sub(needle, "", cmdArgs[match]))
  } else {
    # 'source'd via R console
    return(sys.frames()[[1]]$ofile)
  }
}


if (!interactive()){
    #main()
    suppressWarnings(suppressMessages(main()))
}
