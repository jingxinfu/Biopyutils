#!/usr/bin/Rscript
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 14/02/2020
# Last Modified Date: 19/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

main = function(){
    args <- parse_input()
    exprsn <- fread(args$exprsn_path,data.table=F)
    rownames(exprsn) = exprsn[,1]
    exprsn <- as.matrix(exprsn[,-1])
    # Remove Zero Variance gene
    exprsn <- exprsn[apply(exprsn, 1, var) > 0,]
    meta <- read.csv(args$meta_path,header=T,row.names=1)
    meta <- meta[colnames(exprsn),]
    modcombat = model.matrix(~1, data=meta)
    result <- ComBat(exprsn,batch=meta$Batch,mod=modcombat,par.prior=T,prior.plots=F)
    if(!is.null(args$fig_out)){
        pdf(file=args$fig_out)
        par(mfrow=c(1,2))
        boxplot(as.data.frame(exprsn),main="Original")
        boxplot(as.data.frame(result),main="Batch corrected")
        dev.off()
    }
    write.csv(result,paste0(args$output))
}

parse_input = function(){
    require_pkgs <- c('data.table','sva','argparse')
    scripts_wd <- dirname(thisFile())
    source(file.path(scripts_wd,'utils.R'))
    load_package(require_pkgs)

    parser <- ArgumentParser(description='Execute R bioconductors limma')
    inputs <- parser$add_argument_group('Input Option')
    inputs$add_argument('--exprsn_path',help='Exprsn Path')
    inputs$add_argument('--meta_path',help='Sample info Path')
    outputs <- parser$add_argument_group('Output Option')
    outputs$add_argument('--output',help='Output Path')
    outputs$add_argument('--fig_out',help='Figure output Path')
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
if(!interactive()){
    #main()
    suppressWarnings(suppressMessages(main()))
}
