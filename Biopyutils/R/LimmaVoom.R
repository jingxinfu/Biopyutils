#!/usr/bin/Rscript
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 14/02/2020
# Last Modified Date: 27/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

main = function(){
    args <- parse_input()
    counts <- fread(args$counts_path,data.table=F)
    rownames(counts) = counts[,1]
    counts <- as.matrix(counts[,-1])
    meta <- read.csv(args$meta_path,row.names=1)
    ## Align samples 
    meta <- meta[colnames(counts),]
    group <- meta[,args$group_name]
    contrast <- readLines(args$contrast_path)
    result <- LimmaVoom(counts=counts,group=group,contrast=contrast,batch=meta$Batch,fig_out=args$fig_out)

    for(compare in contrast){
        tmp = result[[compare]]
        write.csv(tmp,paste0(args$output,compare))
    }
}

LimmaVoom = function(counts,group,contrast,batch=NULL,cutoff=1,fig_out=NULL){
    d0 <- DGEList(counts)
    # Calculate normalization factors
    d0 <- calcNormFactors(d0)
    # Filter low-expressed genes
    drop <- which(apply(cpm(d0), 1, max) < cutoff)
    d <- d0[-drop,]
    # Construct the model wihle remove batch effect
    mm <- model.matrix(~ 0+group+batch)
    plot = F
    if(!is.null(fig_out)){
        plot = T
        pdf(file=paste0(fig_out,'.pdf'))
        # Check whether the batch effect has been remove or not
        if(!is.null(batch))plotMDS(d, col = as.numeric(batch))
    }
    y <- voom(d,mm,plot=plot)

    if(plot){
        dev.off()
    }

    fit <- lmFit(y, mm)
    result <- vector(mode='list',length=length(contrast)) 
    names(result) <-contrast
    for( compare in contrast){
        s <-  gsub('Minus',paste0('-group'),compare)
        s <-  paste0('group',s)
        #contr <- makeContrasts(conrasts=s, levels = colnames(coef(fit)))
        contr <- do.call(makeContrasts,c(s,list(levels=colnames(coef(fit)))))
        tmp <- contrasts.fit(fit, contr)
        tmp <- eBayes(tmp)
        top.table <- topTable(tmp, sort.by = "P", n = Inf)
        result[[compare]] <- top.table
    }
    return(result)
}

parse_input = function(){
    require_pkgs <- c('data.table','edgeR','limma','argparse')
    scripts_wd <- dirname(thisFile())
    source(file.path(scripts_wd,'utils.R'))
    load_package(require_pkgs)

    parser <- ArgumentParser(description='Execute R bioconductors Limma Voom')
    inputs <- parser$add_argument_group('Input Option')
    inputs$add_argument('--counts_path',help='Counts Path')
    inputs$add_argument('--meta_path',required=T,help='Sample meta information Path')
    inputs$add_argument('--group_name',required=T,help='Column name in meta dataframe want to compare')
    inputs$add_argument('--contrast_path',required=T,help='Comparison Design Path')
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
