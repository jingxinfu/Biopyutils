#!/usr/bin/Rscript
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 12/02/2020
# Last Modified Date: 27/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>
DESCRIPTION = 'Immune cell abundance deconvolution'
main = function(){
  require_pkgs <- c('parallel','icbi-lab/immunedeconv','sva','data.table','argparse')
  scripts_wd <- dirname(thisFile())
  source(file.path(scripts_wd,'utils.R'))
  load_package(require_pkgs)

  set_cibersort_binary(file.path(scripts_wd,'immuneDeconv_Packed/CIBERSORT.R'))
  set_cibersort_mat(file.path(scripts_wd,"immuneDeconv_Packed/LM22.txt"))

  parser <- ArgumentParser(description=DESCRIPTION)
  parser$add_argument('--exprsn_path',help="[Required] Path to the expression profile (TPM-normalized, not log transformed)")
  parser$add_argument("--output", required=T,type="character", help="[Required] Output path")
  parser$add_argument("--cancer",required=T,type="character",help="[Required] TCGA cancer type")
  args <- parser$parse_args()

  exprsn <- fread(args$exprsn_path,data.table=F)
  exprsn <- exprsn[!is.na(exprsn[,1]),] # remove NA gene symbol
  rownames(exprsn) <- exprsn[,1]
  exprsn <- exprsn[,-1]


  if(args$cancer == 'AUTO'){
    cancerTypeVector = get(load(file.path(scripts_wd,"immuneDeconv_Packed/avgTPM.Rdata")))
    ol_genes = intersect(rownames(cancerTypeVector),rownames(exprsn))
    sp_correlation = na.omit(apply(cor(exprsn[ol_genes,],cancerTypeVector[ol_genes,],method='spearman',use='complete.obs'),
                           1,
                           function(x){
                               if(length(na.omit(x)) > 0){
                                   return(which.max(x))
                               }else{
                                   return(NA)
                               }
                           }))
    max_cnt = table(colnames(cancerTypeVector)[sp_correlation])
    cancer = names(max_cnt)[which.max(max_cnt)]
  }
  list_of_tools=c('quantiseq','timer','cibersort','cibersort_abs','mcp_counter','xcell','epic')
  result <- rbindlist(parallel::mclapply(list_of_tools,function(x) force(deconvAll(exprsn=exprsn,tool=x,cancer=cancer,scripts_wd=scripts_wd)),mc.cores=7))
  fwrite(result,file=args$output,sep=',')
}

TIMER <- function(e, cancer,scripts_wd) {
  load(file.path(scripts_wd,"immuneDeconv_Packed/timer_pureImmune.Rdata"))
  load(file.path(scripts_wd,"immuneDeconv_Packed/timer_geneMarker.Rdata"))
  RemoveBatchEffect <- function(cancer.exp, immune.exp, immune.cellType) {
    ## intersect the gene names of cancer.exp and immune.exp
    tmp.dd <- as.matrix(cancer.exp)
    tmp <- sapply(strsplit(rownames(cancer.exp), '\\|'),
                  function(x) x[[1]])
    rownames(tmp.dd) <- tmp
    tmp.dd <- as.matrix(tmp.dd[which(nchar(tmp)>1), ])
    tmp.ss <- intersect(rownames(tmp.dd), rownames(immune.exp))
    ## bind cancer and immune expression data into one dataframe
    N1 <- ncol(tmp.dd)
    tmp.dd <- cbind(tmp.dd[tmp.ss, ], immune.exp[tmp.ss, ])
    tmp.dd <- as.matrix(tmp.dd)
    mode(tmp.dd) <- 'numeric'
    ## remove batch effects
    N2 <- ncol(immune.exp)
    tmp.batch <- c(rep(1, N1), rep(2, N2))
    tmp.dd0 <- ComBat(tmp.dd, tmp.batch, c())
    ## separate cancer and immune expression data after batch effect removing
    dd.br <- tmp.dd0[, 1:N1]
    immune.exp.br <- tmp.dd0[, (N1+1):(N1+N2)]
    ## a immune category has multiple samples, use the median expression level for a gene
    tmp0 <- c()
    for(kk in unique(names(immune.cellType))){
      tmp.vv <- which(names(immune.cellType)==kk)
      tmp0 <- cbind(tmp0, apply(immune.exp.br[, tmp.vv], 1, median, na.rm=T))
    }
    immune.exp.agg.br <- tmp0
    colnames(immune.exp.agg.br) <- unique(names(immune.cellType))
    return(list(as.matrix(dd.br), immune.exp.br, immune.exp.agg.br))
  }
  GetFractions.Abbas <- function(XX, YY, w=NA){
    ## XX is immune expression data
    ## YY is cancer expression data
    ss.remove=c()
    ss.names=colnames(XX)
    while(T){
      if(length(ss.remove)==0)tmp.XX=XX else{
        if(is.null(ncol(tmp.XX)))return(setNames(rep(0, ncol(XX)), nm = colnames(XX)))
        tmp.XX=tmp.XX[, -ss.remove]
      }
      if(length(ss.remove)>0){
        ss.names=ss.names[-ss.remove]
        if(length(ss.names)==0)return(setNames(rep(0, ncol(XX)), nm = colnames(XX)))
      }
      if(is.na(w[1]))tmp=lsfit(tmp.XX, YY, intercept=F) else tmp=lsfit(tmp.XX, YY, w, intercept=F)
      if(is.null(ncol(tmp.XX)))tmp.beta=tmp$coefficients[1] else tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
      if(min(tmp.beta>0))break
      ss.remove=which.min(tmp.beta)
    }
    tmp.F=rep(0, ncol(XX))
    names(tmp.F)=colnames(XX)
    tmp.F[ss.names]=tmp.beta
    return(tmp.F)
  }
  immune.geneExpression <- immune$genes
  immune.cellTypes <- immune$celltypes
  cancer.expression <- e # replace file
  outlier.genes <- sort(unique(c(as.matrix(apply(cancer.expression, 2, function(x) {rownames(cancer.expression)[tail(order(x), 5)]})))))
  cancer.expression <- cancer.expression[!(rownames(cancer.expression) %in% outlier.genes), , drop=F]
  d.rmBatch <- RemoveBatchEffect(cancer.expression, immune.geneExpression, immune.cellTypes)
  cancer.expNorm <- d.rmBatch[[1]]
  immune.expNormMedian <- d.rmBatch[[3]]
  gene.selected.marker <- intersect(geneMarker[[cancer]], row.names(cancer.expNorm))
  X_immune = immune.expNormMedian[gene.selected.marker, c(-4)]
  Y_cancer = cancer.expNorm[gene.selected.marker, , drop=FALSE]
  as.data.frame(apply(Y_cancer, 2, function(y) {GetFractions.Abbas(X_immune, y)})) %>%  add_rownames(var='cell_type' )%>% as.data.frame()
}


deconvAll = function(exprsn,tool,cancer,scripts_wd){
  
    default_out <- as.data.frame(matrix(nrow=1,ncol=ncol(exprsn)+1),drop=F)
    colnames(default_out) <- c('cell_type',colnames(exprsn))
    out <- tryCatch({
          if (tool == "timer"){
            res = TIMER(exprsn, cancer = cancer,scripts_wd=scripts_wd)
          }else{
            res = deconvolute(exprsn,tool)
            res = as.data.frame(res)
          }
          tool_name = sub('_ABS$','-ABS',toupper(tool))[1]
          res$cell_type = paste(res$cell_type,tool_name,sep='_')
          res
        },error = function(e){
          default_out
      })
   
   return(out)
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
