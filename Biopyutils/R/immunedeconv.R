#!/usr/bin/Rscript
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 12/02/2020
# Last Modified Date: 08/06/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

DESCRIPTION = 'Immune cell abundance deconvolution'
main = function(){
  require_pkgs <- c('immunedeconv','sva','data.table','argparse','tibble','mMCPcounter','dplyr','plyr')
  for(p in require_pkgs){
      suppressMessages(require(p,character.only=T))
  }
  scripts_wd <- dirname(thisFile())
  parser <- ArgumentParser(description=DESCRIPTION)
  parser$add_argument('--exprsn_path',help="[Required] Path to the expression profile (TPM-normalized, not log transformed)")
  parser$add_argument("--output", required=T,type="character", help="[Required] Output path")
  parser$add_argument("--cancer",required=T,type="character",help="[Required] TCGA cancer type")
  parser$add_argument("--species",required=F,type="character",default='Human',help="[Option] Human or Mouse (Human by default)")
  args <- parser$parse_args()

  exprsn <- fread(args$exprsn_path,data.table=F)
  exprsn <- exprsn[!is.na(exprsn[,1]),] # remove NA gene symbol
  rownames(exprsn) <- exprsn[,1]
  exprsn <- exprsn[,-1]
  exprsn[is.na(exprsn)] <- 0
  if(args$species == 'Human') { 
    exprsn_hs <- exprsn
  } else if(args$species == 'Mouse') { 
    exprsn_hs <- mouseToHs(exprsn,scripts_wd=scripts_wd)
    exprsn_mm <- exprsn
  }

  signature <- ifelse(args$cancer == 'AUTO', signatureDetect(exprsn_hs,scripts_wd=scripts_wd), args$cancer)
  result <- list()
  for(tool in c('TIMER','CIBERSORT','CIBERSORT-ABS','QUANTISEQ','XCELL','EPIC')) {
    result[[tool]] <- deconvAll(exprsn=exprsn_hs,tool=tool,cancer=signature,scripts_wd=scripts_wd)
  }
  ### Mouse Data
  if(args$species == 'Human') {
    result[['MCPCOUNTER']] <- deconvAll(exprsn=exprsn_hs,tool='MCPCOUNTER',scripts_wd=scripts_wd)
  } else if(args$species == 'Mouse') {
    result[['MMCPCOUNTER']] <- deconvAll(exprsn=exprsn_mm,tool='MMCPCOUNTER',scripts_wd=scripts_wd)
  }

  est <- rbindlist(result) %>% as.data.frame()
  fwrite(est,file=args$output,sep=',')
}

deconvAll = function(exprsn, tool, cancer,scripts_wd){
  mcpcounter_genes=read.delim(file.path(scripts_wd,"immunedeconv_data/mcpcounter_genes.txt"), stringsAsFactors=FALSE, check.names=F)
  mmcpcounter_cellname_convertion <- c(
    'Fibroblasts' = 'Fibroblast',
    'B derived' = 'B cell',
    'Lymphatics' = 'Lymphoid',
    'Eosinophils' = 'Eosinophil',
    'Basophils' = 'Basophil',
    'Endothelial cells' = 'Endothelial cell',
    'Vessels' = 'Vessels',
    'Mast cells' = 'Mast cell',
    'Monocytes' = 'Monocyte',
    'T cells' = 'T cell',
    'Memory B cells' = 'B cell memory',
    'Monocytes / macrophages' = 'Macrophage/Monocyte',
    'CD8 T cells' = 'T cell CD8+',
    'NK cells' = 'NK cell',
    'Neutrophils' = 'Neutrophil',
    'Granulocytes' = 'Granulocyte'
  )
  TIMERdeconv <- function(e, cancer,scripts_wd) {
    load(file.path(scripts_wd,"immunedeconv_data/timer_pureImmune.Rdata"))
    load(file.path(scripts_wd,"immunedeconv_data/timer_geneMarker.Rdata"))
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
    transName <- c(B_cell = "B cell",
                   T_cell.CD4 = "T cell CD4+",
                   T_cell.CD8 = "T cell CD8+",
                   Neutrophil = "Neutrophil",
                   Macrophage = "Macrophage",
                   DC = "Myeloid dendritic cell")
    as.data.frame(apply(Y_cancer, 2, function(y) {GetFractions.Abbas(X_immune, y)})) %>%  
      dplyr::add_rownames(var='cell_type' ) %>% 
      dplyr::mutate(cell_type=plyr::revalue(cell_type, transName)) %>%
      as.data.frame()
  }
  set_cibersort_binary(file.path(scripts_wd,"immunedeconv_data/CIBERSORT.R"))
  set_cibersort_mat(file.path(scripts_wd,"immunedeconv_data/LM22.txt"))
  tryCatch({
    switch(
      tool,
      TIMER = TIMERdeconv(exprsn, cancer = cancer,scripts_wd=scripts_wd),
      CIBERSORT = deconvolute(exprsn, 'cibersort'),
      `CIBERSORT-ABS` = deconvolute(exprsn, 'cibersort_abs'),
      QUANTISEQ = deconvolute(exprsn, 'quantiseq'),
      MCPCOUNTER = deconvolute(exprsn, 'mcp_counter', probesets=NULL, genes=mcpcounter_genes),
      MMCPCOUNTER = mMCPcounter.estimate(exprsn) %>% tibble::rownames_to_column(var = 'cell_type') %>% mutate(cell_type = revalue(cell_type, mmcpcounter_cellname_convertion)),
      XCELL = deconvolute(exprsn, 'xcell') %>% as.data.frame(),
      EPIC = deconvolute(exprsn, 'epic') %>% as.data.frame()
    ) %>% mutate(cell_type = paste0(cell_type,'_',tool))
  },error = function(e){
    data.frame(matrix(NA, 1, ncol(exprsn)+1, dimnames = list(NULL, c('cell_type',colnames(exprsn)))), check.names = F)
  })
}

signatureDetect <- function(exprsn,scripts_wd) {
  cancerTypeVector = get(load(file.path(scripts_wd,"immunedeconv_data/avgTPM.Rdata")))
  ol_genes = intersect(rownames(cancerTypeVector),rownames(exprsn))
  sp_correlation = apply(cor(exprsn[ol_genes,],cancerTypeVector[ol_genes,],method='spearman',use='complete.obs'),1,which.max)
  max_cnt = table(colnames(cancerTypeVector)[sp_correlation])
  names(max_cnt)[which.max(max_cnt)]
}

mouseToHs <- function(exprsn,scripts_wd) {
# 	human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 	mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# 	load('submodule/humanSymbol.RData')
#   load('submodule/mouseSymbol.RData')
# 	mm2hs <- as.data.table(getLDS(
# 		attributes = c("mgi_symbol"), filters = "mgi_symbol",
# 		values = rownames(exprsn) , mart = mouse, attributesL = c("hgnc_symbol"),
# 		martL = human, uniqueRows=TRUE))
  load(file.path(scripts_wd,'immunedeconv_data/mm2hs.Rdata'))
  dt <- exprsn %>% 
    tibble::rownames_to_column(var='mgi_symbol') %>% 
    left_join(mm2hs) %>% 
    dplyr::select(-mgi_symbol) %>% 
    dplyr::select(hgnc_symbol, everything()) %>% 
    dplyr::filter(!is.na(hgnc_symbol)) %>% 
    ddply(.(hgnc_symbol), function(x) colMeans(x[,-1])) %>% 
    tibble::column_to_rownames(var = 'hgnc_symbol')
  dt
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
