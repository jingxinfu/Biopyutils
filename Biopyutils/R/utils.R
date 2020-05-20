#!/usr/bin/Rscript
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 14/02/2020
# Last Modified Date: 13/04/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

load_package = function(pkgs){
    if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager",repos = "https://cran.us.r-project.org")
        BiocManager::install(version = "3.10")
    }
    if (!requireNamespace("remotes", quietly = TRUE)){
        install.packages("remotes",repos = "https://cran.us.r-project.org")
    }
    #if(!is.element('BiocManager', installed.packages()[,1])){
    #        install.packages('BiocManager',repos = "http://cran.us.r-project.org")
    #}
    for(el in pkgs){
        if (!is.element(el, installed.packages()[,1])){
            if(grepl('/',el)){
                remotes::install_github(el,upgrade ='never')
                el <- strsplit(el,split="/")[[1]][2]
            }else{
                BiocManager::install(el)
            }
        }
        invisible(require(el, character.only=TRUE))
    }
}
