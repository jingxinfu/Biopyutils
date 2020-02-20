#!/usr/bin/Rscript
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 14/02/2020
# Last Modified Date: 20/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

load_package = function(pkgs){
    if(!is.element('BiocManager', installed.packages()[,1])){
            install.packages('BiocManager',repos = "http://cran.us.r-project.org")
    }
    for(el in pkgs){
            if (!is.element(el, installed.packages()[,1])){
                BiocManager::install(el)
            }
            invisible(require(el, character.only=TRUE))
    }
}


