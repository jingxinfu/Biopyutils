#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : GPL3
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 11/02/2020
# Last Modified Date: 12/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

import tempfile
import subprocess
import pandas as pd
from Biopyutils.Comm import Rscript


#from Biopyutils import R_Dir
def batchRemoval():
    pass
def normlizeTPM():
    pass

def gsva(exprsn,geneset,method='ssgsea',kcdf='Guassian'):#exprsn,method,kcdf):
    with tempfile.NamedTemporaryFile() as e_f:
        exprsn.to_csv(e_f.name)
        with tempfile.NamedTemporaryFile() as bio_f:
            geneset.to_csv(bio_f.name)
            cmd = "GSVA.R --kcdf %s --method %s --exprsn %s --gset %s" % ('R/',kcdf,method,e_f.name,bio_f.name)
            result = pd.read_csv(Rscript(cmds))
    return result



def deseq2():
    pass


class RNASeq:
    def __init__(self,tpm,meta,is_count=False):
        self.tpm = tpm
        if is_count:
            self._toTPM()

        self.normTpm = self.tpm
        self.meta = meta

    def _toTPM(self.tpm):
        count_div_length = self.tpm / Comm.idConvert(,map_id='length')
        self.tpm = count_div_length.T * 1e6 )/count_div_length.sum(axis=0)

    def removeBatch(self,batch):
        self.normaTpm = batchRemoval(self.tpm,batch)

    def diff(self,design):
        deseq2(self.tpm,design)

    def gsva(self,method='gsva'):
        gsva(self.exprsn,method=method,kcdf='Guassian')



