#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 11/02/2020
# Last Modified Date: 08/06/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

import tempfile
import scipy
import os
import pandas as pd
import numpy as np
from sklearn import preprocessing,decomposition,linear_model,ensemble
from Biopyutils.Comm import Rscript,lmFit,quantileNorm,fdrAdjust

__all__ = ['bulkRNASeq']
class _BaseExprsn:
    def __init__(self,tpm,counts=None,cancer=None,meta=None,is_normalized=False):
        self.counts = counts
        self.meta = meta
        if is_normalized:
            self.tpm = None
            self.logTpm = None
            self.normTpm = tpm
        else:
            self.tpm = tpm
            self.logTpm = np.log2(1+tpm)
            self.normTpm = self.exprsnNorm()

        self.cancer = cancer
        self.DEG = None
        self.signature = None
        self.infiltras = None

    def exprsnNorm(self):
        return rnaSeqNorm(self.tpm)

    def _toTPM(self):
        pass
        #count_div_length = self.tpm / Comm.idConvert(,map_id='length')
        #self.tpm = count_div_length.T * 1e6 )/count_div_length.sum(axis=0)

    def removeBatch(self,fig_out=None):
        '''This returns an expression matrix, with the same dimensions as your original
        dataset. This new expression matrix has been adjusted for batch'''
        if self.meta is None:
            raise ValueError("Samples meta information is required to do batch effect removal.")
        if not 'Batch' in self.meta.columns:
            raise KeyError("Cannot find 'Batch' column in the meta data frame.")
        self.tpm = batchRemoval(self.tpm,self.meta,fig_out=fig_out)
        fig_out_log = fig_out+'_log' if fig_out else None
        self.logTpm = batchRemoval(self.logTpm,self.meta,fig_out=fig_out_log)



    def getDEG(self,group_name,contrast,force_limma=False,fig_out=None):
        '''Do differential gene expression'''
        if not group_name in self.meta.columns:
            raise KeyError("Cannot find a column named '%s' in the meta dataframe." % group_name)
        if self.tpm.shape[1] >=10 and (not force_limma):
            self.DEG = regDEG(logTpm=self.logTpm,meta=self.meta,group_name=group_name,contrast=contrast)
        elif self.counts is None:
            raise ValueError('Please provide count-base expression profiles.')
        else:
            self.DEG = limmaVoom(self.counts,self.meta,group_name,contrast,fig_out=fig_out)

    def getSignature(self,geneset):
        self.signature = gsva(self.tpm,geneset=geneset,method='ssgsea',kcdf='Gaussian')

    def getInfiltras(self):
        if self.tpm is None:
            return None
        else:
            cancer = 'AUTO' if self.cancer is None else self.cancer
            self.infiltras = immuneDeconv(self.tpm,cancer=cancer)

class bulkRNASeq(_BaseExprsn):
    def __init__(self,tpm,counts=None,cancer=None,meta=None,is_normalized=False):
        super().__init__(tpm,counts,cancer,meta,is_normalized)


def rnaSeqNorm(df):
    """ Normalization steps for RNAseq data
    1. Log2(1+x) transform
    2. Quantile normalization across samples
    3. Centralize distribution of gene expression across samples by subtracting the its average across samples

    Parameters
    ----------
    df : pd.DataFrame
       Raw TPM, indexed by gene name and columned by sample id

    Returns:
    ----------
    pd.DataFrame
        Normalized Expression profile
    """
    df = quantileNorm(np.log2(1+df))
    df = df.subtract(df.mean(axis=1),axis=0)
    return df

def regDEG(logTpm,meta,group_name,contrast):
    """ Using simple Wilcoxon-Test to find differetial Genes

    Parameters
    ----------
    logTpm: pd.DataFrame
        logTPM, indexed by gene name and columned by sample id
    meta : pd.DataFrame
        Meta infomration
    group_name : str
        Column names that has group information
    contrast: list
        A list of comparison

    Returns:
    ----------
    Differential Expression Coefficient and statistic significance
    """
    # Remove low expressed genes and low variant genes
    logTpm = logTpm.loc[logTpm.var(axis=1)> 0.1,:]
    logTpm = logTpm.loc[logTpm.sum(axis=1)>1,:]
    scaler = preprocessing.StandardScaler()

    result = {}
    for statement in contrast:
        # Factorizing contrast to 1,0
        factors = statement.split('Minus')
        factors_map = {factors[0]:1,factors[1]:-1}
        Y  = meta.loc[logTpm.columns,group_name].map(factors_map).dropna(axis=0)
        # replace nan with zero and inf with finite numbers
        rawX = pd.DataFrame(np.nan_to_num(logTpm.loc[:,Y.index]),index=logTpm.index,columns=Y.index)
        logFc =  rawX.apply(lambda v:np.mean(v[Y==1]) - np.mean(v[Y==-1]),axis=1).values
        pvalue = rawX.apply(lambda v: scipy.stats.mannwhitneyu(x=v[Y==1],y=v[Y==-1],alternative='two-sided').pvalue,axis=1).ravel()
        zscore = scaler.fit_transform(logFc.reshape(-1,1)).ravel()
        fdr = fdrAdjust(pvalue)
        result[statement] = pd.DataFrame({
            'logFc':logFc,'Z-Score':zscore,
            'Pvalue':pvalue,'FDR':fdr
        },index=logTpm.index)

    return result

# R Function API
def batchRemoval(exprsn,meta,fig_out=None):
    if meta.columns.intersection(['Batch']).size != 1:
        raise KeyError('Please provide Batch information.')
    with tempfile.TemporaryDirectory(prefix='batchRemovel') as dirpath:
        params= dict(
            exprsn_path = os.path.join(dirpath,'exprsn'),
            meta_path = os.path.join(dirpath,'meta'),
            output = os.path.join(dirpath,'out')
        )
        if not fig_out is None:
            params['fig_out'] = fig_out
        exprsn.to_csv(params['exprsn_path'])
        meta.to_csv(params['meta_path'])
        Rscript(cmd='BatchRemoval.R',params=params)
        result = pd.read_csv(params['output'],index_col=0)

    return result

def gsva(exprsn,geneset,method='ssgsea',kcdf='Gaussian'):
    with tempfile.TemporaryDirectory(prefix='gsva') as dirpath:
        params= dict(
            method=method,
            kcdf=kcdf,
            exprsn_path = os.path.join(dirpath,'exprsn'),
            geneset_path = os.path.join(dirpath,'geneset'),
            output= os.path.join(dirpath,'out')
        )
        exprsn.to_csv(params['exprsn_path'])
        geneset.to_csv(params['geneset_path'],sep='\t',index=False,header=False)

        Rscript(cmd='GSVA.R',params=params)
        result = pd.read_csv(params['output'],index_col=0)

    return result



def limmaVoom(counts,meta,group_name,contrast,fig_out=None):
    with tempfile.TemporaryDirectory(prefix='limmaVoom') as dirpath:
        params= dict(
            counts_path = os.path.join(dirpath,'counts'),
            meta_path = os.path.join(dirpath,'meta'),
            group_name = group_name,
            contrast_path = os.path.join(dirpath,'contrast'),
            output = os.path.join(dirpath,'out')
        )
        if not fig_out is None:
            params['fig_out'] = fig_out
        counts.to_csv(params['counts_path'])
        meta.to_csv(params['meta_path'])
        with open(params['contrast_path'],'w') as f:
            f.write('\n'.join(contrast))

        Rscript(cmd='LimmaVoom.R',params=params)
        result = dict()
        for ct in contrast:
            result[ct] = pd.read_csv(params['output']+ct,index_col=0)

    return result

def immuneDeconv(exprsn,cancer):
    with tempfile.TemporaryDirectory(prefix='immuneDeconv') as dirpath:
        params= dict(
            exprsn_path = os.path.join(dirpath,'exprsn'),
            cancer = cancer,
            output = os.path.join(dirpath,'out')
        )
        exprsn.to_csv(params['exprsn_path'])

        Rscript(cmd='immunedeconv.R',params=params)
        result = pd.read_csv(params['output'],index_col=0)
    return result

