#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : GPL3
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 07/02/2020
# Last Modified Date: 07/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

import os
import logging
from Biopyutils import getGeneRefPath, getSpeciesRefPath

def Keepdtype(func):
    """Keep the original data type"""
    def wrapper(*args,**kwargs):
        df = kwargs['df']
        list_flag = True
        series_flag = True
        if isinstance(df,list):
            list_flag = True
            df = pd.Series([1]*len(df),index=df) # random assign value
        elif isinstance(df,pd.Series):
            series_flag = True
            df = df.to_frame()
        elif isinstance(df,pd.DataFrame):
            pass
        else:
            raise ValueError('Only support df to be pandas.Series, pandas.DataFrame, or a list')
        kwargs['df'] = df
        result = func(*args,**kwargs)
        if series_flag:
            return result.squeeze()
        elif list_flag:
            return result.index.tolist()
        else:
            return result

    return wrapper()

def inferIDsource(df):
    """ inferIDsource infers the name of gene identifiers for input data index
    Parameters
    ----------
    df : pd.DataFrame or pd.Series
        df with index by gene identifiers

    Returns:
    ----------
    str
        The name of gene identifiers: ENSG, ENST, Symbol, or Entrez
    """

    cnt_nonint = sum(map(lambda v: not isinstance(v, (int, float)), df.index)),
    if cnt_nonint == 0:
        return 'Entrez'
    else:
        id_source = pd.Series({
            'ENSG':df.index.str.contains('^ENS*T*',regex=True).sum(),
            'ENST':df.index.str.contains('^ENS*G*',regex=True).sum(),
            'Symbol': (~df.index.str.contains('^ENS*',regex=True).sum(),
            }).idxmax()
        return id_source

def infoMissing(old,new,logger=None):
    """infoMissing checks number of missing entries and throws an error when total entries are missed

    Parameters
    ----------
    old : int or float
        old is the amount of old entries
    new : int or float
        new is the amount of new entries
    logger : None or logging.Logger
        logger is the logiing obj to sending out runing info

    Returns:
    ----------
    None

    Raises:
    ----------
    ValueError: When the amount of new entries is 0
    """
    miss_ratio = 1 - (new/float(old))
    if new == 0:
        raise ValueError('0 out of %d source id has destination id' % old)
    if isinstance(logger,logging.Logger) and miss_ratio > 0.1:
        logger.warn('%d out of %d (%.2f%%) source id do not have destination id to map' % (old,new,miss_ratio*100))

@Keepdtype
def specieConvert(df,from_species,to_species,logger=None):
    """specieConvert converts gene id from one species to another species

    Parameters
    ----------
    df : list, pd.Series, or pd.DataFrame
        df is a list of gene identifiers
        or pd.Series or pd.DataFrame indexed by a list of gene identifiers
    from_species : str
        from_species is the name of species of gene identifiers in df
    to_species : str
        to_species is the name of species you want to convert your gene identifiers to
    logger : None or logging.Logger
        logger is a logging obj to record run info

    Returns:
    ----------
    list, pd.Series, or pd.DataFrame (depend on input data type)
        The <to_species> gene identifiers

    """
    if not all([from_species, to_species] in ['hg','mm'])
        raise ValueError('Only accept conversion between "hg" and "mm"')
    if from_species == to_species:
        return df
    map_id = inferIDsource(df)
    ref = getSpeciesRefPath(map_id=map_id)
    result = df.merge(ref,left_index=True,right_on=from_species)
    infoMissing(old=df.shape[0],new=result.shape[0],logger=logger)

    return result.groupby(to_species).mean()


@Keepdtype
def idConvert(df,species,map_id,logger=None):
    """idConvert converts gene id from a kind of gene identifiers to <map_id> gene identifiers

    Parameters
    ----------
    df : list, pd.Series, or pd.DataFrame
        df is a list of gene identifiers
        or pd.Series or pd.DataFrame indexed by a list of gene identifiers
    species : str
        species is the name of species of gene identifiers in df
    map_id : str
        map_id is the name of gene identifiers you wants to convert to
    logger : None or logging.Logger
        logger is a logging obj to record run info

    Returns:
    ----------
    list, pd.Series, or pd.DataFrame (depend on input data type)
        The <map_id> gene identifiers

    """
    if not species in ['hg','mm']:
        raise ValueError('Only do conversion for "hg" and "mm" species.')
    source_id = inferIDsource(df)
    if source_id == map_id:
        return df

    ref = getGeneRefPath(species=species,map_id=map_id)
    result = df.merge(ref,left_index=True,right_on=source_id)
    infoMissing(old=df.shape[0],new=result.shape[0],logger=logger)
    if source_id == 'ENST':
        result = result.groupby(map_id).sum()
    else:
        result = result.groupby(to_species).mean()

    return result



