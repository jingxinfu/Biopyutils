#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : GPL3
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 10/02/2020
# Last Modified Date: 11/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>
# -*- coding: utf-8 -*-
# Author            : Jingxin Fu <jingxin_fu@outlook.com>
# Date              : 07/02/2020
# Last Modified Date: 10/02/2020
# Last Modified By  : Jingxin Fu <jingxin_fu@outlook.com>

import os
from functools import wraps
from textwrap import dedent
import pandas as pd
import logging
from Biopyutils import getGeneRefPath, getSpeciesRefPath
__all__ = ['idConvert','speciesConvert']

_Comm_docs =dict(
        id_df = dedent("""\
            df: pd.DataFrame, pd.Series, or list
                A Pandas DataFrame/Series indexed by either Entrez Gene ID, Ensembl Gene ID, Ensembl Transcript ID, or Hugo Symbol.
                or a list of Entrez Gene ID, Ensembl Gene ID, Ensembl Transcript ID, or Hugo Symbol.\
            """),
        map_id = dedent("""\
            map_id: str
                Which gene identifier that you want to convert the index of/entries of input `df` to.
                Option: Entrez, Symbol, ENSG
                Here the Entrez stands for Entrez Gene ID; the Symbol stands for Hugo Symbol;
                the ENSG stands for Ensembl Gene ID;\
            """),
        show_num = dedent("""\
            show_num: int, optional
                Number of idenfiers without corresponding `map_id` that will be print on logging stdout (if logger is not None)
                Default is 10.\
                """),
        logger = dedent("""\
            logger: logging.Logger or None, optional
                Logging.Logger to record runtime information. Doesn't record runtime information if it is None.
                Default is None.\
                """)
            )

def Keepdtype(func):
    """Keep the original data type"""
    @wraps(func)
    def wrapper(*args,**kwargs):
        list_flag = False
        series_flag = False
        arg_tuple_flag = False
        if not 'df' in kwargs.keys():
            df = args[0]
            arg_tuple_flag = True
        else:
            df = kwargs['df']
        if isinstance(df,list):
            list_flag = True
            df = pd.Series([1]*len(df),index=df).to_frame() # random assign value
        elif isinstance(df,pd.Series):
            series_flag = True
            df = df.to_frame()
        elif isinstance(df,pd.DataFrame):
            pass
        else:
            raise ValueError('Only support df to be pandas.Series, pandas.DataFrame, or a list')

        # Assign new value
        if arg_tuple_flag:
            args = list(args)
            args[0] = df
            args = tuple(args)
        else:
            kwargs['df'] = df
        result = func(*args,**kwargs)
        if series_flag:
            return result.squeeze()
        elif list_flag:
            return result.index.tolist()
        else:
            return result

    return wrapper

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

    try:
        df.index = df.index.map(int)
        return 'Entrez'
    except:
        id_source = pd.Series({
            'ENSG':df.index.str.contains(r'^ENS\w*G\d',regex=True).sum(),
            'ENST':df.index.str.contains(r'^ENS\w*T\d',regex=True).sum(),
            'Symbol': (~df.index.str.contains('^ENS*',regex=True)).sum(),
            })
        return id_source.idxmax()

def infoMissing(old,new,miss_ins,show_num=10,logger=None):
    """infoMissing checks number of missing entries and throws an error when total entries are missed

    Parameters
    ----------
    old_df : pd.DataFrame or pd.Series
        old is the amount of old entries
    new_df : int or float
        new is the amount of new entries
    miss_ins: list
        list of missing entries in old_df
    show_num: int
        Number of missing entries showing on logger
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
        if show_num > 0:
            show_num = show_num if len(miss_ins) > show_num else len(miss_ins)
            show_text = ', '.join(miss_ins[:show_num])
            logger.warning('%d out of %d (%.2f%%) source id do not have destination id to map. They are %s.' % (old,new,miss_ratio*100,show_text))
        else:
            logger.warning('%d out of %d (%.2f%%) source id do not have destination id to map' % (old,new,miss_ratio*100))

@Keepdtype
def speciesConvert(df,from_species,to_species,logger=None,show_num=10):
    """speciesConvert converts gene id from one species to another species """
    if from_species == to_species:
        return df

    if len(set([from_species, to_species] + ['hg','mm'])) != 2:
        raise ValueError('Only accept conversion between "hg" and "mm"')
    map_id = inferIDsource(df)
    ref = pd.read_pickle(getSpeciesRefPath(map_id=map_id))
    ref = ref.drop_duplicates(subset=[from_species]).set_index(from_species)
    miss_ins = df.index.difference(ref.index).tolist()
    result = df.merge(ref,left_index=True,right_on=from_species)
    infoMissing(old=df.shape[0],new=result.shape[0],miss_ins=miss_ins,show_num=show_num,logger=logger)

    return result.groupby(to_species).mean()

speciesConvert.__doc__ = dedent("""\

        Convert a list of one gene identifiers from one species to the other species.
        If the input `df` is pandas.DataFrame or pandas.Series, all numberic columns will be grouped by the correponding gene identifier of `to_species`  and taken the average (sum, if the original gene identifier is Ensemble transcript ID.)

        Parameters
        ----------
        {id_df}
        from_species : str
            from_species is the name of species that `df`'s gene identifiers belong to.
        to_species : str
            to_species is the name of species you want to convert your gene identifiers to
        {show_num}
        {logger}

        Returns
        ----------
        list, pd.Series, or pd.DataFrame (depend on input data type)
            if it's pd.Series or pd.DataFrame, then the output should be pd.Series or pd.DataFrame that indexed by the corresponding gene identifiers of `to_species`
            if it's a list, the the output will be a list of the corresponding gene identifiers of `to_species`.\


        Examples
        ----------
        Convert Mouse Entrez ID to Human Entrez ID:
            >>> import pandas as pd
            >>> from Biopyutils import Comm
            >>> mouse_df = pd.DataFrame([1]*3,index=[26695,381308,670895])
            >>> hg_df = Comm.speciesConvert(df=mouse_df,from_species='mm',to_species='hg')
            >>> hg_df.head()
                0
            hg
            4332    1
            386672  1
            >>> mouse_list = mouse_df.index.to_list()
            >>> hg_list = Comm.speciesConvert(df=mouse_list,from_species='mm',to_species='hg')
            >>> print(hg_list)
            [4332, 386672]
        """).format(**_Comm_docs)

@Keepdtype
def idConvert(df,species,map_id,logger=None,show_num=10):
    """idConvert converts gene id from a kind of gene identifiers to <map_id> gene identifiers"""
    if not species in ['hg','mm']:
        raise ValueError('Only do conversion for "hg" and "mm" species.')
    source_id = inferIDsource(df)
    if source_id == map_id:
        return df

    ref = pd.read_pickle(getGeneRefPath(species=species,map_id=map_id))[source_id]
    miss_ins = df.index.difference(ref.index).tolist()
    result = df.merge(ref,left_index=True,right_on=source_id)
    infoMissing(old=df.shape[0],new=result.shape[0],miss_ins=miss_ins,show_num=show_num,logger=logger)
    if source_id == 'ENST':
        result = result.groupby(map_id).sum()
    else:
        result = result.groupby(map_id).mean()

    return result

idConvert.__doc__ = dedent("""\
        Convert a list of one gene identifiers to the other one.
        If the input `df` is pandas.DataFrame or pandas.Series, all numberic columns will be grouped by the `map_id` and taken the average (sum, if the original gene identifier is Ensemble transcript ID.)

        Parameters
        ----------
        {id_df}
        species : str
            Which species the input `df` belongs to.
        {map_id}
        {show_num}
        {logger}

        Returns
        ----------
        list, pd.Series, or pd.DataFrame (depend on input data type)
            if it's pd.Series or pd.DataFrame, then the output should be pd.Series or pd.DataFrame that indexed by `map_id`
            if it's a list, the the output will be a list of `map_id` that original list maps to.\


        Examples
        ----------
        Convert Human Hugo Symbol to Human Entrez ID:
            >>> import pandas as pd
            >>> from Biopyutils import Comm
            >>> symbol_df = pd.DataFrame([1]*3,index=['PDCD1','CD274','COP1'])
            >>> entrez_df = Comm.idConvert(df=symbol_df,map_id='Entrez',species='hg')
            >>> entrez_df.head()
                0
            Entrez
            5133.0   1
            29126.0  1
            64326.0  1
        Convert Mouse Ensembl transcript ID to Mouse Entrez ID:
            >>> enst_df = pd.DataFrame([1]*3,index=['ENSMUST00000026432','ENSMUST00000163106','NSMUST0000020529'])
            >>> entrez_df = Comm.idConvert(df=enst_df,map_id='Entrez',species='mm')
            >>> entrez_df.head()
                0
            Entrez
            12565.0      1
            100170401.0  1
            >>> enst_list = enst_df.index.tolist()
            >>> entrez_list = Comm.idConvert(df=enst_list,map_id='Entrez',species='mm')
            >>> print(entrez_list)
            [12565.0, 100170401.0]
        """).format(**_Comm_docs)

