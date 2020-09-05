#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# License           : GPL3
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 10/02/2020
# Last Modified Date: 10/02/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>
# -*- coding: utf-8 -*-
# Author            : Jingxin Fu <jingxin_fu@outlook.com>
# Date              : 09/02/2020
# Last Modified Date: 10/02/2020
# Last Modified By  : Jingxin Fu <jingxin_fu@outlook.com>
__doc__="""
# Stable Table Fetching
Original Data Source:http://useast.ensembl.org/biomart
## Ensl_hg38.gz and Ensl_mm10.gz:
1. Select Ensembl Gene 101
2. Choose Human genes (GRCh38.p13) for hg38 / Mouse genes (GRCm38.p6) for mm10
    GENE Panel:
        - Select Gene stable ID
        - Select Transcript stable ID
        - Select Gene name
    Extenal Reference Panel:
        - Select NCBI Gene ID
        - Select NCBI gene description
Download the result (compressed file CSV , Unique deselected)

## Ensl_mm10_hg38_match.gz:
1. Select Ensembl Gene 101
2. Choose Human genes (GRCh38.p13)
    GENE Panel:
        - Select Gene stable ID
    Extenal Reference Panel:
        - Select NCBI Gene ID
3. Add additional dataset:  Choose Mouse genes (GRCm38.p6)
    GENE Panel:
        - Select Gene stable ID
    Extenal Reference Panel:
        - Select NCBI Gene ID
Download the result (compressed file CSV , Unique deselected)

# Alias Symbols Fetching
Download Alias and Previous gene symbols

## Hugo_hg38.gz and NCBI_hg38.gz
    Source1: wget -O Hugo_hg38 "https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=gd_name_aliases&status=Approved&status=Entry%20Withdrawn&order_by=gd_app_sym_sort&format=text&submit=submit"

    Source2: https://www-ncbi-nlm-nih-gov.ezp-prod1.hul.harvard.edu/gene
        Search (Homo sapiens genome) AND "Homo sapiens"[porgn:__txid9606]
        - send to >> file >> tabular(text)
        rename the file to be NCBI_hg38

    gzip Hugo_hg38
    gzip NCBI_hg38

## NCBI_mm10.gz
    Source: https://www-ncbi-nlm-nih-gov.ezp-prod1.hul.harvard.edu/gene
    Search (Mus musculus genome) AND "Mus musculus"[porgn:__txid10090]
    - send to >> file >> tabular(text)
    rename the file to be NCBI_mm10
    gzip NCBI_mm10
"""

import os
import argparse
import pandas as pd

def getAliasRef(df,alias_col,ref_col,separator=', '):
    """ Explore alias map to unique reference ID (Entrez ID)

    Parameters
    ----------
    df : pd.DataFrame
        df is a dataframe contains the <alias_col> which gene alias are separated by <separator>
    alias_col : str
        alias_col is the name of the column having alias separated by <separator>
    ref_col: str
        ref_col is the name of the column having stable ID to map
    separator : str
        separator is a str that separates a list of alias
    Returns:
    ----------
    pd.DataFrame
       Unique map alias to stable gene ID
    """
    sub_df = df[[alias_col,ref_col]].dropna().copy()
    sub_df[alias_col] = sub_df[alias_col].map(lambda x:x.split(separator))
    sub_df = sub_df.explode(alias_col).drop_duplicates(subset=[alias_col])
    if ref_col == 'Entrez': # convert to int64
        sub_df[ref_col] = sub_df[ref_col].astype('int64')

    return sub_df

def genTestCase(ref,out):
    """genTestCase
    Generate file for test
    Parameters
    ----------
    ref : dict
        df is a dict
        {
        'hg':pd.DataFrame, for human gene id convert
        'mm':pd.DataFrame, for mouse gene id convert
        'match':pd.DataFrame, for mouse to human gene id convert
        }
    out : str
        Output folder
    Returns:
    ----------
    """
    identifiers =  ['ENST','ENSG','Symbol','Name','Entrez']
    id_for_species = ['ENSG','Entrez']
    test_num = 10
    # Gene id convert
    for sp in ['hg','mm']:
        for map_id in identifiers:
            if map_id == 'ENST':
                continue # Don't create <id> to ENST map, since ENST ids are normally mapped to multiple Gene associated ID
            tmp_ref = {}
            for source_id in identifiers:
                if source_id == 'Name':
                    continue
                if map_id == source_id:
                    continue
                tmp = ref[sp][[ source_id,map_id ]].drop_duplicates().dropna().copy() # Remove entire duplicated rows and rows with NA
                if 'Entrez' in [[source_id,map_id]] :
                    tmp['Entrez'] = tmp['Entrez'].astype('int64')
                # Extend map
                if source_id+'Alias' in ref[sp].columns:
                    extend_ref = getAliasRef(ref[sp],alias_col=source_id+'Alias',ref_col=map_id)
                    extend_ref.rename(columns={source_id+'Alias':source_id},inplace=True)
                    tmp = pd.concat([tmp,extend_ref],axis=0,ignore_index=True) # put tmp variable at first to keep official id when there is duplication issue

                ### Testing Files
                prefix = os.path.join(out,sp+'_'+source_id+'_to_'+map_id+'_')
                source_duplicates = tmp.loc[tmp[source_id].duplicated(),source_id].unique()
                tmp = tmp.drop_duplicates(subset=[source_id])
                map_uniq = ~tmp[map_id].duplicated()

                # Generate OneToOne map
                test = tmp.loc[(map_uniq),[source_id,map_id]].iloc[1:test_num,]
                test['Count']  = 1
                test.groupby(source_id)['Count'].mean().to_frame().to_csv(prefix+'OneToOne.input')
                test.groupby(map_id)['Count'].mean().to_frame().to_csv(prefix+'OneToOne.output')

                # Generate OneToMultiple Map
                test = tmp.loc[(tmp[source_id].isin(source_duplicates)) & (map_uniq),[source_id,map_id]].iloc[1:test_num,] # Choose First Index
                test['Count']  = 1
                test.groupby(source_id)['Count'].mean().to_frame().to_csv(prefix+'OneToMultiple.input')
                test.groupby(map_id)['Count'].mean().to_frame().to_csv(prefix+'OneToMultiple.output')

                # Generate MutipleTo One Map
                test = tmp.loc[ (~map_uniq),[source_id,map_id]].iloc[1:test_num,]
                test['Count']=1
                if source_id == 'ENST':
                    test.groupby(source_id)['Count'].mean().to_frame().to_csv(prefix+'MultipleToOne.input')
                    test.groupby(map_id)['Count'].sum().to_frame().to_csv(prefix+'MultipleToOne.output')
                else:
                    test.groupby(source_id)['Count'].mean().to_frame().to_csv(prefix+'MultipleToOne.input')
                    test.groupby(map_id)['Count'].mean().to_frame().to_csv(prefix+'MultipleToOne.output')


    # Species convert
    for gid in id_for_species:
        tmp  = ref['match'][['hg'+gid,'mm'+gid]].drop_duplicates().dropna()
        tmp.columns = tmp.columns.map(lambda x:x.replace(gid,''))
        if gid == 'Entrez':
            tmp = tmp.astype('int64')
        for source_id in ['hg','mm']:
            for map_id in ['hg','mm']:
                if source_id == map_id:
                    continue
                prefix = os.path.join(out,gid+'_'+source_id+'_to_'+map_id+'_')
                source_uniq = ~tmp[source_id].duplicated()
                map_uniq = ~tmp[map_id].duplicated()
                # Generate OneToOne map
                test = tmp.loc[ (source_uniq) & (map_uniq),[source_id,map_id]].iloc[1:test_num,]
                test['Count']  = 1
                test.groupby(source_id)['Count'].mean().to_frame().to_csv(prefix+'OneToOne.input')
                test.groupby(map_id)['Count'].mean().to_frame().to_csv(prefix+'OneToOne.output')
                # Generate OneToMultiple Map
                test = tmp.loc[(~source_uniq) & (map_uniq),[source_id,map_id]].iloc[1:test_num,]
                test['Count']  = 1
                test = test.drop_duplicates(subset=[source_id]) # Choose First Index
                test.groupby(source_id)['Count'].mean().to_frame().to_csv(prefix+'OneToMultiple.input')
                test.groupby(map_id)['Count'].mean().to_frame().to_csv(prefix+'OneToMultiple.output')

                # Generate MutipleTo One Map
                test = tmp.loc[ (source_uniq) & (~map_uniq),[source_id,map_id]].iloc[1:test_num,]
                test['Count']=1
                test.groupby(source_id)['Count'].mean().to_frame().to_csv(prefix+'MultipleToOne.input')
                test.groupby(map_id)['Count'].mean().to_frame().to_csv(prefix+'MultipleToOne.output')





def constructRefDb(ref,out):
    """constructRefDb
    Generate reference database for id mapping
    Parameters
    ----------
    ref : dict
        df is a dict
        {
        'hg':pd.DataFrame, for human gene id convert
        'mm':pd.DataFrame, for mouse gene id convert
        'match':pd.DataFrame, for mouse to human gene id convert
        }
    out : str
        Output folder
    Returns:
    ----------
    """
    identifiers =  ['ENST','ENSG','Symbol','Name','Entrez']
    id_for_species = ['ENSG','Entrez']
    # Gene id convert

    for sp in ['hg','mm']:
        for map_id in identifiers:
            if map_id == 'ENST':
                continue # Don't create <id> to ENST map, since ENST ids are normally mapped to multiple Gene associated ID
            tmp_ref = {}
            for source_id in identifiers:
                if source_id == 'Name':
                    continue
                if map_id == source_id:
                    continue
                tmp = ref[sp][[ source_id,map_id ]].drop_duplicates().dropna().copy() # Remove entire duplicated rows and rows with NA
                if 'Entrez' in [[source_id,map_id]] :
                    tmp['Entrez'] = tmp['Entrez'].astype('int64')
                # Extend map
                if source_id+'Alias' in ref[sp].columns:
                    extend_ref = getAliasRef(ref[sp],alias_col=source_id+'Alias',ref_col=map_id)
                    extend_ref.rename(columns={source_id+'Alias':source_id},inplace=True)
                    tmp = pd.concat([tmp,extend_ref],axis=0,ignore_index=True) # put tmp variable at first to keep official id when there is duplication issue
                tmp.drop_duplicates(subset=[source_id],inplace=True)
                tmp_ref[source_id] = tmp.set_index(source_id)

            pd.Series(tmp_ref).to_pickle(os.path.join(out,'geneId_%s_%s.pickle.gz' % (sp,map_id)))

    # Species convert
    for gid in id_for_species:
        tmp  = ref['match'][['hg'+gid,'mm'+gid]].drop_duplicates().dropna()
        tmp.columns = tmp.columns.map(lambda x:x.replace(gid,''))
        if gid == 'Entrez':
            tmp = tmp.astype('int64')

        tmp.to_pickle(os.path.join(out,'species_hg-mm_%s.pickle.gz' % (gid)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dir",type=str,help="Path to folder contains requied files")
    parser.add_argument("-o", "--output",type=str,required=True,help="Path to output prefix")
    parser.add_argument("-t", "--testOutput",type=str,required=True,help="Path to test case prefix")
    args = parser.parse_args()

    ensl_rename = {
        'Gene stable ID':'ENSG',
        'Transcript stable ID':'ENST',
        # 'NCBI gene ID':'Entrez',
        'NCBI gene (formerly Entrezgene) ID':'Entrez',
        'Gene name':'Symbol',
        # 'NCBI gene description':'Name'
        'NCBI gene (formerly Entrezgene) description':'Name'
    }
    hugo_rename ={
        'Approved symbol':'Symbol',
        'Approved name':'Name',
        'Previous symbols':'Previous',
        'NCBI Gene ID':'Entrez',
        'Ensembl gene ID':'ENSG',
        'Alias names':'SymbolAlias',
    }
    ncbi_rename = {
        'GeneID':'Entrez',
        'Symbol:':'Symbol',
        'Aliases':'SymbolAlias'
    }

    # load data
    hg_ref={
        'ensl':pd.read_csv(os.path.join(args.dir,'Ensl_hg38.gz')).rename(columns=ensl_rename).dropna(subset=['Entrez']),
        'hugo':pd.read_csv(os.path.join(args.dir,'Hugo_hg38.gz'),sep='\t').rename(columns=hugo_rename).dropna(subset=['Entrez']),
        'ncbi':pd.read_csv(os.path.join(args.dir,'NCBI_hg38.gz'),sep='\t')[['GeneID','Symbol','Aliases']].rename(columns=ncbi_rename).dropna(subset=['Entrez'])
    }

    mm_ref ={
        'ncbi': pd.read_csv(os.path.join(args.dir,'NCBI_mm10.gz'),sep='\t')[['GeneID','Symbol','Aliases']].rename(columns=ncbi_rename).dropna(subset=['Entrez']),
        'ensl':pd.read_csv(os.path.join(args.dir,'Ensl_mm10.gz')).rename(columns=ensl_rename).dropna(subset=['Entrez'])
    }

    match_stable = pd.read_csv(os.path.join(args.dir,'Ensl_mm10_hg38_match.gz'))
    match_stable.columns = ['hgENSG','hgEntrez','mmENSG','mmEntrez']

    # Stable reference make sure that the stable version are one-to-one map
    ## Concadinate all source
    hg_ref['hugo']['SymbolAlias'] = hg_ref['hugo']['Previous'] +', '+ hg_ref['hugo']['SymbolAlias']
    hg_ref['hugo'] = hg_ref['hugo'].drop('Previous',axis=1)
    hg_ref = pd.concat(hg_ref.values(),axis=0,ignore_index=True,sort=False)
    mm_ref = pd.concat(mm_ref.values(),axis=0,ignore_index=True,sort=False)

    ref = {
        'hg':hg_ref,
        'mm':mm_ref,
        'match':match_stable
    }
    constructRefDb(ref=ref,out=args.output)
    # Generate Test Case
    genTestCase(ref=ref,out=args.testOutput)




if __name__ == '__main__':
    main()
