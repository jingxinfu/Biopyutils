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

import os
import sys
import pandas as pd
from pandas.testing import assert_frame_equal,assert_series_equal

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import unittest
from Biopyutils import Comm

TestData = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data')


class TestComm(unittest.TestCase):
    def setUp(self):
        species_list = ['hg','mm']
        source_list = ['ENST','ENSG','Symbol','Entrez']
        map_list = ['ENSG','Symbol','Entrez','Name']
        match_id_list = ['ENSG','Entrez']
        self.idConvertCase= [
                (sp,y,os.path.join(TestData,'Comm',sp+'_'+x+'_to_'+ y))
                for sp in species_list for x in source_list for y in map_list if x != y
                ]
        self.speciesCovertCase = [
                (x,y,os.path.join(TestData,'Comm',gid+'_'+x+'_to_'+ y))
                for gid in match_id_list for x in species_list for y in species_list if x != y
                ]
        self.infoMissingCase = [
                os.path.join(TestData,'Comm',x)  for x in os.listdir(os.path.join(TestData,'Comm')) if x.startswith('infoMissing')
                ]
    def test_infoMissing(self):
        import logging
        logging.basicConfig(level=logging.DEBUG)
        logger = logging.getLogger('[InfoMissing]')
        for f in self.infoMissingCase:
            logger.warning('----------Test %s' % f)
            in_df = pd.read_csv(f,index_col=0,sep='\t')
            Comm.idConvert(df=in_df,species='hg',map_id='Entrez',logger=logger,show_num=30)
            Comm.idConvert(df=in_df,species='hg',map_id='ENSG',logger=logger,show_num=30)




    def test_idConvert(self):
        for species,map_id,idcase in self.idConvertCase:
            for test_file in ['OneToOne','OneToMultiple','MultipleToOne']:
                file_prefix = idcase+'_'+ test_file
                in_df = pd.read_csv(file_prefix+'.input',index_col=0)
                if in_df.shape[0] == 0:
                    continue
                out_df = pd.read_csv(file_prefix+'.output',index_col=0)

                # DataFrame Case
                id_result  = Comm.idConvert(in_df,species=species,map_id=map_id)
                print(file_prefix)
                print(id_result)
                print(out_df)
                assert_frame_equal(id_result.reset_index(),out_df.reset_index())

                # Series Case
                in_series = in_df.squeeze()
                id_result  = Comm.idConvert(df=in_series,species=species,map_id=map_id).to_frame().reset_index()
                assert_frame_equal(id_result,out_df.reset_index())

                ## list Case
                in_list = in_df.index.tolist()
                out_list = out_df.index.tolist()
                id_result  = Comm.idConvert(df=in_list,species=species,map_id=map_id)
                self.assertEqual(id_result,out_list)

    def test_speciesConvert(self):
        for from_species,to_species,idcase in self.speciesCovertCase:
            for test_file in ['MultipleToOne','OneToOne','MultipleToOne']:
                file_prefix = idcase+'_'+ test_file
                in_df = pd.read_csv(file_prefix+'.input',index_col=0)
                if in_df.shape[0] == 0:
                    continue
                out_df = pd.read_csv(file_prefix+'.output',index_col=0)

                # DataFrame Case
                id_result  = Comm.speciesConvert(df=in_df,from_species=from_species,to_species=to_species)
                assert_frame_equal(id_result.reset_index(),out_df.reset_index())

                ## Series Case
                in_series = in_df.squeeze()
                id_result  = Comm.speciesConvert(df=in_series,from_species=from_species,to_species=to_species).to_frame().reset_index()
                assert_frame_equal(id_result,out_df.reset_index())


                ## list Case
                in_list = in_df.index.tolist()
                out_list = out_df.index.tolist()
                id_result  = Comm.speciesConvert(df=in_list,from_species=from_species,to_species=to_species)
                self.assertEqual(id_result,out_list)




