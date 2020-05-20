#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author            : Jingxin Fu <jingxinfu.tj@gmail.com>
# Date              : 18/02/2020
# Last Modified Date: 20/05/2020
# Last Modified By  : Jingxin Fu <jingxinfu.tj@gmail.com>

import os
import sys
import pandas as pd
from pandas.testing import assert_frame_equal
import unittest
from Biopyutils.Exprsn import bulkRNASeq

TestData = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data')


class TestbulkRNASeq(unittest.TestCase):
    def setUp(self):
        '''The original TPM and counts are from TCGA LUAD samples'''
        self.Tpm = pd.read_csv(os.path.join(TestData,'Exprsn','bulkRNASeq_tpm.csv'),index_col=0)
        self.counts = pd.read_csv(os.path.join(TestData,'Exprsn','bulkRNASeq_counts.csv'),index_col=0)
        self.meta = pd.read_csv(os.path.join(TestData,'Exprsn','bulkRNASeq_meta.csv'),index_col=0)
        # Parameters
        self.geneset = pd.read_csv(os.path.join(TestData,'Exprsn','bulkRNASeq_geneset.txt'),sep='\t')

        self.contrast = []
        with open(os.path.join(TestData,'Exprsn','bulkRNASeq_contrast.csv'),'r') as f:
            for l in f:
                self.contrast.append(l.strip())

        # Results
        self.noBatchTpm = pd.read_csv(os.path.join(TestData,'Exprsn','bulkRNASeq_noBatchTpm.csv'),index_col=0)
        self.DEG = pd.read_csv(os.path.join(TestData,'Exprsn','bulkRNASeq_DEG.TumorMinusNormal'),index_col=0)
        self.infiltras = pd.read_csv(os.path.join(TestData,'Exprsn','bulkRNASeq_infiltras.csv'),index_col=0)
        self.signature = pd.read_csv(os.path.join(TestData,'Exprsn','bulkRNASeq_signature.csv'),index_col=0)

    def test_removeBatch(self):
        testIns = bulkRNASeq(tpm=self.Tpm,counts=self.counts,meta=self.meta)
        testIns.removeBatch(fig_out='test_Batch')
        assert_frame_equal(testIns.tpm,self.noBatchTpm)

    def test_getDEG(self):
       testIns = bulkRNASeq(tpm=self.Tpm,counts=self.counts,meta=self.meta)
       testIns.getDEG(group_name='Source',contrast=self.contrast,force_limma=True)
       assert_frame_equal(testIns.DEG['TumorMinusNormal'],self.DEG)

       # Test the linear one
       testIns.getDEG(group_name='Source',contrast=self.contrast)

    def test_getSignature(self):
       testIns = bulkRNASeq(tpm=self.Tpm,counts=self.counts,meta=self.meta)
       testIns.getSignature(geneset=self.geneset)
       assert_frame_equal(testIns.signature,self.signature)

    def test_getInfiltras(self):
        testIns = bulkRNASeq(tpm=self.Tpm,counts=self.counts,meta=self.meta)
        testIns.getInfiltras()
        assert_frame_equal(testIns.infiltras, self.infiltras)

