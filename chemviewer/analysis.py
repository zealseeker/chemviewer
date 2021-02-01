import json
import os
from threading import Thread
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
from chemviewer.common import progress
import logging

logger = logging.Logger('analysis')
logger.setLevel('DEBUG')
formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
ch = logging.StreamHandler()
ch.setLevel('DEBUG')
ch.setFormatter(formatter)
logger.addHandler(ch)

def calc_fingerprint(obj, smi_col, fp='ECFP4'):
    fps = []
    indices = []
    df = obj.df
    l = len(df)
    crt_i = 0
    logger.debug("Calculating fingerprints...")
    for i,smi in enumerate(df[smi_col]):
        try:
            m = Chem.MolFromSmiles(smi)
            bit = AllChem.GetMorganFingerprintAsBitVect(m,2, nBits=2048).ToBitString()
            fp_output = list(map(int,list(bit)))
            indices.append(df.index[i])
            fps.append(fp_output)
        except:
            pass
        percentile = i * 100 // l
        if percentile > crt_i:
            crt_i = percentile
            progress[obj.jobid] = crt_i
    obj.fps = {fp: {'matrix': fps, 'idx': indices}}
    logger.debug("Calculated!")
    obj.save()
    progress[obj.jobid] = 100
    return fps, indices

class Analysis:

    def __init__(self, jobid, df=None):
        self.jobid = jobid
        self.smi_col = None
        self.df = df
        self.fps = {}
        self.df_path = '.tmp/{}-df.csv'.format(self.jobid)
        self.fp_path = '.tmp/{}-fp.json'.format(self.jobid)
        self.profile_path = '.tmp/{}.cv'.format(self.jobid)
        self.pca_path = '.tmp/{}-pca.json'.format(self.jobid)
        progress[jobid] = 0
        # if df is None:
        #     self.load()
        # else:
        #     self.df = df

    def calc_fp(self, smi_col=None, fp='ECFP4'):
        if self.df is None:
            self.load('df')
        if smi_col is None:
            smi_col = self.smi_col
        if smi_col is None:
            raise ValueError('No smi_col')
        df = self.df
        x = Thread(target=calc_fingerprint, args=(self, smi_col, fp))
        x.start()
        return None
        

    def load(self, property='all'):
        if os.path.exists('.tmp/{}-df.csv'.format(self.jobid)):
            self.df = pd.read_csv(self.df_path, index_col=0)
        if property == 'all' or property == 'fp':
            if os.path.exists(self.fp_path):
                with open(self.fp_path) as fp:
                    self.fps = json.load(fp)

    def save(self, property='all'):
        if not os.path.exists('.tmp'):
            os.mkdir('.tmp')
        self.df.to_csv(self.df_path)
        with open(self.fp_path, 'w') as fp:
            json.dump(self.fps, fp)
        with open(self.profile_path, 'w') as fp:
            json.dump(self.profile, fp)

    @property
    def profile(self):
        if self.df is None:
            with open(self.profile_path) as fp:
                res = json.load(fp)
        else:
            res = {
                'columns': self.df.columns.to_list(),
                'length': len(self.df),
                'jobid': self.jobid
            }
        return res

    def calc_pca(self, fp='ECFP4'):
        self.load()
        pca = PCA(4)
        X = self.fps[fp]['matrix']
        idx = self.fps[fp]['idx']
        res = pca.fit_transform(X)
        for i in range(4):
            self.df.loc[idx, 'PCA_{}'.format(i)] = res[:, i]
            self.df['PCA_{}'.format(i)].fillna(0, inplace=True)
        self.save('df')
        res = self.df.loc[idx, ['PCA_0', 'PCA_1']].reset_index().to_dict(orient='records')
        with open(self.pca_path, 'w') as fp:
            json.dump(res, fp)
        return res

    def load_pca(self):
        with open(self.pca_path) as fp:
            res = json.load(fp)
        return res
