import json
import os
import random
import string
from threading import Thread
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.decomposition import PCA
import pandas as pd
from chemviewer.common import progress
from chemviewer.draw import draw_molecule
from chemviewer.common import logger
from chemviewer.utils import NumpyEncoder


def calc_fingerprint(obj: 'Analysis', smi_col, fp='ECFP4'):
    fps = []
    indices = []
    df = obj.df
    l = len(df)
    crt_i = 0
    logger.debug("Calculating fingerprints...")
    for i, smi in enumerate(df[smi_col]):
        try:
            m = Chem.MolFromSmiles(smi)
            bit = AllChem.GetMorganFingerprintAsBitVect(
                m, 2, nBits=2048).ToBitString()
            fp_output = list(map(int, list(bit)))
            idx = df.index[i]
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


def draw_structuers(obj, smi_col):
    svgs = {}
    crt_i = 0
    l = len(obj.df)
    for i, smi in enumerate(obj.df[smi_col]):
        try:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            svg = draw_molecule(m)
            svgs[i] = svg
        except Exception as e:
            logger.debug(e)
        percentile = i * 100 // l
        if percentile > crt_i:
            crt_i = percentile
            progress[obj.jobid] = crt_i
    obj.svgs = svgs
    logger.debug('Drawed!')
    obj.save()
    progress[obj.jobid] = 100
    return svgs


class Analysis:

    def __init__(self, jobid=None, df=None):
        if jobid is None or jobid == 'undefined':
            jobid = ''.join(random.choices(
                string.ascii_uppercase + string.digits, k=6))
        self.jobid = jobid
        self.smi_col = None
        self.df = df
        self.fps = {}
        self.svgs = {}
        self.df_path = '.tmp/{}-df.csv'.format(self.jobid)
        self.fp_path = '.tmp/{}-fp.json'.format(self.jobid)
        self.profile_path = '.tmp/{}.cv'.format(self.jobid)
        self.pca_path = '.tmp/{}-pca.json'.format(self.jobid)
        self.svg_path = '.tmp/{}-svg.json'.format(self.jobid)
        progress[jobid] = 0

    def calc_fp(self, smi_col=None, fp='ECFP4'):
        if self.df is None:
            self.load('df')
        if smi_col is None:
            smi_col = self.smi_col
            self.smi_col = smi_col
        if smi_col is None or smi_col=='null' or smi_col not in self.df:
            raise KeyError(f'SMILES columns {smi_col} is not available')
        x = Thread(target=calc_fingerprint, args=(self, smi_col, fp))
        x.start()
        return None

    def load(self, property='all'):
        if property in ['all', 'df']:
            if os.path.exists('.tmp/{}-df.csv'.format(self.jobid)):
                self.df = pd.read_csv(self.df_path, index_col=0)
        if property in ['all', 'svg']:
            with open(self.svg_path) as fp:
                self.svgs = json.load(fp)
        if property == 'all' or property == 'fp':
            if os.path.exists(self.fp_path):
                with open(self.fp_path) as fp:
                    self.fps = json.load(fp)

    def save(self, property='all'):
        if not os.path.exists('.tmp'):
            os.mkdir('.tmp')
        if property == 'all' or property == 'df':
            self.df.to_csv(self.df_path)
        if property in ['all', 'fp']:
            with open(self.fp_path, 'w') as fp:
                json.dump(self.fps, fp, cls=NumpyEncoder)
        if property in ['all', 'svg']:
            with open(self.svg_path, 'w') as fp:
                json.dump(self.svgs, fp)
        if property in ['all', 'profile']:
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
                'jobid': self.jobid,
                'smi_col': self.smi_col
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
        res = self.df.loc[idx, ['PCA_0', 'PCA_1']
                          ].reset_index().to_dict(orient='records')
        with open(self.pca_path, 'w') as fp:
            json.dump(res, fp)
        return res

    def load_pca(self):
        with open(self.pca_path) as fp:
            res = json.load(fp)
        return res

    def get_structures(self, smi_col=None):
        if self.df is None:
            self.load('df')
        if smi_col is None:
            smi_col = self.smi_col
            self.smi_col = smi_col
        if smi_col is None:
            raise ValueError('No smi_col')
        if smi_col not in self.df:
            raise KeyError(f'SMILES column error: {smi_col} is not in the data')
        x = Thread(target=draw_structuers, args=(self, smi_col))
        x.start()
        return None
