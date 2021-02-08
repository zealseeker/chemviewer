import os
import pandas as pd
from flask import jsonify
from flask import render_template
from flask import request
from flask import abort
from flask import Response
from chemviewer.core import index
from chemviewer.core import app
from chemviewer.core import ana_table
from chemviewer.core import generate_key
from chemviewer.analysis import Analysis
from chemviewer.common import progress
from chemviewer.common import logger
from chemviewer import __version__


@app.route('/')
def entrence():
    return index(type='standalone')


@app.route('/dir')
def list_dir():
    # TODO make it safer
    path = request.args.get('path')
    filenames = os.listdir(path)
    res = {'dir': [], 'file': []}
    for each in filenames:
        if os.path.isdir(os.path.join(path, each)):
            res['dir'].append(each)
        else:
            if each.split('.')[-1] in set(['csv', 'xlsx', 'txt', 'xls', 'cv']):
                res['file'].append(each)
    return jsonify(res)


@app.route('/file')
def select_file():
    path = request.args.get('path')
    suffix = path.split('.')[-1]
    if suffix == 'cv':
        jobid = path.split('/')[-1].split('.')[0]
        analysis = Analysis(jobid)
        return jsonify(analysis.profile)
    if suffix not in set(['csv', 'xlsx', 'txt', 'xls']):
        return abort(404)
    try:
        df = {'csv': pd.read_csv, 'xlsx': pd.read_excel,
              'txt': pd.read_table, 'xls': pd.read_excel}[suffix](path)
    except Exception as e:
        abort(Response('Cannot parse the file'))
    load_mode = request.args.get('mode')
    if load_mode is None:
        type = None
        for col in df.columns:
            if str(col).lower() in set(['smiles', 'smi', 'cannonical_smiles']):
                type = 'compound'
            if str(col).lower() in set(['reaction']):
                type = 'reaction'
            if type:
                smi_col = col
                break
        else:
            smi_col = None
        return jsonify({'tableData': ana_table(df, smi_col, type),
                        'queryNo': generate_key(),
                        'type': type})
    elif load_mode == 'pre':
        res = {
            'columns': df.columns.to_list(),
            'length': len(df)
        }
        return jsonify(res)


@app.route('/analysis')
def analysis_view():
    return render_template('analysis.html', standalone=(type == 'standalone'),
                           version=__version__)


@app.route('/analysis/calcfp')
def calc_fp():
    path = request.args.get('path')
    fp = request.args.get('fp')
    jobid = request.args.get('jobid')
    smi_col = request.args.get('smi_col')
    suffix = path.split('.')[-1]
    if jobid:
        analysis = Analysis(jobid)
    elif path:
        if suffix not in set(['csv', 'xlsx', 'txt', 'xls']):
            return abort(404)
        try:
            df = {'csv': pd.read_csv, 'xlsx': pd.read_excel,
                'txt': pd.read_table, 'xls': pd.read_excel}[suffix](path)
        except Exception as e:
            return jsonify({'message': 'Cannot parse the file', 'status': 500})
        analysis = Analysis(jobid, df)
    else:
        abort(Response("Wrong parameters"))
    try:
        analysis.calc_fp(smi_col, fp)
    except KeyError as e:
        return jsonify({'message': str(e).strip("'"), 'status': 500})
    return jsonify({'jobid': analysis.jobid, 'status': 200})


@app.route('/analysis/pca')
def pca_fp():
    fp = request.args.get('fp')
    jobid = request.args.get('jobid')
    analysis = Analysis(jobid)
    if request.args.get('load'):
        res = analysis.load_pca()
    else:
        res = analysis.calc_pca(fp)
    return jsonify(res)


@app.route('/analysis/draw-structure')
def get_structure():
    smi_col = request.args.get('smi_col')
    jobid = request.args.get('jobid')
    analysis = Analysis(jobid)
    try:
        res = analysis.get_structures(smi_col)
    except KeyError as e:
        return jsonify({'message': str(e).strip("'"), 'status': 500})
    return jsonify({'jobid': analysis.jobid, 'status': 200})


@app.route('/analysis/fetch-structure')
def fetch_structure():
    jobid = request.args.get('jobid')
    analysis = Analysis(jobid)
    analysis.load('svg')
    return jsonify(analysis.svgs)


@app.route('/progress')
def get_progress():
    jobid = request.args.get('jobid')
    if jobid in progress:
        return str(progress[jobid])
    else:
        return '0'


if __name__ == '__main__':
    app.run(debug=True)
    logger.setLevel('DEBUG')
