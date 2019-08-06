from flask import jsonify
from flask import render_template
from flask import redirect, url_for
from flask import request
from flask import abort
from flask import Response
from chemviewer.core import index
from chemviewer.core import upload
from chemviewer.core import app
from chemviewer.core import ana_table
from chemviewer.core import generate_key
import os
import pandas as pd

@app.route('/')
def entrence():
    return index(type='standalone')

@app.route('/dir')
def list_dir():
    # TODO make it safer
    path = request.args.get('path')
    filenames = os.listdir(path)
    res = {'dir':[],'file':[]}
    for each in filenames:
        if os.path.isdir(os.path.join(path, each)):
            res['dir'].append(each)
        else:
            if each.split('.')[-1] in set(['csv','xlsx','txt','xls']):
                res['file'].append(each)
    return jsonify(res)

@app.route('/file')
def select_file():
    path = request.args.get('path')
    suffix = path.split('.')[-1]
    if suffix not in set(['csv', 'xlsx', 'txt', 'xls']):
        return abort(404)
    try:
        df = {'csv': pd.read_csv, 'xlsx': pd.read_excel,
              'txt': pd.read_table, 'xls': pd.read_excel}[suffix](path)
    except Exception as e:
        abort(Response('Cannot parse the file'))
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
                    'queryNo':generate_key(),
                    'type': type})


if __name__ == '__main__':
    app.run(debug=True)
