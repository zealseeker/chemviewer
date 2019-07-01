from flask import render_template
from flask import Flask
from flask import jsonify
from flask import request
from werkzeug.utils import secure_filename
import random
import pandas as pd
from six import BytesIO
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from chemviewer import __version__

app = Flask(__name__)
FILEMAXSIZE = 100

def index(type='standalone', **kwargs):
    return render_template('index.html', standalone=(type=='standalone'),
                           version=__version__, **kwargs)

def upload(request):
    if request.method == 'POST':
        if request.files:
            f = request.files['file']
            suffix = secure_filename(f.filename).split('.')[-1]
            try:
                df = {'csv': pd.read_csv, 'xlsx': pd.read_excel,
                      'txt': pd.read_table, 'xls': pd.read_excel}[suffix](f)
            except:
                return jsonify('Err')
        elif 'data' in request.form:
            data = request.form['data']
            strIO = BytesIO()
            strIO.write(data.encode('utf8'))
            strIO.seek(0)
            df = pd.read_table(strIO, header=None, sep='[\t|\s+,]', engine='python')
            df = df[[0,1]]
            df.columns = ['SMILES', 'Name'][:len(df.columns)]
        for col in df.columns:
            if str(col).lower() in set(['smiles', 'smi', 'cannonical_smiles']):
                smi_col = col
                break
        else:
            smi_col = None
        return jsonify({'tableData': ana_table(df, smi_col),'queryNo':generate_key()})
    else:
        print('Not POST')
        return jsonify('Error')

def generate_key():
    return ''.join(random.sample('zyxwvutsrqponmlkjihgfedcba',8))

def draw_without_label(m,highlightAtoms=None,size=(200,200)):
    rdDepictor.Compute2DCoords(m)
    try:
        Chem.Kekulize(m)
    except:
        pass
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    drawer.DrawMolecule(m,highlightAtoms = highlightAtoms)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    return svg

def ana_table(df, col='SMILES'):
    # TODO find SMILES column
    data = df.to_dict(orient='records')
    if col:
        for item in data:
            # It's not the best way to use try
            try:
                item['_svg'] = draw_without_label(Chem.MolFromSmiles(item[col]),size=(150,150))
            except:
                item['_svg'] = ''
    return data

@app.route('/upload', methods=['POST', 'GET'])
def upload_wrapper():
    return upload(request)
