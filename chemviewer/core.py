from flask import render_template
from flask import Flask
from flask import jsonify
from flask import request
from werkzeug.utils import secure_filename
import random
import pandas as pd
from six import BytesIO
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor, Draw
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
            df = pd.read_csv(strIO, header=None, sep='\t|\s+,', engine='python')
            if re.findall('>.*>', df[0][0]):
                df.columns = ['Reaction', 'Name'][:len(df.columns)]
            else:
                df.columns = ['SMILES', 'Name'][:len(df.columns)]
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
    else:
        print('Not POST')
        return jsonify('Error')

def generate_key():
    return ''.join(random.sample('zyxwvutsrqponmlkjihgfedcba',8))

def draw_without_label(m,highlightAtoms=None,size=(200,200)):
    rdDepictor.Compute2DCoords(m)
    try:
        Chem.Kekulize(m)
        mc = rdMolDraw2D.PrepareMolForDrawing(m, kekulize=True)
    except:
        mc = rdMolDraw2D.PrepareMolForDrawing(m, kekulize=False)
        pass
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    drawer.DrawMolecule(mc, highlightAtoms=highlightAtoms)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    return svg

def draw_reaction(rxn, subImgSize=(150, 150)):
    r = rxn
    width = subImgSize[0] * (rxn.GetNumReactantTemplates() + rxn.GetNumProductTemplates() + 1)
    d = rdMolDraw2D.MolDraw2DSVG(width, subImgSize[1])
    atlab = d.drawOptions().atomLabels
    gi = 0  # Global index
    for mol in r.GetReactants():
        rdDepictor.Compute2DCoords(mol)
        m = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
        i_gi_map = []
        for i, atom in enumerate(m.GetAtoms()):
            if atom.GetSymbol()=='C':
                atlab[gi] = str(atom.GetAtomMapNum())
            i_gi_map.append(gi)
            gi += 1
    covered_idx = set()
    for mol in r.GetProducts():
        rdDepictor.Compute2DCoords(mol)
        m = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
        i_gi_map = []
        for i, atom in enumerate(m.GetAtoms()):
            if atom.GetSymbol()=='C':
                atlab[gi] = str(atom.GetAtomMapNum())
            i_gi_map.append(gi)
            gi += 1
    d.DrawReaction(r)
    d.FinishDrawing()
    return d.GetDrawingText()

def ana_table(df, col='SMILES', type=None):
    data = df.to_dict(orient='records')
    if type == 'compound':
        for item in data:
            # TODO It's not the best way to use try
            try:
                item['_svg'] = draw_without_label(Chem.MolFromSmiles(item[col]),size=(200,200))
            except:
                item['_svg'] = ''
    elif type == 'reaction':
        for item in data:
            try:
                rxn = AllChem.ReactionFromSmarts(item[col], useSmiles=True)
                svg = draw_reaction(rxn, subImgSize=(150, 150))
                item['_svg'] = svg
            except Exception as e:
                item['_svg'] = ''
    return data

@app.route('/upload', methods=['POST', 'GET'])
def upload_wrapper():
    return upload(request)
