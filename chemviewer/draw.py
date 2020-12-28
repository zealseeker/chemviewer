from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
from rdkit.Geometry import Point2D
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ChemicalReaction

svg_temp = """<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='{}px' height='{}px' viewBox='0 0 {} {}'>{}</svg>"""


def mol_to_svg(mol, size=(150, 150), *args, **kwargs):
    """Draw a molecule into svg (without header)

    Args:
        mol: RDKit molecule
        size: Width and height of the size
        *args, **kwargs: Used in `rdMolDraw2D.DrawMolecule`

    Returns:
        str: An SVG string without header
    """
    d = rdMolDraw2D.MolDraw2DSVG(*size)
    if mol == '->':
        d.DrawArrow(Point2D(0, size[1]/2), Point2D(size[0], size[1]/2))
    elif mol == '+':
        d.SetFontSize(30)
        d.DrawString('+', Point2D(size[0]/2, size[1]/2))
    else:
        Chem.SanitizeMol(mol)
        m = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
        d.DrawMolecule(m, *args, **kwargs)
    s = d.GetDrawingText().split('<!-- END OF HEADER -->\n')[1]
    return s

def extract_from_svg(svg: str):
    "Remove header of the original SVG string"
    if s[-6:] == '</svg>':
        s = s[:-6]
    s = s.partition('<!-- END OF HEADER -->\n')[2]
    return s


def group_svg(svg, offset: float = 0):
    """Wrap <g> tag for the svg
    Args:
        offset: The width offset in the whole reaction

    Returns:
        str: Grouped SVG string (without header)
    """
    if offset == 0:
        first_part = ''
    else:
        first_part = ' transform="translate({})"'.format(offset)
    return """<g{}>{}</g>""".format(first_part, svg)


def draw_reaction(r: ChemicalReaction, sub_size=(150, 150)):
    """A new implementation of reaction draw

    Args:
        r: The chemical reaction object of RDKit, e.g. from `AllChem.ReactionFromSmarts`
        sub_size: The size of each molecule in the reaction

    Returns:
        str: An SVG string with header (can be renderred directly)
    """
    sw = sub_size[0]
    sh = sub_size[1]
    svg = ''
    plus_w = 20
    offset = 0
    for i, reactant in enumerate(r.GetReactants()):
        if i != 0:
            svg += group_svg(mol_to_svg('+', size=(plus_w, sh)), offset=offset)
            offset += plus_w
        svg += group_svg(mol_to_svg(reactant, size=sub_size), offset=offset)
        offset += sw
    svg += group_svg(mol_to_svg('->', size=(50, sh)), offset=offset)
    offset += 50
    for i, product in enumerate(r.GetProducts()):
        if i != 0:
            svg += group_svg(mol_to_svg('+', size=(plus_w, sh)), offset=offset)
            offset += plus_w
        svg += group_svg(mol_to_svg(product, size=sub_size), offset=offset)
        offset += sw
    width = offset
    height = sh
    return svg_temp.format(width, height, width, height, svg)

def draw_molecule(m, highlightAtoms: list, size=(200, 200)) -> str:
    """Draw molecule with highlight atoms

    Args:
        m: RDKit molecule object
        size: The size of the SVG figure

    Returns:
        str: An SVG string with header
    """
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
