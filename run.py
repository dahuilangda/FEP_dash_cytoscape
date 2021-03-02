import dash
import dash_cytoscape as cyto
import dash_html_components as html
from dash.dependencies import Input, Output
import dash_core_components as dcc
# enable svg export
cyto.load_extra_layouts()

import networkx as nx
from networkx.readwrite import cytoscape_data
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)
import os
from urllib import parse

# Get FEP values
from perses.analysis.load_simulations import Simulation
import pickle
import sys
import os
import numpy as np

def smi2svg(smi, ref_smi):
    mol = Chem.MolFromSmiles(smi)
    ref_mol = Chem.MolFromSmiles(ref_smi)
    try:
        Chem.rdmolops.Kekulize(mol)
        Chem.rdmolops.Kekulize(ref_mol)
    except:
        pass
    
    mcs = rdFMCS.FindMCS([mol, ref_mol])
    mcs_query = Chem.MolFromSmarts(mcs.smartsString)
    AllChem.Compute2DCoords(mcs_query)
    AllChem.GenerateDepictionMatching2DStructure(mol,mcs_query)
    match = mol.GetSubstructMatch(mcs_query)
#     unmatch = list(map(lambda x: x[0]-x[1], zip(list(range(mol.GetNumAtoms())), list(match))))
    unmatch = []
    for x in range(mol.GetNumAtoms()):
        if not x in match:
            unmatch.append(x)
    
    drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
    AllChem.Compute2DCoords(mol)
    opts = drawer.drawOptions()
    opts.clearBackground=False
    drawer.DrawMolecule(mol, highlightAtoms=unmatch)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")
    return svg

FEP_list = []
for root, dirs, files in os.walk("mm_data", topdown=False):
    for name in files:
        if 'pi' == name[-2:]:
            FEP_list.append(os.path.join(root, name))

def extractDeltaG(FEP_file):
    try:
        results = np.load(FEP_file, allow_pickle=True)
        dg = round(-1*(results.bindingdg/results.bindingdg.unit),2)
        ddg = round(results.bindingddg/results.bindingddg.unit,2)
    except:
        dg = 999
        ddg = 999
    return dg, ddg


dg_dict = {}
ddg_dict = {}
for FEP_file in FEP_list: 
    dg, ddg = extractDeltaG(FEP_file)
    dg_dict[FEP_file] = dg
    ddg_dict[FEP_file] = ddg


def smi2image(smi, ref_smi):
    svg_string = smi2svg(smi, ref_smi)
    impath = 'data:image/svg+xml;charset=utf-8,' + parse.quote(svg_string, safe="")
    return impath

g = nx.Graph()
# add node

sdf = Chem.SDMolSupplier('mm_data/inputs/Tyk2_ligands_shifted.sdf')
mols = [m for m in sdf]
for mol in mols:
    Chem.rdmolops.Kekulize(mol)

# 
# name = mols[0].GetProp('_Name')
# g.add_node(ref_smi, img=smi2image(ref_smi, ref_smi), name=name)
# for mol, dg, ddg in zip(mols[1:], dg_list, ddg_list):

ref_smi = Chem.MolToSmiles(mols[0])
for mol in mols:
    smi = Chem.MolToSmiles(mol)
    name = mol.GetProp('_Name')
    g.add_node(smi, img=smi2image(smi, ref_smi), name=name)

for FEP_file in FEP_list:
    i, j = FEP_file.split('/')[1].split('.')[0][3:].split('to')
    i, j = int(i), int(j)

    smi_i = Chem.MolToSmiles(mols[i])
    smi_j = Chem.MolToSmiles(mols[j])
    dg = dg_dict[FEP_file]
    ddg = ddg_dict[FEP_file]

    edge_label = f'{dg}±{ddg}'
    if edge_label == '999±999':
        edge_label = 'nan'
    # name = mols[j].GetProp('_Name')
    # g.add_node(smi, img=smi2image(smi, ref_smi), name=name)
    if dg == 999:
        weight = 5
        color = 'red'
    elif dg >= -1:
        weight = 5
        color = 'black'
    else:
        weight = -1*dg*5
        color = 'lightblue'
    g.add_edge(smi_j, smi_i, name=edge_label, weight=weight, color=color)

    app = dash.Dash(__name__)
server = app.server

simple_elements = cytoscape_data(g)['elements']

styles = {
    'output': {
        'overflow-y': 'scroll',
        'overflow-wrap': 'break-word',
        'height': 'calc(100% - 25px)',
        'border': 'thin lightgrey solid'
    },
    'tab': {'height': 'calc(98vh - 115px)'}
}

stylesheet=[{'selector': 'node',
  'css': {'shape': 'rectangle',
   'width': '2000',
   'height': '1200',
   'border-color': 'rgb(0,0,0)',
   'border-opacity': 1.0,
   'border-width': 0.0,
   'color': '#4579e8',
   'background-image': 'data(img)',
   'background-fit': 'contain',
   'content': 'data(name)',
   'font-size': '150px',
   'background-opacity': '0'}},
 {'selector': 'edge',
  'style': {
   'font-size': '120px',
   'line-color': 'data(color)',
   'width': 'data(weight)',
#   'target-arrow-color': 'black',
   'target-arrow-shape': 'triangle',
   'arrow-scale': '5/data(weight)',
   'curve-style': 'bezier',
   'label': 'data(name)'}}]

app.layout = html.Div([
    html.Div([
    cyto.Cytoscape(
        id='cytoscape',
        elements=simple_elements,
        stylesheet=stylesheet,
        style={
            'height': '95vh',
            'width': 'calc(100% - 200px)',
            'float': 'left'}
    )
    ]),

    html.Div(className='four columns', children=[
        dcc.Tabs(id='tabs-image-export', children=[
            dcc.Tab(label='generate jpg', value='jpg'),
            dcc.Tab(label='generate png', value='png')
        ]),
        html.Div(style=styles['tab'], children=[
            html.Div(
                id='image-text',
                children='image data will appear here',
                style=styles['output']
            )
        ]),
        html.Div('Download graph:'),
        html.Button("as jpg", id="btn-get-jpg"),
        html.Button("as png", id="btn-get-png"),
        html.Button("as svg", id="btn-get-svg")
    ])
])

@app.callback(
    Output('image-text', 'children'),
    Input('cytoscape', 'imageData'),
    )
def put_image_string(data):
    return data


@app.callback(
    Output("cytoscape", "generateImage"),
    [
        Input('tabs-image-export', 'value'),
        Input("btn-get-jpg", "n_clicks"),
        Input("btn-get-png", "n_clicks"),
        Input("btn-get-svg", "n_clicks"),
    ])

def get_image(tab, get_jpg_clicks, get_png_clicks, get_svg_clicks):

    # File type to output of 'svg, 'png', 'jpg', or 'jpeg' (alias of 'jpg')
    ftype = tab

    # 'store': Stores the image data in 'imageData' !only jpg/png are supported
    # 'download'`: Downloads the image as a file with all data handling
    # 'both'`: Stores image data and downloads image as file.
    action = 'store'

    ctx = dash.callback_context
    if ctx.triggered:
        input_id = ctx.triggered[0]["prop_id"].split(".")[0]

        if input_id != "tabs":
            action = "download"
            ftype = input_id.split("-")[-1]

    return {
        'type': ftype,
        'action': action
        }

if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')
