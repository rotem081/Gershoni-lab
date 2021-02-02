import plotly.express as px
import pandas as pd
import flask
from flask import request, jsonify
import plotly.express as px

app = flask.Flask(__name__)
app.config["DEBUG"] = True


def remove_space(lst):
    lst = lst[0].split('|')
    x=' '.join(lst)
    return x.split()

def gene_name(lst):
    return lst[3]

def inbreed(x):
    return round(x,3)

def effect(lst):
    return lst[2]

# Create some test data for our catalog in the form of a list of dictionaries.
df_exome = pd.read_pickle('df_exome')
df_exome['mutation'] =df_exome['mutation'].apply(remove_space)
df_exome['gene'] =df_exome['mutation'].apply(gene_name)
df_exome['effect'] =df_exome['mutation'].apply(effect)
df_exome['Inbreeding'] =df_exome['Inbreeding'].apply(inbreed)
df_exome.drop('mutation',axis='columns', inplace=True)



@app.route('/', methods=['GET'])
def home():
    return '''<h1>Distant Reading Archive</h1>
        <p>A prototype API for distant reading of science fiction novels.</p>'''

@app.route('/api/dfbovine/search/genes', methods=['GET'])
def api_gene():
    # Check if an ID was provided as part of the URL.
    # If ID is provided, assign it to a variable.
    # If no ID is provided, display an error in the browser.
    if 'gene' in request.args:
        x = df_exome.where(df_exome['gene'] == request.args['gene']).dropna()
        return jsonify(x.to_dict())
    
    if 'pos' in request.args:
        request_chr, request_pos = request.args['pos'].split(':')
        x= df_exome.where((df_exome['Position']==int(request_pos)) & (df_exome['Chromosome']==request_chr)).dropna()
        return jsonify(x.to_dict())
    
    if 'allgene' in request.args:
        lst = df_exome['gene'].unique()
        lst=pd.DataFrame(lst)
        return jsonify(lst.to_dict())
    
    if 'allpos' in request.args:
        lst = df_exome['Position'].unique()
        lst=pd.DataFrame(lst)
        return jsonify(lst.to_dict())
    
    if 'allchr' in request.args:
        lst = df_exome['Chromosome'].unique()
        lst=pd.DataFrame(lst)
        return jsonify(lst.to_dict())
    
    if 'columns' in request.args:
        lst = df_exome.columns
        lst=pd.DataFrame(lst)
        return jsonify(lst.to_dict())
    
    return ''

app.run()
