# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 18:30:58 2021

@author: Rotem
"""

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from skimage import io
import dash_table
import numpy as np
import plotly.express as px
import pandas as pd
import urllib
import base64
import os

global df

#read df from url and complete the path (allgene, column)
def get_from_url(parameters):
    ip = '10.26.88.106:5000'
    url =fr'http://{ip}/api/dfbovine/search/genes?{parameters}'
    f = urllib.request.urlopen(url)
    myfile = f.read()
    df_result = pd.read_json(myfile, orient='records')
    return df_result

gene_list = get_from_url('allgene=1')[:1000]
gene_list= [{'label': x, 'value': x} for x in gene_list[0]]

col_list = get_from_url('columns=1')
col_list= [{'label': x, 'value': x} for x in col_list[0]]

#list for the order of the column
new_col_list=[]
for i in col_list:
    new_col_list.append([x for x in i.values()][0])

col_dict=[{'id': x, 'name': x} for x in new_col_list]
df_empty = pd.DataFrame(columns=new_col_list)
df_empty= df_empty.append(['None'])
df = df_empty


def inbreed(x):
    return round(x,3)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
#app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app = dash.Dash(__name__,suppress_callback_exceptions=True,
    external_stylesheets=[dbc.themes.BOOTSTRAP]
)
 
#server = app.server
# HEX to rgb
colors = {
    'background': '#FFFFFF',
    'text': '#000000'
}
white_button_style = {'background-color': '#FFFFFF',
                      'color': 'black',
                      'height': '50px'}

red_button_style = {'background-color':'#0000FF',
                    'color': 'white',
                    'height': '50px'}

image_filename = 'count_2.png'# replace with your own image
encoded_image = base64.b64encode(open(image_filename, 'rb').read()).decode('ascii')

# the style arguments for the sidebar. We use position:fixed and a fixed width
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

# the styles for the main content position it to the right of the sidebar and
# add some padding.
CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

url_bar_and_content_div = html.Div([
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])

#let's start
sidebar = html.Div([
        html.H2("BExAD", className="display-4"),
        html.Hr(),
        html.P(
            "A Bovine Exomes Aggregate Database", className="lead"
        ),
        dbc.Nav(
            [
                dbc.NavLink("Home", href="/", active="exact"),
                dbc.NavLink("Find Gene", href="/page-1", active="exact"),
                #dbc.NavLink("Page 2", href="/page-2", active="exact"),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    style=SIDEBAR_STYLE,
)

content = html.Div(id="page-content", style=CONTENT_STYLE)

app.layout = html.Div([dcc.Location(id="url"), sidebar, content])
       
user_page= html.P([
    dbc.Label('Choose which information to show:', style = {'fontSize' : '20px',
                                                            'fontWeight': 'bold'}
              ),  
    dbc.Checklist(
        id="all-or-none",
        options=[{"label": "Select All", "value": "All"}],
        inline=True,
        value=[],     
    ), 
    dbc.Checklist(
        id='my-checklist',
        options=col_list,
        inline=True,
        value=[],
    ),
    html.Hr(),
    dbc.Label("Enter Gene Name:", style = {'fontSize' : '20px','fontWeight': 'bold'}),
    dcc.Dropdown(
        options= gene_list,
        id='gene_name', style= {'fontSize' : '20px','width': '400px'}),
    
    dbc.Button(children='Submit', id='button_1',style=white_button_style,
               color="secondary", className="mr-1"),
    html.Br(),
    dbc.Label("Or enter Gene Location: ",style = {'fontSize' : '20px','fontWeight': 'bold'}),
    html.Br(),
    dcc.Input(id='position', style = {'fontSize' : '20px','width': '400px'}),
    html.Br(),
    dbc.Button(children='Submit', id='button_2', style=white_button_style, 
               color="secondary", className="mr-1"),
    html.Hr(),
    html.Div(
        dash_table.DataTable(
        id='Table',
        data=df_empty.to_dict('records'),
        columns= col_dict,
        editable=True,
        page_current=0,
        page_size=15,
        style_cell= {'font_size': '16px','whiteSpace': 'normal','textAlign': 'left'},
        style_cell_conditional=[
            {
                'if': {'column_id': c},
                'textAlign': 'left',
            } for c in ['Date', 'Region']
        ],
        style_data_conditional=[
            {
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(248, 248, 248)'
            }
        ],
        style_header={
            'backgroundColor': 'rgb(230, 230, 230)',
            'fontWeight': 'bold'
        },
            page_action='custom'
            ),
    )
])
Home_page=html.Div(children=[
    html.H1(
        children='Bovine Exomes Aggregate Database (BExAD)',
        style={'padding': '30px', 'backgroundColor': '#C0C0C0',
            'textAlign': 'center', 'width': '100%',
            'color': colors['text']
        }
    ),
    
    dbc.Label("This DB contains exomes bovine variants of 269 cows from the NCBI web.\n"
           " You can search for allele frequency or other information\n"
           " by writing Gene name or, Chrmosome and Position.\n",
        style={'padding': '30px', 'textAlign': 'left',
               'width': '100%', 'fontWeight': 'bold', 'font_size': '26px'}),

    html.P(dcc.Markdown(''' In order to better interpret and annotate identified
           genetic variations, we established a bovine whole-exome aggregate database.
           We obtained Fastq files of bovine short reads from the Sequences
           Reads Archive (SRA) at the NCBI server, a publicly available repository
           of raw sequencing data ([https://www.ncbi.nlm.nih.gov/sra](/)).
           We used these sequences to build a database of the bovine exome genetic
           variation. Overall, we obtained over 250 bovine samples. Exome sequences
           were analyzed as detailed in the WGS analysis section
           (GitHub: Gershoni-lab/ WGS analysis). This pipeline produced a per-sample
           GVCF file, which contains the following information: variant position,
           reference and alternative allele, quality metrics and genotype.
           Next, VCF files were joined to form a multi-sample file, which included
           the population allele frequencies, the total alleles number per position
           and the inbreeding co-efficiency. The joint VCF was then annotated by
           Variant Effect Predictor (VEP). Then, filtering was done using a dedicated
           Python script that uses the VCF parser PyVCF, as detailed in our GitHub:
           [Gershoni-lab/WGS analysis/exome database](/)''',
           style={'padding': '20px', 'textAlign': 'left',
                  'width': '100%'}),),
        
    html.Img(src='data:image/png;base64,{}'.format(encoded_image),
             style={'textAlign':'left','width': '100%'})
             
])
    

# "complete" layout
app.validation_layout = html.Div([
    Home_page,
    user_page
])

@app.callback(
Output("page-content", "children"),
[Input("url", "pathname")])

def render_page_content(pathname):
    if pathname == "/":
        return Home_page
    elif pathname == "/page-1":
        return user_page
   # elif pathname == "/page-2":
    #    return html.P("Oh cool, this is page 2!")
    # If the user tries to reach a different page, return a 404 message
    return dbc.Jumbotron(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ]
    ) 


@app.callback(
Output("my-checklist", "value"),
[Input("all-or-none", "value")],
[State("my-checklist", "options")],
)
def select_all_none(all_selected, options):
    all_or_none = []
    all_or_none = [option["value"] for option in options if all_selected]
    return all_or_none



@app.callback(
    Output('Table', 'data'),
    Output('Table','page_count'),
    Output('Table','columns'),
    Input('Table', "page_current"),
    Input('Table', "page_size"),
    Input('button_1', 'n_clicks'),
    Input('button_2', 'n_clicks'),
    State('gene_name','value'),
    State('position','value'),
    [State('my-checklist', 'value')],
        
    )
def update_output_name(page_current,page_size, button1, 
                       button2, input_value_gene, input_value_loc,
                       selected_check):
    global df
    #which buttom was pressed on
    ctx = dash.callback_context 
    if ctx.triggered:
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
        if button_id=='button_1':
            if input_value_gene:
                df = get_from_url(f'gene={input_value_gene}')
            
        else: #button 2 pressed
            if input_value_loc:
                df = get_from_url(f'pos={input_value_loc}')
        
        if not df.equals(df_empty):
            df=df[new_col_list]
            df=df[selected_check].dropna()
            if 'Allele frequency' in df.columns:
                df['Allele frequency'] =df['Allele frequency'].apply(inbreed)
            if 'Inbreeding' in df.columns:
                df['Inbreeding'] =df['Inbreeding'].apply(inbreed)
            df['index'] = range(1, len(df) + 1)

    return df.iloc[
        page_current*page_size:(page_current+ 1)*page_size
    ].to_dict('records'), np.ceil(len(df)/page_size),[{'id': x, 'name': x} for x in selected_check]
    

#button style
@app.callback(Output('button_1', 'style'), [Input('button_1', 'n_clicks')])
def change_button_style(n_clicks):

    if n_clicks and n_clicks > 0:

        return red_button_style

    else:

        return white_button_style          
    
@app.callback(Output('button_2', 'style'), [Input('button_2', 'n_clicks')])
def change_button_style_2(n_clicks):

    if n_clicks and n_clicks > 0:

        return red_button_style

    else:

        return white_button_style          
          
        


#style={'columnCount': 2}))])
if __name__ == '__main__':
    app.run_server(debug=True,host='0.0.0.0')
    
