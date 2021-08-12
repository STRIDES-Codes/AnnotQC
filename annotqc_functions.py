import xml.etree.ElementTree as ET
import pandas as pd 
import numpy as np
import matplotlib as plt
import plotly.express as px
import math
import ipywidgets as widgets

def get_featcounts(annot_xml):
    tree = ET.parse(annot_xml)
    root = tree.getroot()
    dictionary_attiribute_feature={}
    tem_d={}
    for c in root.iter():
        if c.tag == "FeatureCounts":
            for c2 in c:
                if c2.tag =="OverallCounts":
                    for i in c2:
                        for feature in i:
                            if len(feature.attrib.values())>=1:
                                tem_d[list(feature.attrib.values())[0]]=int(feature.text)
                        dictionary_attiribute_feature[list(i.attrib.values())[0]] = tem_d
                        tem_d={}
    return (dictionary_attiribute_feature)

def get_longreadcounts(annot_xml):
    tree = ET.parse(annot_xml)
    root = tree.getroot()
    long_read_reports={}
    for c in root.iter():
        if c.tag == "LongReadAlignReport":
            for c2 in c:
                if c2.tag =="Assembly":
                    for i in c2:
                        if i.attrib["run"] == 'ALL': continue
                        long_read_reports[(i.attrib)["run"]]=i.attrib
                        for feature in i:
                            if feature.tag =="PostFilter":
                                for j in feature:
                                    long_read_reports[(i.attrib)["run"]][j.tag]=j.text
                                long_read_reports[(i.attrib)["run"]].pop('run', None)
    return (long_read_reports)

def create_barplot(df, name, c):
    # Creating the bar plot 
    fig = px.bar(df, orientation='h',
                 height=250,
                 title='Counts')
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    fig.update_layout(title_text = f'{name}, total count: {c}', title_x=0.1, title_y = 0.75)
    fig.update_traces(hovertemplate= 0)
    fig.update_traces(hovertemplate= '%{x:.0f}<br>')
    fig.update_layout(

        paper_bgcolor = "white",
        plot_bgcolor = "white",

        legend = dict(
            yanchor = "top",
            y = 0.90,
            xanchor = "right",
            x = 1.3),

        showlegend=False,
        )
    fig.show()

# Defining width and background color of overall table
def get_table_id(width, bcolor):
    table_id = """<table width={w}px style=\"background-color:{c}\" border="1" class="center">""".format(
        w = width, c = bcolor
        )
    return table_id

# Returns html for table header (column names)
def get_table_header(df):
    table="<tr>"
    table += """<th>Rows</th>"""
    for col in df.columns:
        table += """<th>{c}</th>""".format(c=col)
    table += "</tr>"
    return table

def perc(x): # TODO may need to adjust depending on data
    x = nan_to_b(x)
    if isinstance(x, str):
        return "0%"
    return str(x*100)+"%"

def nan_to_b(x):
    if math.isnan(x):
        return ""
    return x

# Returns html for all data in table
def get_table_data(df):
    data = ""
    for index, row in df.iterrows():
        data += "<tr>"
        
        data += """<th scope="row">{row}</th>""".format(row = index)
        
        for i in range(len(df.columns)):
            data+="""
                <td style="text-align: center; background-size: {width} 100%" class="bar">{x}</td>
            """.format(width=perc(row[i]), x=nan_to_b(row[i]))
        data += "</tr>"
    data += "</table>" # end table
    return data

# Defines style of bars in cells
def get_style(r=9, g=220, b=250, a=1): # TODO adjust colors?
    style = "<style>"
    style += """
        td {
            background-image: linear-gradient(to right, rgba(9, 220, 250, 1) 0%, rgba(9, 220, 250, 1) 100%);
            background-repeat: no-repeat;
        }
    """
    style += "</style>"
    return style

# Wrapper function that calls all subfunctions and returns complete html for plot
def get_html(df, width=750, bcolor="#FFFFFF"):
    html = get_table_id(width, bcolor)
    html += get_table_header(df)
    html += get_table_data(df)
    html += get_style()
    return html



## functions for plots 
def create_genes_plot(annot_xml):
    feat_counts = get_featcounts(annot_xml)

    df = pd.DataFrame.from_dict(feat_counts)
    featsubtypes = ['other', 'non-transcribed pseudo', 'transcribed pseudo', 
                    'Ig TCR segment', 'non coding', 'protein coding']
    df = df.filter(items = featsubtypes, axis = 0) # filter by feat subtypes 
    df = df.filter(['genes'], axis = 1) # keep only gene counts 
    df = df.sort_values(by = 'genes', ascending=False) # sort by gene counts 
    df = df.transpose() # transpose for plot 

    c = int(df.sum(axis=1))
    display(df)
    create_barplot(df, "Genes", c)
    
def create_transcripts_plots(annot_xml):
    feat_counts = get_featcounts(annot_xml)
    df = pd.DataFrame.from_dict(feat_counts)

    featsubtypes = ['model RefSeq', 'known RefSeq']
    df = df.filter(items = featsubtypes, axis = 0)
    display(df.transpose().dropna())
    for i in ['mRNAs', 'non-coding RNAs', 'pseudo transcripts', 'CDSs']:
        dfi = df.filter([i], axis = 1) 
        dfi = dfi.transpose()
        c = int(dfi.sum(axis=1))
        create_barplot(dfi, i, c)
        
def tabulate_genes_attributes(annot_xml):
    feat_counts = get_featcounts(annot_xml)
    genesubtypes = ['other', 'non-transcribed pseudo', 'transcribed pseudo', 
                'Ig TCR segment', 'non coding', 'protein coding']
    genes_total = sum([v for k,v in feat_counts['genes'].items() if k in genesubtypes])
    d = {k:v for k,v in feat_counts['genes'].items() if k not in genesubtypes}
    display(pd.DataFrame.from_dict(d, orient='index', columns=['genes']).transpose())

    d = {k:np.round(v/genes_total, 4) for k,v in feat_counts['genes'].items() if k not in genesubtypes}
    
    return pd.DataFrame.from_dict(d, orient='index', columns=['genes']).transpose()

def tabulate_tx_attributes(annot_xml):
    feattypes = ['mRNAs', 'non-coding RNAs', 'pseudo transcripts', 'CDSs']

    d = get_featcounts(annot_xml)
    d.pop('genes', None)

    def add_total(d):
        for k,v in d.items():
            v['total'] = v['known RefSeq'] + v['model RefSeq']
            v.pop('known RefSeq', None)
            v.pop('model RefSeq', None)

    add_total(d)
    display(pd.DataFrame.from_dict(d, orient='index'))

    for k,v in d.items():
        c = v['total']
        v.pop('total')
        for i,j in v.items():
            v[i] = np.round(j/c, 4)

    return pd.DataFrame.from_dict(d, orient='index')

def tabulate_longread_aligns(annot_xml):
    longread_counts = get_longreadcounts(annot_xml)
    display(pd.DataFrame.from_dict(longread_counts, orient='index'))

    for k,v in longread_counts.items():
        s1 = set(v.keys())
        s2 = set(['AlignedReadsPct', 'PctCoverage', 'PctIdentity'])
        for i in s1 - s2:
            v.pop(i, None)
        for i in s2:
            v[i] = np.round(float(v[i])/100, 3)

    return pd.DataFrame.from_dict(longread_counts, orient='index')

