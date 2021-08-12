import xml.etree.ElementTree as ET
import pandas as pd 
import numpy as np
import matplotlib as plt
import plotly.express as px
import math
import gffutils
import ipywidgets as widgets

######################################
## enter file paths (user interaction)
######################################

def acceptFiles():
    instructions = widgets.HTML('Enter the .gff3 and .xml file links or paths into the spaces below. Hit "Enter" after entering each:')
    display(instructions)

    gff3Text = widgets.Text(
        value='',
        placeholder='Insert file path',
        description='.GFF3 File Path or Link',
        disabled=False,
        style= {'description_width': 'initial'}

    )
    display(gff3Text)

    def gff3Submit(placeHolder):
        gff3_url = gff3Text.value
        display(widgets.HTML(".gff3 file location entered."))
        return gff3_url
    gff3Text.on_submit(gff3Submit)

    xmlText = widgets.Text(
        value='',
        placeholder='Insert file path',
        description='.XML File Path or Link',
        disabled=False,
        style= {'description_width': 'initial'}

    )
    display(xmlText)

    def xmlSubmit(placeHolder):
        xml_url = xmlText.value
        display(widgets.HTML(".xml file location entered."))
        return xml_url
    xmlText.on_submit(xmlSubmit)

######################################
## parsing xml 
######################################

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
                                tem_d["Total"]=int(i[0].text)
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

def get_shortreadcounts(annot_xml, tagname):
    tree = ET.parse(annot_xml)
    root = tree.getroot()
    RNA_seq_data={}
    for c in root.iter():
        if c.tag == tagname:
            for c2 in c:
                if c2.tag == "AssemblyUnit":
                    for i in c2:
                        if i.tag == "RunStats":
                            if i.attrib["run"] == 'ALL': continue
                            RNA_seq_data[(i.attrib)["run"]] =i.attrib
                            for j in i:
                                RNA_seq_data[(i.attrib)["run"]][j.tag]=(j.text)
                            RNA_seq_data[(i.attrib)["run"]].pop('run', None)
    return (RNA_seq_data)

######################################
## functions generating plots 
######################################

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
def get_style(r=9, g=220, b=250, a=1):
    # RGB (9, 220, 250) corresponds to light blue
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

# Wrapper function which returns array of html for transcripts,genes, CDS total/longest/shortest summary respectively
def summ_dict(length_dict):
    # transcripts
    transcripts_html = html_total(length_dict['total transcript:'], "Transcripts")
    longkey, shortkey = get_keys(length_dict, 'transcript')
    transcripts_html += html_dict(longkey, shortkey, length_dict)
    
    # genes
    genes_html = html_total(length_dict['total genes:'], "Genes")
    longkey, shortkey = get_keys(length_dict, 'gene')
    genes_html += html_dict(longkey, shortkey, length_dict)
    
    # CDS
    CDS_html = html_total(length_dict['total CDS:'], "CDS")
    longkey, shortkey = get_keys(length_dict, 'CDS')
    CDS_html += html_dict(longkey, shortkey, length_dict)
    return [transcripts_html, genes_html, CDS_html]

# Returns html for header listing total
def html_total(tot, t):
    totstring = t + " (" + str(tot) + " total)"
    header = """
    <header>
        <h3>{total}</h3>
    </header>
    """.format(total = totstring)
    return header

# Returns which keys correspond to longest/shortest
def get_keys(length_dict, type):
    larg, sarg = 'longest '+type, 'shortest '+type
    for key in length_dict.keys():
        if larg in key:
            longkey = key
        if sarg in key:
            shortkey = key
    return longkey, shortkey

# Returns html for table summarizing longest/shortest
def html_dict(longkey, shortkey, length_dict):
    longname, shortname = longkey.split(": ")[1], shortkey.split(": ")[1]
    long, short = length_dict[longkey], length_dict[shortkey]
    html = """
    <table width=500px>
        <tr>
            <th>Longest</th>
            <th>Shortest</th>
        </tr>
        <tr>
            <td>{lname}</td>
            <td>{sname}</td>
        </tr>
        <tr>
            <td>{longest}</td>
            <td>{shortest}</td>
        </tr>
    </table>
    """.format(lname=longname, sname=shortname, longest=str(long)+" bp", shortest=str(short)+" bp")
    return html

######################################
## parsing gff3
######################################

def create_gff3_db(annot_gff, gff3_db):
    gffutils.create_db(
        annot_gff,
        dbfn=gff3_db,
        force=True,
        keep_order=True, 
        sort_attribute_values=True, 
        merge_strategy = "merge")

def get_length_dict(dbfn):
    db = gffutils.FeatureDB(dbfn)
    stats_dict = {}
    
    #get gene info from gff file
    gene_count = 0
    longest_gene = 0
    shortest_gene = 0
    l_gene = ""
    s_gene = ""
    exon_count = 0
    longest_exon = 0
    shortest_exon = 0
    l_exon = ""
    s_exon = ""
    
    for gene in db.features_of_type('gene'):
        gene_count += 1
        gene_name = gene.attributes.get('Name', None)[0]
        gene_range = (gene.end - gene.start)
        if(longest_gene < gene_range):
            longest_gene = gene_range
            l_gene = gene_name
        if(shortest_gene > gene_range or shortest_gene == 0):
            shortest_gene = gene_range
            s_gene = gene_name

         #get exon info   
        exons = list(db.children(gene, featuretype='exon'))
        for exon in exons:
            exon_count += 1
            exon_name = gene.attributes.get('ID', None)[0]
        exon_range = (exon.end - exon.start)
        if(longest_exon < exon_range):
            longest_exon = exon_range
            l_exon = exon_name
        if(shortest_exon > exon_range or shortest_exon == 0):
            shortest_exon = exon_range
            s_exon = exon_name

    #add gene info to dict
    l_gene = "longest gene: " + str(l_gene)
    s_gene = "shortest gene: " + str(s_gene)
    stats_dict[l_gene] = longest_gene
    stats_dict[s_gene] = shortest_gene
    stats_dict["total genes:"] = gene_count

    #add exon info to dict
    l_exon = "longest exon: " + str(l_exon)
    s_exon = "shortest exon: " + str(s_exon)
    stats_dict[l_exon] = longest_exon
    stats_dict[s_exon] = shortest_exon
    stats_dict["total exons:"] = exon_count

    #get transcript info from gff file
    transcript_count = 0
    longest_transcript = 0
    shortest_transcript = 0
    l_transcript = ""
    s_transcript = ""
    
    for transcript in db.features_of_type('transcript'):
        transcript_count += 1
        transcript_name = transcript.attributes.get('Name', None)[0]
        transcript_len = 0
        exons = list(db.children(transcript, featuretype='exon'))
        for exon in exons:
            exon_range = (exon.end - exon.start)
            transcript_len += exon_range
        if(longest_transcript < transcript_len):
            longest_transcript = transcript_len
            l_transcript = transcript_name
        if(shortest_transcript > transcript_len or shortest_transcript == 0):
            shortest_transcript = transcript_len
            s_transcript = transcript_name

    #add transcript info to dict
    l_transcript = "longest transcript: " + str(l_transcript)
    s_transcript = "shortest transcript: " + str(s_transcript)
    stats_dict[l_transcript] = longest_transcript
    stats_dict[s_transcript] = shortest_transcript
    stats_dict["total transcript:"] = transcript_count

    #get CDS info from gff file
    CDS_count = 0
    longest_CDS = 0
    shortest_CDS = 0
    l_CDS = ""
    s_CDS = ""
    prot_dict = {}
    for CDS in db.features_of_type('CDS'):
        CDS_count += 1
        CDS_name = None
        if CDS.attributes.get('Name', None):
            CDS_name = CDS.attributes.get('Name', None)[0]
        prot_acc = None
        if CDS.attributes.get('protein_id', None):
            prot_acc = CDS.attributes.get('protein_id', None)[0]
        CDS_range = (CDS.end - CDS.start)

        if prot_acc in prot_dict.keys():
            prot_dict[prot_acc] = prot_dict[prot_acc] + CDS_range
        else:
            prot_dict[prot_acc] = CDS_range

        if(longest_CDS < prot_dict[prot_acc]):
            longest_CDS = prot_dict[prot_acc]
            l_CDS = prot_acc
        if(shortest_CDS > prot_dict[prot_acc] or shortest_CDS == 0):
            shortest_CDS = prot_dict[prot_acc]
            s_CDS = prot_acc

    #add transcript info to dict
    l_CDS = "longest CDS: " + str(l_CDS)
    s_CDS = "shortest CDS: " + str(s_CDS)
    stats_dict[l_CDS] = longest_CDS
    stats_dict[s_CDS] = shortest_CDS
    stats_dict["total CDS:"] = CDS_count
    
    return stats_dict

######################################
## executing steps in notebook 
######################################

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

def tabulate_shortread_aligns(annot_xml):
    rnaseq_stats = get_shortreadcounts(annot_xml, 'RnaseqAlignReport')
    
    for k,v in rnaseq_stats.items():
        s1 = set(v.keys())
        s2 = set(['AlignedReadsPct', 'SplicedReadPct'])
        for i in s1 - s2:
            v.pop(i, None)
        for i in s2:
            v[i] = np.round(float(v[i])/100, 3)

    return pd.DataFrame.from_dict(rnaseq_stats, orient='index')

def tabulate_cage_aligns(annot_xml):
    rnaseq_stats = get_shortreadcounts(annot_xml, 'CageAlignReport')
    
    for k,v in rnaseq_stats.items():
        s1 = set(v.keys())
        s2 = set(['AlignedReadsPct', 'SplicedReadPct'])
        for i in s1 - s2:
            v.pop(i, None)
        for i in s2:
            v[i] = np.round(float(v[i])/100, 3)

    return pd.DataFrame.from_dict(rnaseq_stats, orient='index')
