import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
from datetime import datetime
import base64
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import inch
import io
import math
import csv
import os

# Initialize CSV storage
def init_csv_storage():
    """Initialize or load the CSV storage file"""
    csv_file = "patient_data.csv"
    headers = [
        'timestamp', 'patient_id', 'age', 'gender',
        # Lab values
        'glucose', 'lactate', 'pyruvate', 'urea', 'creatinine', 'bun', 'ammonia',
        'alt', 'ast', 'bilirubin', 'vitamin_d', 'pth', 'tsh', 'cortisol',
        'insulin', 'c_peptide', 'estrogen', 'progesterone', 'testosterone',
        'apo_a', 'apo_b', 'total_cholesterol', 'ldl', 'hdl', 'triglycerides',
        'phosphorus', 'magnesium', 'uric_acid', 'anti_tpo', 'ana', 'anti_ccp',
        'crp', 'fgf_23', 'homocysteine', 'glutamine', 'glutamate', 'gsh', 'mda',
        # Analysis results
        'top_pathway', 'top_score', 'directly_affected_count', 'at_risk_count'
    ]
    
    if not os.path.exists(csv_file):
        with open(csv_file, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=headers)
            writer.writeheader()
    
    return csv_file

def save_to_csv(csv_file, patient_info, lab_values, analysis_results):
    """Save patient data and analysis results to CSV"""
    # Prepare data row
    data_row = {
        'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'patient_id': patient_info.get('patient_id', ''),
        'age': patient_info.get('age', ''),
        'gender': patient_info.get('gender', ''),
    }
    
    # Add lab values (fill missing with empty string)
    all_lab_params = [
        'glucose', 'lactate', 'pyruvate', 'urea', 'creatinine', 'bun', 'ammonia',
        'alt', 'ast', 'bilirubin', 'vitamin_d', 'pth', 'tsh', 'cortisol',
        'insulin', 'c_peptide', 'estrogen', 'progesterone', 'testosterone',
        'apo_a', 'apo_b', 'total_cholesterol', 'ldl', 'hdl', 'triglycerides',
        'phosphorus', 'magnesium', 'uric_acid', 'anti_tpo', 'ana', 'anti_ccp',
        'crp', 'fgf_23', 'homocysteine', 'glutamine', 'glutamate', 'gsh', 'mda'
    ]
    
    for param in all_lab_params:
        data_row[param] = lab_values.get(param, '')
    
    # Add analysis results
    if analysis_results and analysis_results.get('sorted_pathways'):
        data_row.update({
            'top_pathway': analysis_results['sorted_pathways'][0][1]['pathway']['name'],
            'top_score': analysis_results['sorted_pathways'][0][1]['total_score'],
            'directly_affected_count': len(analysis_results.get('directly_affected', [])),
            'at_risk_count': len(analysis_results.get('at_risk', []))
        })
    else:
        data_row.update({
            'top_pathway': '',
            'top_score': '',
            'directly_affected_count': '',
            'at_risk_count': ''
        })
    
    # Write to CSV
    with open(csv_file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=data_row.keys())
        writer.writerow(data_row)

def load_csv_data(csv_file):
    """Load all data from CSV into a DataFrame"""
    if os.path.exists(csv_file):
        df = pd.read_csv(csv_file)
        # Convert timestamp to datetime for better sorting
        if 'timestamp' in df.columns:
            df['timestamp'] = pd.to_datetime(df['timestamp'])
            df = df.sort_values('timestamp', ascending=False)
        return df
    return pd.DataFrame()

# Initialize CSV storage
csv_file = init_csv_storage()

# Remove Streamlit branding and customize appearance
st.set_page_config(
    page_title="AIKEGG Metabolic Pathway Analyzer", 
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Hide Streamlit default elements with CSS
hide_streamlit_style = """
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
.stDeployButton {display:none;}
.stDecoration {display:none;}
[data-testid="stToolbar"] {display: none;}
[data-testid="stHeader"] {display: none;}
[data-testid="stStatusWidget"] {display: none;}
#root > div:nth-child(1) > div > div > div > div > section > div {padding-top: 0rem;}
</style>
"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

@st.cache_data  
def get_metabolic_database():
    return {
        'pathways': {
            'hsa00010': {
                'name': 'Glycolysis / Gluconeogenesis',
                'compounds': ['C00031', 'C00103', 'C00111', 'C00236', 'C00022', 'C00118', 'C00668'],  # Glucose, fructose, pyruvate, ATP
                'compound_names': ['D-Glucose', 'D-Fructose', 'Glycerone phosphate', '3-Phospho-D-glycerate', 'Pyruvate', 'D-Glyceraldehyde 3-phosphate', 'ATP'],
                'enzymes': ['EC:2.7.1.1', 'EC:5.3.1.9', 'EC:2.7.1.11', 'EC:1.2.1.12', 'EC:4.1.2.13', 'EC:2.7.2.3'],
                'enzyme_names': ['Hexokinase', 'Glucose-6-phosphate isomerase', '6-phosphofructokinase', 'Glyceraldehyde-3-phosphate dehydrogenase', 'Aldolase', 'Phosphoglycerate kinase'],
                'reactions': ['R00299', 'R00771', 'R01068', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Glucose phosphorylation', 'Glucose-6-phosphate isomerization', 'Fructose-6-phosphate phosphorylation', 'Pyruvate formation', 'Aldol cleavage', 'ATP generation'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'reversible', 'reversible'],
                'deltaG': -73.3,
                'disease': 'Diabetes Mellitus',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00010',
                'clinical_significance': 'Primary glucose metabolism pathway',
                'affected_organs': ['Liver', 'Muscle', 'Brain', 'Red blood cells'],
                'biomarkers': ['Glucose', 'Lactate', 'Pyruvate', 'ATP/ADP ratio']
            },
            'hsa00020': {
                'name': 'TCA Cycle',
                'compounds': ['C00022', 'C00158', 'C00042', 'C00036', 'C00122', 'C00417', 'C00026'],  # Pyruvate, citrate, oxaloacetate, NADH
                'compound_names': ['Pyruvate', 'Citrate', 'Oxaloacetate', 'Succinate', 'Fumarate', 'cis-Aconitate', 'Oxoglutarate'],
                'enzymes': ['EC:1.2.4.1', 'EC:4.2.1.3', 'EC:1.1.1.37', 'EC:1.3.5.1', 'EC:4.2.1.2', 'EC:1.2.4.2'],
                'enzyme_names': ['Pyruvate dehydrogenase', 'Aconitase', 'Malate dehydrogenase', 'Succinate dehydrogenase', 'Fumarase', 'Oxoglutarate dehydrogenase'],
                'reactions': ['R00209', 'R00351', 'R00709', 'R02164', 'R01082', 'R00243'],
                'reaction_names': ['Pyruvate decarboxylation', 'Citrate isomerization', 'Malate oxidation', 'Succinate oxidation', 'Fumarate hydration', 'Oxoglutarate decarboxylation'],
                'reaction_types': ['irreversible', 'reversible', 'reversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -62.8,
                'disease': 'Mitochondrial Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00020',
                'clinical_significance': 'Central energy production pathway',
                'affected_organs': ['All tissues', 'Mitochondria'],
                'biomarkers': ['Pyruvate', 'Lactate', 'Citrate', 'NADH/NAD+ ratio']
            },
            'hsa00220': {
                'name': 'Urea Cycle',
                'compounds': ['C00062', 'C00169', 'C00077', 'C00327', 'C00152', 'C00049', 'C00025'],  # Urea, ornithine, arginine, glutamine
                'compound_names': ['Urea', 'Ornithine', 'Arginine', 'Citrulline', 'Alanine', 'Aspartate', 'Glutamate'],
                'enzymes': ['EC:6.3.4.5', 'EC:2.1.3.3', 'EC:3.5.3.1', 'EC:2.6.1.13', 'EC:6.3.1.2', 'EC:2.7.3.3'],
                'enzyme_names': ['Carbamoyl phosphate synthetase', 'Ornithine transcarbamylase', 'Arginase', 'Ornithine aminotransferase', 'Glutamine synthetase', 'Carbamate kinase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R00355', 'R00259', 'R01063'],
                'reaction_names': ['Ammonia detoxification', 'Citrulline formation', 'Urea formation', 'Ornithine transamination', 'Glutamine synthesis', 'Carbamoyl phosphate synthesis'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -38.2,
                'disease': 'Urea Cycle Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00220',
                'clinical_significance': 'Nitrogen waste elimination',
                'affected_organs': ['Liver', 'Kidney'],
                'biomarkers': ['Urea', 'Ammonia', 'BUN', 'Citrulline', 'Glutamine']
            },
            'hsa00250': {
                'name': 'Amino Acid Metabolism',
                'compounds': ['C00077', 'C00123', 'C00152', 'C00049', 'C00122', 'C00025', 'C00026'],  # Arginine, leucine, alanine, aspartate
                'compound_names': ['Arginine', 'Leucine', 'Alanine', 'Aspartate', 'Fumarate', 'Glutamate', 'Oxoglutarate'],
                'enzymes': ['EC:2.6.1.1', 'EC:2.6.1.2', 'EC:2.6.1.6', 'EC:4.3.1.1', 'EC:1.4.1.2', 'EC:1.4.1.3'],
                'enzyme_names': ['Alanine transaminase', 'Aspartate transaminase', 'Branched-chain-amino-acid transaminase', 'Argininosuccinate lyase', 'Glutamate dehydrogenase (NAD+)', 'Glutamate dehydrogenase (NADP+)'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01015', 'R00259', 'R00355'],
                'reaction_names': ['Transamination', 'Deamination', 'Decarboxylation', 'Argininosuccinate cleavage', 'Glutamate synthesis', 'Ornithine transamination'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'reversible'],
                'deltaG': -45.7,
                'disease': 'Aminoacidopathies',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00250',
                'clinical_significance': 'Protein metabolism and nitrogen balance',
                'affected_organs': ['Liver', 'Muscle', 'Kidney'],
                'biomarkers': ['Creatinine', 'ALT', 'AST', 'Ammonia', 'Glutamine']
            },
            'hsa00860': {
                'name': 'Porphyrin and Chlorophyll Metabolism',
                'compounds': ['C00032', 'C00037', 'C00049', 'C00251', 'C01024', 'C00091', 'C00500'],  # Heme precursors
                'compound_names': ['Heme', 'Glycine', 'Succinyl-CoA', 'Porphobilinogen', 'Coproporphyrinogen III', 'Uroporphyrinogen III', 'Protoporphyrin IX'],
                'enzymes': ['EC:4.2.1.24', 'EC:4.99.1.1', 'EC:1.3.3.3', 'EC:2.1.1.107', 'EC:4.1.1.37', 'EC:1.3.3.4'],
                'enzyme_names': ['Porphobilinogen deaminase', 'Ferrochelatase', 'Coproporphyrinogen oxidase', 'Uroporphyrinogen III methyltransferase', 'Hydroxymethylbilane synthase', 'Protoporphyrinogen oxidase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Heme biosynthesis', 'Porphyrin formation', 'Bilirubin production', 'Uroporphyrinogen decarboxylation', 'Coproporphyrinogen oxidation', 'Protoporphyrinogen oxidation'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible'],
                'deltaG': -52.3,
                'disease': 'Porphyria',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00860',
                'clinical_significance': 'Heme and bilirubin metabolism',
                'affected_organs': ['Liver', 'Bone marrow', 'Spleen'],
                'biomarkers': ['Bilirubin', 'Porphobilinogen', 'Protoporphyrin']
            },
            'hsa00140': {
                'name': 'Steroid Hormone Biosynthesis',
                'compounds': ['C00410', 'C00735', 'C00533', 'C01953', 'C00280', 'C01176', 'C01227'],  # Steroid hormones
                'compound_names': ['Progesterone', 'Cortisol', 'Testosterone', 'Estradiol', 'Androstenedione', '17-Hydroxyprogesterone', 'Dehydroepiandrosterone'],
                'enzymes': ['EC:1.14.14.1', 'EC:1.14.15.4', 'EC:1.1.1.62', 'EC:1.14.14.19', 'EC:1.1.1.146', 'EC:1.1.1.239'],
                'enzyme_names': ['Aromatase', '21-hydroxylase', '17β-hydroxysteroid dehydrogenase', '17α-hydroxylase', '3β-hydroxysteroid dehydrogenase', '11β-hydroxysteroid dehydrogenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Steroidogenesis', 'Hydroxylation', 'Reduction', 'Side-chain cleavage', 'Aromatization', 'Dehydrogenation'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -68.9,
                'disease': 'Endocrine Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00140',
                'clinical_significance': 'Steroid hormone production',
                'affected_organs': ['Adrenal glands', 'Gonads', 'Placenta'],
                'biomarkers': ['Cortisol', 'Estrogen', 'Progesterone', 'Testosterone', 'Vitamin D', 'DHEA-S']
            },
            'hsa04910': {
                'name': 'Insulin Signaling Pathway',
                'compounds': ['C00031', 'C00267', 'C00089', 'C00103', 'C00668', 'C00008', 'C00009'],  # Glucose, insulin-related
                'compound_names': ['D-Glucose', 'α-D-Glucose', 'Sucrose', 'D-Fructose', 'ATP', 'ADP', 'Phosphate'],
                'enzymes': ['EC:2.7.1.1', 'EC:2.7.11.1', 'EC:3.1.3.48', 'EC:2.7.1.153', 'EC:2.7.11.2', 'EC:3.1.3.16'],
                'enzyme_names': ['Hexokinase', 'Protein kinase B', 'Protein-tyrosine-phosphatase', 'Phosphoinositide 3-kinase', 'Protein kinase C', 'Phosphorylase phosphatase'],
                'reactions': ['R00299', 'R01786', 'R02736', 'R02740', 'R02741', 'R02742'],
                'reaction_names': ['Glucose phosphorylation', 'Insulin receptor activation', 'Signal transduction', 'PI3K activation', 'GLUT4 translocation', 'Glycogen synthase activation'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -85.3,
                'disease': 'Diabetes Mellitus',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04910',
                'clinical_significance': 'Glucose homeostasis and insulin sensitivity',
                'affected_organs': ['Pancreas', 'Liver', 'Muscle', 'Adipose tissue'],
                'biomarkers': ['Glucose', 'Insulin', 'C-peptide', 'HbA1c']
            },
            'hsa04930': {
                'name': 'Type II Diabetes Mellitus',
                'compounds': ['C00031', 'C00159', 'C00267', 'C00369', 'C00668', 'C00008', 'C00009'],  # Glucose, mannose, starch
                'compound_names': ['D-Glucose', 'D-Mannose', 'α-D-Glucose', 'Starch', 'ATP', 'ADP', 'Phosphate'],
                'enzymes': ['EC:2.7.1.1', 'EC:3.1.3.48', 'EC:3.1.4.4', 'EC:2.7.11.1', 'EC:2.7.1.153', 'EC:3.1.3.16'],
                'enzyme_names': ['Hexokinase', 'Protein-tyrosine-phosphatase', 'Phospholipase C', 'Protein kinase B', 'Phosphoinositide 3-kinase', 'Phosphorylase phosphatase'],
                'reactions': ['R00299', 'R02736', 'R02740', 'R02741', 'R02742', 'R02743'],
                'reaction_names': ['Glucose phosphorylation', 'Insulin resistance', 'Signal transduction', 'PI3K activation', 'GLUT4 translocation', 'Glycogen synthase inhibition'],
                'reaction_types': ['irreversible', 'irreversible', 'reversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -68.7,
                'disease': 'Diabetes Mellitus',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04930',
                'clinical_significance': 'Insulin resistance and metabolic dysfunction',
                'affected_organs': ['Pancreas', 'Liver', 'Muscle', 'Adipose tissue'],
                'biomarkers': ['Glucose', 'Insulin', 'HbA1c', 'Leptin']
            },
            'hsa04920': {
                'name': 'Adipocytokine Signaling Pathway',
                'compounds': ['C00031', 'C00162', 'C00422', 'C02530', 'C00668', 'C00008', 'C00009'],  # Glucose, fatty acids, adipokines
                'compound_names': ['D-Glucose', 'Fatty acid', 'Triacylglycerol', 'Adiponectin', 'ATP', 'ADP', 'Phosphate'],
                'enzymes': ['EC:2.7.1.1', 'EC:2.7.11.1', 'EC:3.1.1.3', 'EC:2.7.11.27', 'EC:2.7.11.2', 'EC:3.1.3.16'],
                'enzyme_names': ['Hexokinase', 'AMP-activated protein kinase', 'Hormone-sensitive lipase', 'JAK-STAT kinase', 'Protein kinase C', 'Phosphorylase phosphatase'],
                'reactions': ['R00299', 'R01786', 'R02736', 'R02740', 'R02741', 'R02742'],
                'reaction_names': ['Lipolysis', 'Fatty acid oxidation', 'Insulin sensitization', 'Adiponectin signaling', 'Leptin signaling', 'Resistin signaling'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -72.5,
                'disease': 'Metabolic Syndrome',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04920',
                'clinical_significance': 'Adipose tissue signaling and metabolism',
                'affected_organs': ['Adipose tissue', 'Liver', 'Muscle'],
                'biomarkers': ['Leptin', 'Adiponectin', 'Resistin', 'Free fatty acids']
            },
            'hsa04922': {
                'name': 'Glucagon Signaling Pathway',
                'compounds': ['C00031', 'C00159', 'C00267', 'C00369', 'C00668', 'C00008', 'C00009'],  # Glucose, glycogen
                'compound_names': ['D-Glucose', 'D-Mannose', 'α-D-Glucose', 'Starch', 'ATP', 'ADP', 'Phosphate'],
                'enzymes': ['EC:2.7.1.1', 'EC:2.7.11.1', 'EC:3.1.3.16', 'EC:2.7.11.11', 'EC:2.7.11.19', 'EC:3.1.3.15'],
                'enzyme_names': ['Hexokinase', 'Protein kinase A', 'Phosphorylase phosphatase', 'Glycogen synthase kinase', 'Calcium/calmodulin-dependent protein kinase', 'Phosphoprotein phosphatase'],
                'reactions': ['R00299', 'R01786', 'R02736', 'R02740', 'R02741', 'R02742'],
                'reaction_names': ['Glycogenolysis', 'Gluconeogenesis', 'Lipolysis', 'PKA activation', 'Calcium signaling', 'Glycogen synthase inhibition'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -58.9,
                'disease': 'Diabetes Mellitus',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04922',
                'clinical_significance': 'Counter-regulatory hormone signaling',
                'affected_organs': ['Liver', 'Adipose tissue'],
                'biomarkers': ['Glucose', 'Glucagon', 'cAMP']
            },
            'hsa04950': {
                'name': 'Maturity Onset Diabetes of the Young (MODY)',
                'compounds': ['C00031', 'C00103', 'C00267', 'C00668', 'C00008', 'C00009', 'C00020'],  # Glucose, ATP
                'compound_names': ['D-Glucose', 'D-Fructose', 'α-D-Glucose', 'ATP', 'ADP', 'Phosphate', 'AMP'],
                'enzymes': ['EC:2.7.1.1', 'EC:2.7.1.2', 'EC:3.1.3.48', 'EC:2.7.11.1', 'EC:2.7.11.27', 'EC:3.1.3.16'],
                'enzyme_names': ['Hexokinase', 'Glucokinase', 'Protein-tyrosine-phosphatase', 'Protein kinase B', 'JAK-STAT kinase', 'Phosphorylase phosphatase'],
                'reactions': ['R00299', 'R01786', 'R02736', 'R02740', 'R02741', 'R02742'],
                'reaction_names': ['Glucose sensing', 'Insulin secretion', 'Signal transduction', 'PI3K activation', 'JAK-STAT signaling', 'Glycogen synthase activation'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -78.5,
                'disease': 'Diabetes Mellitus',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04950',
                'clinical_significance': 'Monogenic diabetes with early onset',
                'affected_organs': ['Pancreatic beta cells', 'Liver'],
                'biomarkers': ['Glucose', 'C-peptide', 'Proinsulin']
            },
            'hsa04979': {
                'name': 'Cholesterol Metabolism',
                'compounds': ['C00187', 'C00300', 'C00410', 'C01154', 'C00270', 'C01100', 'C05381'],  # Cholesterol, bile acids
                'compound_names': ['Cholesterol', 'Creatine', 'Progesterone', 'Cholesteryl ester', 'Mevalonate', '7-Dehydrocholesterol', 'Bile acid'],
                'enzymes': ['EC:1.1.1.34', 'EC:1.14.13.70', 'EC:2.3.1.26', 'EC:1.14.14.25', 'EC:1.3.1.21', 'EC:1.14.18.9'],
                'enzyme_names': ['HMG-CoA reductase', 'Cholesterol 7α-hydroxylase', 'Cholesterol acyltransferase', 'Lanosterol 14α-demethylase', '7-Dehydrocholesterol reductase', 'Squalene monooxygenase'],
                'reactions': ['R00351', 'R01878', 'R02540', 'R02541', 'R02542', 'R02543'],
                'reaction_names': ['Cholesterol biosynthesis', 'Bile acid synthesis', 'Cholesterol esterification', 'Lanosterol demethylation', 'Squalene epoxidation', 'Bile acid conjugation'],
                'reaction_types': ['irreversible', 'reversible', 'reversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -142.5,
                'disease': 'Dyslipidemia',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04979',
                'clinical_significance': 'Central to lipid homeostasis and cardiovascular risk',
                'affected_organs': ['Liver', 'Intestine', 'Adrenal glands'],
                'biomarkers': ['Total cholesterol', 'LDL', 'HDL', 'Apo B', 'Triglycerides']
            },
            'hsa03320': {
                'name': 'PPAR Signaling Pathway',
                'compounds': ['C00187', 'C00410', 'C02530', 'C06427', 'C00162', 'C00219', 'C00249'],  # Lipids, fatty acids
                'compound_names': ['Cholesterol', 'Progesterone', 'Fatty acid', 'Prostaglandin', 'Arachidonic acid', 'Eicosapentaenoic acid', 'Docosahexaenoic acid'],
                'enzymes': ['EC:1.13.11.52', 'EC:1.14.19.3', 'EC:1.3.1.22', 'EC:1.14.19.2', 'EC:1.13.11.34', 'EC:1.13.11.33'],
                'enzyme_names': ['Lipoxygenase', 'Fatty acid desaturase', 'Acyl-CoA dehydrogenase', 'Prostaglandin-endoperoxide synthase', 'Cyclooxygenase', 'Lipoxygenase'],
                'reactions': ['R01878', 'R02540', 'R04779', 'R04780', 'R04781', 'R04782'],
                'reaction_names': ['Lipid oxidation', 'Fatty acid metabolism', 'Beta-oxidation', 'Prostaglandin synthesis', 'Leukotriene synthesis', 'Lipoxin synthesis'],
                'reaction_types': ['reversible', 'reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -118.9,
                'disease': 'Dyslipidemia',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa03320',
                'clinical_significance': 'Regulates lipid metabolism and inflammation',
                'affected_organs': ['Liver', 'Adipose tissue', 'Muscle'],
                'biomarkers': ['Triglycerides', 'Free fatty acids', 'Apo A1']
            },
            'hsa00561': {
                'name': 'Glycerolipid Metabolism',
                'compounds': ['C00044', 'C00116', 'C00162', 'C00422', 'C00157', 'C00160', 'C00681'],  # Glycerol, triglycerides
                'compound_names': ['GTP', 'Glycerol', 'Fatty acid', 'Triacylglycerol', 'Phosphatidate', 'Diacylglycerol', 'CDP-diacylglycerol'],
                'enzymes': ['EC:2.3.1.15', 'EC:3.1.1.3', 'EC:3.1.1.4', 'EC:2.3.1.51', 'EC:2.7.1.107', 'EC:3.1.3.4'],
                'enzyme_names': ['Glycerol-3-phosphate acyltransferase', 'Triacylglycerol lipase', 'Phospholipase A2', '1-acylglycerol-3-phosphate O-acyltransferase', 'Diacylglycerol kinase', 'Phosphatidate phosphatase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Glycerol phosphorylation', 'Triacylglycerol synthesis', 'Lipolysis', 'Phosphatidic acid synthesis', 'Diacylglycerol synthesis', 'CDP-diacylglycerol synthesis'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -92.4,
                'disease': 'Dyslipidemia',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00561',
                'clinical_significance': 'Triglyceride synthesis and breakdown',
                'affected_organs': ['Liver', 'Adipose tissue', 'Muscle'],
                'biomarkers': ['Triglycerides', 'Glycerol', 'Free fatty acids']
            },
            'hsa03320': {
                'name': 'PPAR Signaling Pathway',
                'compounds': ['C00187', 'C00410', 'C02530', 'C06427', 'C00162', 'C00219', 'C00249'],  # Lipids, fatty acids
                'compound_names': ['Cholesterol', 'Progesterone', 'Fatty acid', 'Prostaglandin', 'Arachidonic acid', 'Eicosapentaenoic acid', 'Docosahexaenoic acid'],
                'enzymes': ['EC:1.13.11.52', 'EC:1.14.19.3', 'EC:1.3.1.22', 'EC:1.14.19.2', 'EC:1.13.11.34', 'EC:1.13.11.33'],
                'enzyme_names': ['Lipoxygenase', 'Fatty acid desaturase', 'Acyl-CoA dehydrogenase', 'Prostaglandin-endoperoxide synthase', 'Cyclooxygenase', 'Lipoxygenase'],
                'reactions': ['R01878', 'R02540', 'R04779', 'R04780', 'R04781', 'R04782'],
                'reaction_names': ['Lipid oxidation', 'Fatty acid metabolism', 'Beta-oxidation', 'Prostaglandin synthesis', 'Leukotriene synthesis', 'Lipoxin synthesis'],
                'reaction_types': ['reversible', 'reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -118.9,
                'disease': 'Dyslipidemia',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa03320',
                'clinical_significance': 'Regulates lipid metabolism and inflammation',
                'affected_organs': ['Liver', 'Adipose tissue', 'Muscle'],
                'biomarkers': ['Triglycerides', 'Free fatty acids', 'Apo A1']
            },
            'hsa00561': {
                'name': 'Glycerolipid Metabolism',
                'compounds': ['C00044', 'C00116', 'C00162', 'C00422', 'C00157', 'C00160', 'C00681'],  # Glycerol, triglycerides
                'compound_names': ['GTP', 'Glycerol', 'Fatty acid', 'Triacylglycerol', 'Phosphatidate', 'Diacylglycerol', 'CDP-diacylglycerol'],
                'enzymes': ['EC:2.3.1.15', 'EC:3.1.1.3', 'EC:3.1.1.4', 'EC:2.3.1.51', 'EC:2.7.1.107', 'EC:3.1.3.4'],
                'enzyme_names': ['Glycerol-3-phosphate acyltransferase', 'Triacylglycerol lipase', 'Phospholipase A2', '1-acylglycerol-3-phosphate O-acyltransferase', 'Diacylglycerol kinase', 'Phosphatidate phosphatase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Glycerol phosphorylation', 'Triacylglycerol synthesis', 'Lipolysis', 'Phosphatidic acid synthesis', 'Diacylglycerol synthesis', 'CDP-diacylglycerol synthesis'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -92.4,
                'disease': 'Dyslipidemia',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00561',
                'clinical_significance': 'Triglyceride synthesis and breakdown',
                'affected_organs': ['Liver', 'Adipose tissue', 'Muscle'],
                'biomarkers': ['Triglycerides', 'Glycerol', 'Free fatty acids']
            },
            'hsa04960': {
                'name': 'Aldosterone-Regulated Sodium Reabsorption',
                'compounds': ['C00009', 'C00080', 'C01330', 'C00788', 'C00238', 'C00076', 'C00288'],  # Electrolytes
                'compound_names': ['Orthophosphate', 'H+', 'Sodium', 'L-Histidine', 'Potassium', 'Calcium', 'Bicarbonate'],
                'enzymes': ['EC:1.6.3.1', 'EC:3.6.3.9', 'EC:4.2.1.1', 'EC:3.6.3.8', 'EC:1.14.15.4', 'EC:1.1.1.146'],
                'enzyme_names': ['NAD(P)H oxidase', 'Na+/K+-ATPase', 'Carbonic anhydrase', 'Ca2+-ATPase', 'Aldosterone synthase', '11β-hydroxysteroid dehydrogenase'],
                'reactions': ['R00114', 'R00115', 'R00116', 'R00117', 'R00118', 'R00119'],
                'reaction_names': ['Sodium transport', 'Potassium exchange', 'pH regulation', 'Calcium transport', 'Aldosterone synthesis', 'Cortisol metabolism'],
                'reaction_types': ['irreversible', 'reversible', 'reversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -45.2,
                'disease': 'Chronic Kidney Disease',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04960',
                'clinical_significance': 'Electrolyte balance and blood pressure regulation',
                'affected_organs': ['Kidney', 'Adrenal cortex'],
                'biomarkers': ['Sodium', 'Potassium', 'Aldosterone', 'Renin']
            },
            'hsa04961': {
                'name': 'Endocrine and Other Factor-Regulated Calcium Reabsorption',
                'compounds': ['C00076', 'C01353', 'C00005', 'C00009', 'C00238', 'C00288', 'C00080'],  # Calcium, phosphate
                'compound_names': ['Calcium', 'Calcitriol', 'NADPH', 'Orthophosphate', 'Potassium', 'Bicarbonate', 'H+'],
                'enzymes': ['EC:3.6.3.8', 'EC:1.1.1.37', 'EC:4.2.1.1', 'EC:1.14.15.18', 'EC:1.14.14.16', 'EC:1.17.7.1'],
                'enzyme_names': ['Ca2+-ATPase', 'Malate dehydrogenase', 'Carbonic anhydrase', '25-hydroxyvitamin D3 1α-hydroxylase', 'Vitamin D 25-hydroxylase', 'Vitamin K epoxide reductase'],
                'reactions': ['R00114', 'R00115', 'R00116', 'R00117', 'R00118', 'R00119'],
                'reaction_names': ['Calcium transport', 'Energy metabolism', 'Acid-base balance', 'Vitamin D activation', 'Vitamin D hydroxylation', 'Vitamin K cycle'],
                'reaction_types': ['irreversible', 'reversible', 'reversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -38.7,
                'disease': 'Chronic Kidney Disease',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04961',
                'clinical_significance': 'Calcium homeostasis and bone metabolism',
                'affected_organs': ['Kidney', 'Parathyroid', 'Bone'],
                'biomarkers': ['Calcium', 'Phosphate', 'PTH', 'Vitamin D']
            },
            'hsa04966': {
                'name': 'Collecting Duct Acid Secretion',
                'compounds': ['C00014', 'C00058', 'C00238', 'C01353', 'C00080', 'C00288', 'C00076'],  # Ammonia, potassium
                'compound_names': ['Ammonia', 'Formate', 'Potassium', 'Calcitriol', 'H+', 'Bicarbonate', 'Calcium'],
                'enzymes': ['EC:3.6.3.1', 'EC:4.2.1.1', 'EC:1.1.1.37', 'EC:1.4.1.2', 'EC:1.4.1.3', 'EC:1.4.1.4'],
                'enzyme_names': ['H+-ATPase', 'Carbonic anhydrase', 'Malate dehydrogenase', 'Glutamate dehydrogenase (NAD+)', 'Glutamate dehydrogenase (NADP+)', 'Glutamate dehydrogenase (FAD+)'],
                'reactions': ['R00114', 'R00115', 'R00116', 'R00117', 'R00118', 'R00119'],
                'reaction_names': ['Proton secretion', 'Bicarbonate formation', 'Acid-base regulation', 'Glutamate deamination', 'Ammonia production', 'Glutamine synthesis'],
                'reaction_types': ['irreversible', 'reversible', 'reversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -52.1,
                'disease': 'Chronic Kidney Disease',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04966',
                'clinical_significance': 'Renal acid-base balance and pH homeostasis',
                'affected_organs': ['Kidney collecting duct'],
                'biomarkers': ['Ammonia', 'Bicarbonate', 'Urine pH']
            },
            'hsa00240': {
                'name': 'Pyrimidine Metabolism',
                'compounds': ['C00106', 'C00262', 'C00178', 'C00286', 'C00299', 'C00311', 'C00385'],  # Uracil, uridine
                'compound_names': ['Uracil', 'Uridine', 'Thymine', 'Deoxyuridine', 'Cytidine', 'Deoxycytidine', 'Orotate'],
                'enzymes': ['EC:2.4.2.3', 'EC:3.5.4.1', 'EC:1.3.1.2', 'EC:2.1.1.45', 'EC:2.7.4.9', 'EC:3.5.4.12'],
                'enzyme_names': ['Uridine phosphorylase', 'Cytidine deaminase', 'Dihydroorotate dehydrogenase', 'Thymidylate synthase', 'Nucleoside-diphosphate kinase', 'Deoxycytidylate deaminase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Pyrimidine salvage', 'Deamination', 'Reduction', 'Methylation', 'Phosphorylation', 'Deamination'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -48.3,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00240',
                'clinical_significance': 'Nucleotide metabolism',
                'affected_organs': ['Liver', 'Bone marrow'],
                'biomarkers': ['Uric acid', 'Orotic acid']
            },
            'hsa00230': {
                'name': 'Purine Metabolism',
                'compounds': ['C00262', 'C00366', 'C00147', 'C00242', 'C00130', 'C00212', 'C00362'],  # Purines
                'compound_names': ['Uridine', 'Hypoxanthine', 'Adenine', 'Guanine', 'IMP', 'XMP', 'GMP'],
                'enzymes': ['EC:2.4.2.1', 'EC:3.5.4.4', 'EC:1.17.3.2', 'EC:2.7.4.3', 'EC:6.3.4.4', 'EC:2.1.1.45'],
                'enzyme_names': ['Purine nucleoside phosphorylase', 'Adenosine deaminase', 'Xanthine oxidase', 'Adenylate kinase', 'FGAM synthetase', 'Guanosine monophosphate synthase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Purine salvage', 'Deamination', 'Oxidation', 'Phosphorylation', 'Synthesis', 'Methylation'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'reversible', 'irreversible', 'irreversible'],
                'deltaG': -55.7,
                'disease': 'Gout',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00230',
                'clinical_significance': 'Purine nucleotide metabolism',
                'affected_organs': ['Liver', 'Kidney'],
                'biomarkers': ['Uric acid', 'Hypoxanthine', 'Xanthine']
            },
            'hsa00480': {
                'name': 'Glutathione Metabolism',
                'compounds': ['C00051', 'C00025', 'C00127', 'C00669', 'C00167', 'C00097', 'C00094'],  # Glutathione
                'compound_names': ['Glutathione', 'Glutamate', 'Cysteine', 'Oxidized glutathione', 'Ascorbate', 'Cystine', 'Sulfite'],
                'enzymes': ['EC:1.11.1.9', 'EC:1.8.1.7', 'EC:2.5.1.18', 'EC:1.8.5.1', 'EC:1.8.4.2', 'EC:1.8.3.2'],
                'enzyme_names': ['Glutathione peroxidase', 'Glutathione reductase', 'Glutathione transferase', 'Glutathione dehydrogenase', 'Glutathione-disulfide reductase', 'Sulfite oxidase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Glutathione oxidation', 'Reduction', 'Conjugation', 'Detoxification', 'Disulfide reduction', 'Sulfite oxidation'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -42.8,
                'disease': 'Oxidative Stress',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00480',
                'clinical_significance': 'Antioxidant defense system',
                'affected_organs': ['Liver', 'Red blood cells'],
                'biomarkers': ['GSH', 'MDA', 'Oxidized glutathione']
            },
            'hsa04964': {
                'name': 'Proximal Tubule Bicarbonate Reclamation',
                'compounds': ['C00014', 'C00288', 'C00076', 'C01353', 'C00080', 'C00238', 'C00025'],  # Bicarbonate
                'compound_names': ['Ammonia', 'HCO3-', 'Calcium', 'Calcitriol', 'H+', 'Potassium', 'Glutamate'],
                'enzymes': ['EC:4.2.1.1', 'EC:3.6.3.10', 'EC:1.1.1.37', 'EC:1.4.1.2', 'EC:1.4.1.3', 'EC:1.4.1.4'],
                'enzyme_names': ['Carbonic anhydrase', 'Na+/HCO3- cotransporter', 'Malate dehydrogenase', 'Glutamate dehydrogenase (NAD+)', 'Glutamate dehydrogenase (NADP+)', 'Glutamate dehydrogenase (FAD+)'],
                'reactions': ['R00114', 'R00115', 'R00116', 'R00117', 'R00118', 'R00119'],
                'reaction_names': ['Bicarbonate reabsorption', 'Acid secretion', 'Ammoniagenesis', 'Glutamate deamination', 'Ammonia production', 'Glutamine synthesis'],
                'reaction_types': ['irreversible', 'reversible', 'reversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -35.9,
                'disease': 'Chronic Kidney Disease',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04964',
                'clinical_significance': 'Renal acid-base handling',
                'affected_organs': ['Kidney proximal tubule'],
                'biomarkers': ['Bicarbonate', 'Ammonia', 'Urine pH']
            },
            'hsa04970': {
                'name': 'Salivary Secretion',
                'compounds': ['C00076', 'C01330', 'C00238', 'C00080', 'C00288', 'C00025', 'C00014'],  # Electrolytes
                'compound_names': ['Calcium', 'Sodium', 'Potassium', 'H+', 'Bicarbonate', 'Glutamate', 'Ammonia'],
                'enzymes': ['EC:3.6.3.8', 'EC:3.6.3.9', 'EC:4.2.1.1', 'EC:3.4.21.1', 'EC:3.2.1.1', 'EC:3.2.1.20'],
                'enzyme_names': ['Ca2+-ATPase', 'Na+/K+-ATPase', 'Carbonic anhydrase', 'Trypsin', 'α-Amylase', 'β-Galactosidase'],
                'reactions': ['R00114', 'R00115', 'R00116', 'R00117', 'R00118', 'R00119'],
                'reaction_names': ['Ion transport', 'Fluid secretion', 'pH regulation', 'Protein digestion', 'Carbohydrate digestion', 'Mucin breakdown'],
                'reaction_types': ['irreversible', 'reversible', 'reversible', 'irreversible', 'irreversible', 'irreversible'],
                'deltaG': -28.7,
                'disease': 'Sjögren Syndrome',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04970',
                'clinical_significance': 'Exocrine gland function',
                'affected_organs': ['Salivary glands', 'Lacrimal glands'],
                'biomarkers': ['ANA', 'Anti-SSA/Ro', 'Anti-SSB/La']
            },
            'hsa04630': {
                'name': 'JAK-STAT Signaling Pathway',
                'compounds': ['C00031', 'C00162', 'C00422', 'C02530', 'C00008', 'C00009', 'C00020'],  # Signaling molecules
                'compound_names': ['D-Glucose', 'Fatty acid', 'Triacylglycerol', 'Cytokine', 'ADP', 'Phosphate', 'AMP'],
                'enzymes': ['EC:2.7.10.2', 'EC:2.7.11.1', 'EC:3.1.3.48', 'EC:2.7.11.27', 'EC:2.7.11.2', 'EC:3.1.3.16'],
                'enzyme_names': ['Janus kinase', 'Signal transducer and activator of transcription', 'Protein-tyrosine-phosphatase', 'JAK-STAT kinase', 'Protein kinase C', 'Phosphorylase phosphatase'],
                'reactions': ['R00299', 'R01786', 'R02736', 'R02740', 'R02741', 'R02742'],
                'reaction_names': ['Signal transduction', 'Gene regulation', 'Feedback inhibition', 'Cytokine signaling', 'Inflammatory response', 'Immune regulation'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -65.3,
                'disease': 'Autoimmune Diseases',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04630',
                'clinical_significance': 'Cytokine signaling and immune regulation',
                'affected_organs': ['Immune cells', 'Liver'],
                'biomarkers': ['CRP', 'Anti-CCP', 'IL-6']
            },
            'hsa04660': {
                'name': 'T Cell Receptor Signaling Pathway',
                'compounds': ['C00031', 'C00162', 'C00422', 'C02530', 'C00008', 'C00009', 'C00020'],  # Signaling molecules
                'compound_names': ['D-Glucose', 'Fatty acid', 'Triacylglycerol', 'Antigen', 'ADP', 'Phosphate', 'AMP'],
                'enzymes': ['EC:2.7.10.2', 'EC:2.7.11.1', 'EC:3.1.3.48', 'EC:2.7.11.7', 'EC:2.7.11.10', 'EC:3.1.3.16'],
                'enzyme_names': ['Lck', 'ZAP-70', 'Protein-tyrosine-phosphatase', 'Protein kinase C theta', 'IKK-beta', 'Phosphorylase phosphatase'],
                'reactions': ['R00299', 'R01786', 'R02736', 'R02740', 'R02741', 'R02742'],
                'reaction_names': ['Signal transduction', 'Gene regulation', 'Feedback inhibition', 'NF-κB activation', 'Inflammatory response', 'Immune regulation'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -58.2,
                'disease': 'Autoimmune Diseases',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04660',
                'clinical_significance': 'T cell activation and immune response',
                'affected_organs': ['Thymus', 'Lymph nodes'],
                'biomarkers': ['Anti-CCP', 'ANA', 'TNF-α']
            },
            'hsa04060': {
                'name': 'Cytokine-cytokine Receptor Interaction',
                'compounds': ['C00031', 'C00162', 'C00422', 'C02530', 'C00008', 'C00009', 'C00020'],  # Cytokines
                'compound_names': ['D-Glucose', 'Fatty acid', 'Triacylglycerol', 'Cytokine', 'ADP', 'Phosphate', 'AMP'],
                'enzymes': ['EC:2.7.10.2', 'EC:2.7.11.1', 'EC:3.1.3.48', 'EC:2.7.11.27', 'EC:2.7.11.2', 'EC:3.1.3.16'],
                'enzyme_names': ['Receptor tyrosine kinase', 'JAK', 'Protein-tyrosine-phosphatase', 'JAK-STAT kinase', 'Protein kinase C', 'Phosphorylase phosphatase'],
                'reactions': ['R00299', 'R01786', 'R02736', 'R02740', 'R02741', 'R02742'],
                'reaction_names': ['Signal transduction', 'Gene regulation', 'Feedback inhibition', 'Cytokine signaling', 'Inflammatory response', 'Immune regulation'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -62.4,
                'disease': 'Inflammatory Diseases',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa04060',
                'clinical_significance': 'Cytokine signaling and inflammation',
                'affected_organs': ['Immune cells', 'Various tissues'],
                'biomarkers': ['CRP', 'Anti-TPO', 'IL-1β']
            },
            'hsa00590': {
                'name': 'Arachidonic Acid Metabolism',
                'compounds': ['C00219', 'C00208', 'C00696', 'C00251', 'C00627', 'C00639', 'C00640'],  # Arachidonic acid
                'compound_names': ['Arachidonic acid', 'Prostaglandin H2', 'Leukotriene B4', 'Porphobilinogen', 'Thromboxane A2', 'Prostacyclin', 'Prostaglandin E2'],
                'enzymes': ['EC:1.14.14.1', 'EC:1.13.11.34', 'EC:1.13.11.31', 'EC:5.3.99.1', 'EC:1.1.1.189', 'EC:1.1.1.141'],
                'enzyme_names': ['Cyclooxygenase', 'Lipoxygenase', 'Prostaglandin-endoperoxide synthase', 'Prostaglandin E synthase', 'Prostaglandin dehydrogenase', '15-hydroxyprostaglandin dehydrogenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Arachidonic acid release', 'Prostaglandin synthesis', 'Leukotriene synthesis', 'Thromboxane synthesis', 'Prostacyclin synthesis', 'Lipoxin synthesis'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible'],
                'deltaG': -78.5,
                'disease': 'Inflammatory Diseases',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00590',
                'clinical_significance': 'Inflammatory mediator production',
                'affected_organs': ['Immune cells', 'Platelets'],
                'biomarkers': ['Prostaglandins', 'Leukotrienes', 'CRP']
            },
            'hsa00620': {
                'name': 'Pyruvate Metabolism',
                'compounds': ['C00022', 'C00122', 'C00026', 'C00036', 'C00042', 'C00024', 'C00149'],  # Pyruvate
                'compound_names': ['Pyruvate', 'Fumarate', 'Oxoglutarate', 'Succinate', 'Oxaloacetate', 'Acetyl-CoA', 'Malate'],
                'enzymes': ['EC:1.2.4.1', 'EC:1.1.1.37', 'EC:2.3.1.12', 'EC:4.1.1.31', 'EC:1.1.1.40', 'EC:2.7.9.1'],
                'enzyme_names': ['Pyruvate dehydrogenase', 'Malate dehydrogenase', 'Pyruvate carboxylase', 'Phosphoenolpyruvate carboxykinase', 'Lactate dehydrogenase', 'Pyruvate kinase'],
                'reactions': ['R00209', 'R00709', 'R01015', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Pyruvate decarboxylation', 'Malate oxidation', 'Oxaloacetate formation', 'Phosphoenolpyruvate formation', 'Lactate formation', 'Acetyl-CoA formation'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -52.8,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00620',
                'clinical_significance': 'Central metabolic junction',
                'affected_organs': ['Liver', 'Muscle', 'Brain'],
                'biomarkers': ['Pyruvate', 'Lactate', 'Alanine']
            },
            'hsa00640': {
                'name': 'Propanoate Metabolism',
                'compounds': ['C00163', 'C00122', 'C00042', 'C00149', 'C00026', 'C00036', 'C00022'],  # Propanoate
                'compound_names': ['Propanoate', 'Fumarate', 'Oxaloacetate', 'Malate', 'Oxoglutarate', 'Succinate', 'Pyruvate'],
                'enzymes': ['EC:1.2.1.3', 'EC:1.1.1.35', 'EC:1.1.1.37', 'EC:4.2.1.2', 'EC:1.3.5.1', 'EC:1.2.4.1'],
                'enzyme_names': ['Propionaldehyde dehydrogenase', 'Methylmalonate-semialdehyde dehydrogenase', 'Malate dehydrogenase', 'Fumarase', 'Succinate dehydrogenase', 'Pyruvate dehydrogenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Propanoate activation', 'Methylmalonyl-CoA formation', 'Succinyl-CoA formation', 'TCA cycle integration', 'Energy production', 'Odd-chain fatty acid metabolism'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -45.7,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00640',
                'clinical_significance': 'Odd-chain fatty acid and amino acid metabolism',
                'affected_organs': ['Liver', 'Muscle'],
                'biomarkers': ['Propionate', 'Methylmalonate']
            },
            'hsa00650': {
                'name': 'Butanoate Metabolism',
                'compounds': ['C00246', 'C00162', 'C00036', 'C00122', 'C00042', 'C00026', 'C00149'],  # Butanoate
                'compound_names': ['Butanoate', 'Fatty acid', 'Succinate', 'Fumarate', 'Oxaloacetate', 'Oxoglutarate', 'Malate'],
                'enzymes': ['EC:1.3.8.1', 'EC:1.3.8.7', 'EC:1.3.99.2', 'EC:4.2.1.2', 'EC:1.1.1.37', 'EC:1.2.4.1'],
                'enzyme_names': ['Butyryl-CoA dehydrogenase', 'Short-chain acyl-CoA dehydrogenase', 'Glutaryl-CoA dehydrogenase', 'Fumarase', 'Malate dehydrogenase', 'Pyruvate dehydrogenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Butanoate activation', 'β-oxidation', 'Acetyl-CoA formation', 'TCA cycle integration', 'Energy production', 'Ketone body formation'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -68.3,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00650',
                'clinical_significance': 'Short-chain fatty acid metabolism',
                'affected_organs': ['Liver', 'Colon'],
                'biomarkers': ['Butyrate', 'Acetoacetate']
            },
            'hsa00071': {
                'name': 'Fatty Acid Degradation',
                'compounds': ['C00162', 'C00036', 'C00122', 'C00042', 'C00026', 'C00149', 'C00022'],  # Fatty acids
                'compound_names': ['Fatty acid', 'Succinate', 'Fumarate', 'Oxaloacetate', 'Oxoglutarate', 'Malate', 'Pyruvate'],
                'enzymes': ['EC:1.3.8.1', 'EC:1.3.8.7', 'EC:1.3.99.2', 'EC:4.2.1.2', 'EC:1.1.1.37', 'EC:1.2.4.1'],
                'enzyme_names': ['Acyl-CoA dehydrogenase', 'Enoyl-CoA hydratase', '3-hydroxyacyl-CoA dehydrogenase', 'Fumarase', 'Malate dehydrogenase', 'Pyruvate dehydrogenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Fatty acid activation', 'β-oxidation', 'Acetyl-CoA formation', 'TCA cycle integration', 'Energy production', 'Ketone body formation'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -85.2,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00071',
                'clinical_significance': 'Fatty acid breakdown and energy production',
                'affected_organs': ['Liver', 'Muscle', 'Adipose tissue'],
                'biomarkers': ['Free fatty acids', 'Ketones', 'Carnitine']
            },
            'hsa00061': {
                'name': 'Fatty Acid Biosynthesis',
                'compounds': ['C00024', 'C00083', 'C00162', 'C00246', 'C00422', 'C00641', 'C00681'],  # Fatty acids
                'compound_names': ['Acetyl-CoA', 'Malonyl-CoA', 'Fatty acid', 'Butanoate', 'Triacylglycerol', 'Palmitate', 'Stearate'],
                'enzymes': ['EC:2.3.1.85', 'EC:1.1.1.100', 'EC:1.3.1.9', 'EC:2.3.1.16', 'EC:2.3.1.20', 'EC:2.3.1.75'],
                'enzyme_names': ['Acetyl-CoA carboxylase', 'Fatty acid synthase', 'Enoyl-ACP reductase', '3-oxoacyl-ACP synthase', '3-oxoacyl-ACP reductase', 'Malonyl-CoA-ACP transacylase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Acetyl-CoA carboxylation', 'Malonyl-CoA formation', 'Fatty acid elongation', 'Desaturation', 'Esterification', 'Triacylglycerol synthesis'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -92.7,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00061',
                'clinical_significance': 'Fatty acid synthesis and storage',
                'affected_organs': ['Liver', 'Adipose tissue'],
                'biomarkers': ['Free fatty acids', 'Triglycerides', 'VLDL']
            },
            'hsa00630': {
                'name': 'Glyoxylate and Dicarboxylate Metabolism',
                'compounds': ['C00048', 'C00042', 'C00122', 'C00026', 'C00036', 'C00149', 'C00022'],  # Glyoxylate
                'compound_names': ['Glyoxylate', 'Oxaloacetate', 'Fumarate', 'Oxoglutarate', 'Succinate', 'Malate', 'Pyruvate'],
                'enzymes': ['EC:4.1.3.1', 'EC:1.1.1.26', 'EC:1.1.1.37', 'EC:4.2.1.2', 'EC:1.2.4.1', 'EC:2.6.1.4'],
                'enzyme_names': ['Isocitrate lyase', 'Glyoxylate reductase', 'Malate dehydrogenase', 'Fumarase', 'Pyruvate dehydrogenase', 'Alanine-glyoxylate transaminase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Glyoxylate formation', 'Glycolate formation', 'Oxalate formation', 'TCA cycle integration', 'Energy production', 'Amino acid metabolism'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -38.9,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00630',
                'clinical_significance': 'Specialized metabolic bypass pathway',
                'affected_organs': ['Liver', 'Kidney'],
                'biomarkers': ['Glyoxylate', 'Oxalate']
            },
            'hsa00670': {
                'name': 'One Carbon Pool by Folate',
                'compounds': ['C00101', 'C00019', 'C00037', 'C00041', 'C00064', 'C00135', 'C00234'],  # Folate
                'compound_names': ['Tetrahydrofolate', 'Glycine', 'Formate', 'Alanine', 'Methionine', 'Histidine', 'Serine'],
                'enzymes': ['EC:1.5.1.5', 'EC:2.1.2.1', 'EC:2.1.2.3', 'EC:2.1.2.5', 'EC:2.1.2.10', 'EC:2.1.2.14'],
                'enzyme_names': ['Methylenetetrahydrofolate reductase', 'Glycine hydroxymethyltransferase', 'Formiminotransferase', 'Aminomethyltransferase', 'Methionine synthase', 'Serine hydroxymethyltransferase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Methyl group transfer', 'Formate activation', 'Histidine metabolism', 'Methionine synthesis', 'Purine synthesis', 'Thymidylate synthesis'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible'],
                'deltaG': -28.5,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00670',
                'clinical_significance': 'One-carbon metabolism and nucleotide synthesis',
                'affected_organs': ['Liver', 'Bone marrow'],
                'biomarkers': ['Homocysteine', 'Folate', 'Vitamin B12']
            },
            'hsa00260': {
                'name': 'Glycine, Serine and Threonine Metabolism',
                'compounds': ['C00037', 'C00065', 'C00188', 'C00143', 'C00041', 'C00064', 'C00101'],  # Glycine
                'compound_names': ['Glycine', 'Serine', 'Threonine', 'Sarcosine', 'Alanine', 'Methionine', 'Tetrahydrofolate'],
                'enzymes': ['EC:2.1.2.1', 'EC:4.2.1.20', 'EC:1.4.3.3', 'EC:2.6.1.44', 'EC:2.1.2.10', 'EC:4.1.2.5'],
                'enzyme_names': ['Glycine hydroxymethyltransferase', 'Serine dehydratase', 'D-amino acid oxidase', 'Serine-pyruvate transaminase', 'Methionine synthase', 'Threonine aldolase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Glycine synthesis', 'Serine synthesis', 'Threonine metabolism', 'One-carbon transfer', 'Methionine synthesis', 'Purine synthesis'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -32.7,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00260',
                'clinical_significance': 'Amino acid metabolism and one-carbon transfer',
                'affected_organs': ['Liver', 'Kidney'],
                'biomarkers': ['Glycine', 'Serine', 'Homocysteine']
            },
            'hsa00270': {
                'name': 'Cysteine and Methionine Metabolism',
                'compounds': ['C00041', 'C00073', 'C00019', 'C00065', 'C00037', 'C00155', 'C00051'],  # Cysteine
                'compound_names': ['Methionine', 'Homocysteine', 'Glycine', 'Serine', 'Cysteine', 'Glutathione', 'Glutathione disulfide'],
                'enzymes': ['EC:2.5.1.47', 'EC:2.1.1.13', 'EC:4.4.1.1', 'EC:2.6.1.3', 'EC:2.5.1.18', 'EC:1.8.1.7'],
                'enzyme_names': ['Cystathionine γ-synthase', 'Methionine synthase', 'Cystathionine β-lyase', 'Cysteine transaminase', 'Glutathione transferase', 'Glutathione reductase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Transsulfuration', 'Methionine salvage', 'Cysteine synthesis', 'Glutathione synthesis', 'Detoxification', 'Antioxidant defense'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'reversible'],
                'deltaG': -45.3,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00270',
                'clinical_significance': 'Sulfur amino acid metabolism and antioxidant defense',
                'affected_organs': ['Liver', 'Kidney'],
                'biomarkers': ['Homocysteine', 'Cystathionine', 'Glutathione']
            },
            'hsa00280': {
                'name': 'Valine, Leucine and Isoleucine Degradation',
                'compounds': ['C00183', 'C00123', 'C00407', 'C00141', 'C00036', 'C00026', 'C00022'],  # Branched-chain amino acids
                'compound_names': ['Valine', 'Leucine', 'Isoleucine', '3-Methyl-2-oxobutanoate', 'Succinate', 'Oxoglutarate', 'Pyruvate'],
                'enzymes': ['EC:2.6.1.42', 'EC:1.2.4.4', 'EC:1.3.99.3', 'EC:2.3.1.9', 'EC:1.1.1.35', 'EC:4.2.1.17'],
                'enzyme_names': ['Branched-chain-amino-acid transaminase', 'Branched-chain α-keto acid dehydrogenase', 'Acyl-CoA dehydrogenase', 'Acetyl-CoA acyltransferase', '3-hydroxyacyl-CoA dehydrogenase', 'Enoyl-CoA hydratase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Transamination', 'Oxidative decarboxylation', 'β-oxidation', 'TCA cycle integration', 'Energy production', 'Ketone body formation'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'reversible', 'irreversible', 'reversible'],
                'deltaG': -68.2,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00280',
                'clinical_significance': 'Branched-chain amino acid catabolism',
                'affected_organs': ['Liver', 'Muscle'],
                'biomarkers': ['Valine', 'Leucine', 'Isoleucine']
            },
            'hsa00300': {
                'name': 'Lysine Biosynthesis',
                'compounds': ['C00047', 'C00134', 'C00449', 'C00026', 'C00042', 'C00122', 'C00036'],  # Lysine
                'compound_names': ['Lysine', 'Diaminopimelate', 'Tetrahydrodipicolinate', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:1.1.1.86', 'EC:2.6.1.17', 'EC:4.1.1.20', 'EC:2.3.1.117', 'EC:1.2.1.31', 'EC:1.5.1.8'],
                'enzyme_names': ['Dihydrodipicolinate reductase', 'Diaminopimelate epimerase', 'Diaminopimelate decarboxylase', 'Tetrahydrodipicolinate acetyltransferase', 'Succinyl-diaminopimelate transaminase', 'Saccharopine dehydrogenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Aspartate semialdehyde formation', 'Dihydrodipicolinate synthesis', 'Tetrahydrodipicolinate synthesis', 'Diaminopimelate synthesis', 'Lysine synthesis', 'Saccharopine pathway'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -52.8,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00300',
                'clinical_significance': 'Essential amino acid biosynthesis',
                'affected_organs': ['Liver', 'Bacteria'],
                'biomarkers': ['Lysine', 'Pipecolate']
            },
            'hsa00310': {
                'name': 'Lysine Degradation',
                'compounds': ['C00047', 'C00449', 'C00134', 'C00026', 'C00042', 'C00122', 'C00036'],  # Lysine
                'compound_names': ['Lysine', 'Saccharopine', '2-Aminoadipate', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:1.5.1.8', 'EC:1.5.1.9', 'EC:2.6.1.39', 'EC:1.2.1.31', 'EC:4.1.1.20', 'EC:1.2.1.95'],
                'enzyme_names': ['Saccharopine dehydrogenase', 'Pipecolate oxidase', '2-Aminoadipate transaminase', 'Glutaryl-CoA dehydrogenase', '2-Aminoadipate decarboxylase', '2-Oxoadipate dehydrogenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Saccharopine pathway', 'Pipecolate pathway', '2-Aminoadipate pathway', 'Glutaryl-CoA formation', 'Acetyl-CoA formation', 'TCA cycle integration'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -48.3,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00310',
                'clinical_significance': 'Lysine catabolism and acetyl-CoA production',
                'affected_organs': ['Liver', 'Brain'],
                'biomarkers': ['Lysine', 'Saccharopine']
            },
            'hsa00330': {
                'name': 'Arginine and Proline Metabolism',
                'compounds': ['C00077', 'C00148', 'C00149', 'C00026', 'C00042', 'C00122', 'C00036'],  # Arginine
                'compound_names': ['Arginine', 'Proline', 'Glutamate', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:1.14.13.39', 'EC:1.5.1.2', 'EC:4.2.1.11', 'EC:3.5.3.1', 'EC:2.6.1.13', 'EC:1.2.1.41'],
                'enzyme_names': ['Nitric oxide synthase', 'Pyrroline-5-carboxylate reductase', 'Pyrroline-5-carboxylate dehydrogenase', 'Arginase', 'Ornithine aminotransferase', 'Glutamate-5-semialdehyde dehydrogenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Nitric oxide synthesis', 'Proline synthesis', 'Ornithine synthesis', 'Glutamate synthesis', 'TCA cycle integration', 'Urea cycle connection'],
                'reaction_types': ['irreversible', 'reversible', 'irreversible', 'irreversible', 'reversible', 'reversible'],
                'deltaG': -42.7,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00330',
                'clinical_significance': 'Amino acid metabolism and nitric oxide production',
                'affected_organs': ['Liver', 'Kidney'],
                'biomarkers': ['Arginine', 'Citrulline', 'Ornithine']
            },
            'hsa00340': {
                'name': 'Histidine Metabolism',
                'compounds': ['C00135', 'C00785', 'C00164', 'C00026', 'C00042', 'C00122', 'C00036'],  # Histidine
                'compound_names': ['Histidine', 'Urocanate', 'Glutamate', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:4.3.1.3', 'EC:1.4.3.16', 'EC:3.5.2.7', 'EC:2.6.1.38', 'EC:1.1.1.23', 'EC:4.2.1.49'],
                'enzyme_names': ['Histidine ammonia-lyase', 'Urocanate hydratase', 'Imidazolonepropionase', 'Formiminoglutamate transferase', 'Formiminoglutamate dehydrogenase', 'Glutamate formiminotransferase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Histidine degradation', 'Urocanate formation', 'Formiminoglutamate formation', 'Glutamate formation', 'TCA cycle integration', 'One-carbon transfer'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'reversible'],
                'deltaG': -38.5,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00340',
                'clinical_significance': 'Histidine catabolism and one-carbon transfer',
                'affected_organs': ['Liver', 'Kidney'],
                'biomarkers': ['Histidine', 'Formiminoglutamate']
            },
            'hsa00350': {
                'name': 'Tyrosine Metabolism',
                'compounds': ['C00082', 'C00078', 'C00135', 'C00026', 'C00042', 'C00122', 'C00036'],  # Tyrosine
                'compound_names': ['Tyrosine', 'Phenylalanine', 'Homogentisate', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:1.13.11.27', 'EC:4.2.1.5', 'EC:1.13.11.5', 'EC:2.6.1.5', 'EC:1.10.3.1', 'EC:1.14.16.2'],
                'enzyme_names': ['Tyrosine aminotransferase', 'Homogentisate 1,2-dioxygenase', '4-Hydroxyphenylpyruvate dioxygenase', 'Aromatic-amino-acid transaminase', 'Tyrosinase', 'Tyrosine hydroxylase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Tyrosine transamination', 'Homogentisate formation', 'Maleylacetoacetate formation', 'Fumarate formation', 'TCA cycle integration', 'Catecholamine synthesis'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -45.2,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00350',
                'clinical_significance': 'Tyrosine catabolism and neurotransmitter synthesis',
                'affected_organs': ['Liver', 'Brain'],
                'biomarkers': ['Tyrosine', 'Homogentisate']
            },
            'hsa00360': {
                'name': 'Phenylalanine Metabolism',
                'compounds': ['C00078', 'C00082', 'C00166', 'C00026', 'C00042', 'C00122', 'C00036'],  # Phenylalanine
                'compound_names': ['Phenylalanine', 'Tyrosine', 'Phenylpyruvate', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:1.14.16.1', 'EC:4.2.1.91', 'EC:1.13.11.27', 'EC:2.6.1.5', 'EC:1.2.1.20', 'EC:1.3.1.12'],
                'enzyme_names': ['Phenylalanine hydroxylase', 'Phenylpyruvate tautomerase', 'Tyrosine aminotransferase', 'Aromatic-amino-acid transaminase', 'Phenylpyruvate dehydrogenase', 'Dihydropteridine reductase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Phenylalanine hydroxylation', 'Tyrosine formation', 'Phenylpyruvate formation', 'Fumarate formation', 'TCA cycle integration', 'Tetrahydrobiopterin regeneration'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'reversible'],
                'deltaG': -42.8,
                'disease': 'Phenylketonuria',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00360',
                'clinical_significance': 'Phenylalanine catabolism and tyrosine synthesis',
                'affected_organs': ['Liver', 'Brain'],
                'biomarkers': ['Phenylalanine', 'Phenylpyruvate']
            },
            'hsa00380': {
                'name': 'Tryptophan Metabolism',
                'compounds': ['C00078', 'C00082', 'C00166', 'C00026', 'C00042', 'C00122', 'C00036'],  # Tryptophan
                'compound_names': ['Tryptophan', 'Kynurenine', 'Serotonin', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:1.13.11.11', 'EC:4.2.1.20', 'EC:1.14.13.9', 'EC:2.6.1.7', 'EC:4.1.1.45', 'EC:1.13.11.52'],
                'enzyme_names': ['Tryptophan 2,3-dioxygenase', 'Kynureninase', 'Indoleamine 2,3-dioxygenase', 'Kynurenine-oxoglutarate transaminase', 'Aromatic-L-amino-acid decarboxylase', 'Tryptophan 5-hydroxylase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Kynurenine pathway', 'Serotonin synthesis', 'Melatonin synthesis', 'NAD+ synthesis', 'TCA cycle integration', 'Neurotransmitter synthesis'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -38.7,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00380',
                'clinical_significance': 'Tryptophan catabolism and neurotransmitter synthesis',
                'affected_organs': ['Liver', 'Brain'],
                'biomarkers': ['Tryptophan', 'Kynurenine']
            },
            'hsa00400': {
                'name': 'Phenylalanine, Tyrosine and Tryptophan Biosynthesis',
                'compounds': ['C00078', 'C00082', 'C00078', 'C00026', 'C00042', 'C00122', 'C00036'],  # Aromatic amino acids
                'compound_names': ['Phenylalanine', 'Tyrosine', 'Tryptophan', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:2.5.1.19', 'EC:4.2.3.4', 'EC:4.2.1.51', 'EC:1.1.1.25', 'EC:2.6.1.57', 'EC:4.1.3.27'],
                'enzyme_names': ['Chorismate mutase', 'Anthranilate synthase', 'Prephenate dehydratase', 'Arogenate dehydrogenase', 'Aromatic-amino-acid transaminase', 'Chorismate synthase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Shikimate pathway', 'Chorismate formation', 'Prephenate formation', 'Arogenate formation', 'Phenylalanine formation', 'Tyrosine formation'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible'],
                'deltaG': -52.3,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00400',
                'clinical_significance': 'Essential amino acid biosynthesis',
                'affected_organs': ['Liver', 'Bacteria'],
                'biomarkers': ['Phenylalanine', 'Tyrosine']
            },
            'hsa00410': {
                'name': 'β-Alanine Metabolism',
                'compounds': ['C00099', 'C00136', 'C00152', 'C00026', 'C00042', 'C00122', 'C00036'],  # β-Alanine
                'compound_names': ['β-Alanine', 'Malonate semialdehyde', 'Uracil', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:3.5.1.6', 'EC:1.2.1.18', 'EC:2.6.1.19', 'EC:4.1.1.15', 'EC:1.2.1.11', 'EC:1.1.1.61'],
                'enzyme_names': ['β-Ureidopropionase', 'Malonate-semialdehyde dehydrogenase', 'β-Alanine-pyruvate transaminase', 'Aspartate 1-decarboxylase', 'Acetaldehyde dehydrogenase', '3-Hydroxyisobutyrate dehydrogenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Uracil degradation', 'Malonate semialdehyde formation', 'β-Alanine formation', 'Carnosine synthesis', 'Pantothenate synthesis', 'TCA cycle integration'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -35.8,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00410',
                'clinical_significance': 'β-Alanine metabolism and carnosine synthesis',
                'affected_organs': ['Liver', 'Muscle'],
                'biomarkers': ['β-Alanine', 'Carnosine']
            },
            'hsa00430': {
                'name': 'Taurine and Hypotaurine Metabolism',
                'compounds': ['C00245', 'C00519', 'C00127', 'C00026', 'C00042', 'C00122', 'C00036'],  # Taurine
                'compound_names': ['Taurine', 'Hypotaurine', 'Cysteine', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:1.13.11.20', 'EC:1.8.1.3', 'EC:2.6.1.3', 'EC:2.3.1.30', 'EC:4.2.1.65', 'EC:1.14.11.20'],
                'enzyme_names': ['Cysteine dioxygenase', 'Hypotaurine dehydrogenase', 'Cysteine transaminase', 'Cysteamine dioxygenase', 'Cysteine sulfinate decarboxylase', 'Taurine dioxygenase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Cysteine oxidation', 'Hypotaurine formation', 'Taurine formation', 'Bile acid conjugation', 'Neurotransmitter function', 'Osmolyte synthesis'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -28.9,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00430',
                'clinical_significance': 'Taurine metabolism and bile acid conjugation',
                'affected_organs': ['Liver', 'Brain'],
                'biomarkers': ['Taurine', 'Hypotaurine']
            },
            'hsa00450': {
                'name': 'Selenocompound Metabolism',
                'compounds': ['C00065', 'C00097', 'C00123', 'C00026', 'C00042', 'C00122', 'C00036'],  # Selenium
                'compound_names': ['Selenocysteine', 'Selenomethionine', 'Selenite', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:2.9.1.1', 'EC:2.9.1.2', 'EC:1.8.3.4', 'EC:1.8.1.4', 'EC:2.5.1.146', 'EC:4.4.1.16'],
                'enzyme_names': ['Selenocysteine synthase', 'Selenophosphate synthetase', 'Selenocysteine lyase', 'Thioredoxin reductase', 'Selenocysteine methyltransferase', 'Selenocysteine β-lyase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Selenocysteine synthesis', 'Selenoprotein synthesis', 'Selenomethionine metabolism', 'Antioxidant function', 'Detoxification', 'Thyroid hormone metabolism'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -32.5,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00450',
                'clinical_significance': 'Selenium metabolism and selenoprotein synthesis',
                'affected_organs': ['Liver', 'Thyroid'],
                'biomarkers': ['Selenium', 'Glutathione peroxidase']
            },
            'hsa00460': {
                'name': 'Cyanoamino Acid Metabolism',
                'compounds': ['C00065', 'C00097', 'C00123', 'C00026', 'C00042', 'C00122', 'C00036'],  # Cyanoamino acids
                'compound_names': ['β-Cyanoalanine', 'Cyanide', 'Nitrilase', 'Oxoglutarate', 'Oxaloacetate', 'Fumarate', 'Succinate'],
                'enzymes': ['EC:4.4.1.9', 'EC:3.5.5.1', 'EC:1.4.3.5', 'EC:2.6.1.64', 'EC:4.2.1.65', 'EC:4.3.1.9'],
                'enzyme_names': ['β-Cyanoalanine synthase', 'Nitrilase', 'Cyanide dioxygenase', 'β-Cyanoalanine transaminase', 'Cyanide hydratase', 'Cyanase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Cyanide detoxification', 'β-Cyanoalanine formation', 'Asparagine formation', 'Ammonia release', 'TCA cycle integration', 'Cyanate degradation'],
                'reaction_types': ['irreversible', 'irreversible', 'irreversible', 'irreversible', 'irreversible', 'reversible'],
                'deltaG': -25.8,
                'disease': 'Metabolic Disorders',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00460',
                'clinical_significance': 'Cyanide detoxification and asparagine synthesis',
                'affected_organs': ['Liver', 'Kidney'],
                'biomarkers': ['Cyanide', 'Thiocyanate']
            },
            'hsa00480': {
                'name': 'Glutathione Metabolism',
                'compounds': ['C00051', 'C00025', 'C00127', 'C00669', 'C00167', 'C00097', 'C00094'],  # Glutathione
                'compound_names': ['Glutathione', 'Glutamate', 'Cysteine', 'Oxidized glutathione', 'Ascorbate', 'Cystine', 'Sulfite'],
                'enzymes': ['EC:1.11.1.9', 'EC:1.8.1.7', 'EC:2.5.1.18', 'EC:1.8.5.1', 'EC:1.8.4.2', 'EC:1.8.3.2'],
                'enzyme_names': ['Glutathione peroxidase', 'Glutathione reductase', 'Glutathione transferase', 'Glutathione dehydrogenase', 'Glutathione-disulfide reductase', 'Sulfite oxidase'],
                'reactions': ['R00256', 'R01046', 'R01049', 'R01061', 'R01070', 'R01512'],
                'reaction_names': ['Glutathione oxidation', 'Reduction', 'Conjugation', 'Detoxification', 'Disulfide reduction', 'Sulfite oxidation'],
                'reaction_types': ['reversible', 'irreversible', 'irreversible', 'irreversible', 'reversible', 'irreversible'],
                'deltaG': -42.8,
                'disease': 'Oxidative Stress',
                'kegg_map_url': 'https://www.genome.jp/kegg-bin/show_pathway?hsa00480',
                'clinical_significance': 'Antioxidant defense system',
                'affected_organs': ['Liver', 'Red blood cells'],
                'biomarkers': ['GSH', 'MDA', 'Oxidized glutathione']
            }
        }
    }

@st.cache_data
def get_reference_ranges():
    return {
        'glucose': {'min': 70, 'max': 100, 'unit': 'mg/dL'},
        'lactate': {'min': 0.5, 'max': 2.2, 'unit': 'mmol/L'},
        'pyruvate': {'min': 0.03, 'max': 0.08, 'unit': 'mmol/L'},
        'urea': {'min': 10, 'max': 50, 'unit': 'mg/dL'},
        'creatinine': {'min': 0.7, 'max': 1.3, 'unit': 'mg/dL'},
        'bun': {'min': 6, 'max': 24, 'unit': 'mg/dL'},
        'ammonia': {'min': 15, 'max': 45, 'unit': 'μg/dL'},
        'alt': {'min': 7, 'max': 56, 'unit': 'U/L'},
        'ast': {'min': 10, 'max': 40, 'unit': 'U/L'},
        'bilirubin': {'min': 0.1, 'max': 1.2, 'unit': 'mg/dL'},
        'vitamin_d': {'min': 30, 'max': 100, 'unit': 'ng/mL'},
        'pth': {'min': 10, 'max': 65, 'unit': 'pg/mL'},
        'tsh': {'min': 0.4, 'max': 4.0, 'unit': 'mIU/L'},
        'cortisol': {'min': 5, 'max': 25, 'unit': 'μg/dL'},
        'insulin': {'min': 2.6, 'max': 24.9, 'unit': 'μIU/mL'},
        'c_peptide': {'min': 1.1, 'max': 4.4, 'unit': 'ng/mL'},
        'apo_a': {'min': 110, 'max': 180, 'unit': 'mg/dL'},
        'apo_b': {'min': 60, 'max': 130, 'unit': 'mg/dL'},
        'total_cholesterol': {'min': 0, 'max': 200, 'unit': 'mg/dL'},
        'ldl': {'min': 0, 'max': 100, 'unit': 'mg/dL'},
        'hdl': {'min': 40, 'max': 200, 'unit': 'mg/dL'},
        'triglycerides': {'min': 0, 'max': 150, 'unit': 'mg/dL'},
        'lipoprotein_a': {'min': 0, 'max': 30, 'unit': 'mg/dL'},
        'magnesium': {'min': 1.7, 'max': 2.2, 'unit': 'mg/dL'},
        'phosphorus': {'min': 2.5, 'max': 4.5, 'unit': 'mg/dL'},
        'uric_acid': {'min': 3.4, 'max': 7.0, 'unit': 'mg/dL'},
        'anti_tpo': {'min': 0, 'max': 9, 'unit': 'IU/mL'},
        'ana': {'min': 0, 'max': 1, 'unit': 'titer'},
        'anti_ccp': {'min': 0, 'max': 20, 'unit': 'U/mL'},
        'fgf_23': {'min': 0, 'max': 100, 'unit': 'RU/mL'},
        'homocysteine': {'min': 4, 'max': 15, 'unit': 'μmol/L'},
        'glutamine': {'min': 390, 'max': 650, 'unit': 'μmol/L'},
        'glutamate': {'min': 10, 'max': 50, 'unit': 'μmol/L'},
        'arginine': {'min': 50, 'max': 150, 'unit': 'μmol/L'},
        'citrulline': {'min': 15, 'max': 50, 'unit': 'μmol/L'},
        'zinc': {'min': 70, 'max': 120, 'unit': 'μg/dL'},
        'selenium': {'min': 70, 'max': 150, 'unit': 'μg/L'},
        'copper': {'min': 70, 'max': 140, 'unit': 'μg/dL'},
        'chromium': {'min': 0.05, 'max': 0.5, 'unit': 'μg/L'},
        'manganese': {'min': 0.6, 'max': 2.3, 'unit': 'μg/L'},
        'gsh': {'min': 3.8, 'max': 5.5, 'unit': 'μmol/L'},
        'mda': {'min': 0, 'max': 1.5, 'unit': 'μmol/L'},
        'free_t3': {'min': 2.3, 'max': 4.2, 'unit': 'pg/mL'},
        'free_t4': {'min': 0.8, 'max': 1.8, 'unit': 'ng/dL'},
        'estrogen': {'min': 15, 'max': 350, 'unit': 'pg/mL'},
        'progesterone': {'min': 0.1, 'max': 25, 'unit': 'ng/mL'},
        'testosterone': {'min': 2.5, 'max': 8.5, 'unit': 'ng/mL'},
        'crp': {'min': 0, 'max': 3, 'unit': 'mg/L'},
        'il_6': {'min': 0, 'max': 5, 'unit': 'pg/mL'},
        'tumor_necrosis_factor_alpha': {'min': 0, 'max': 8.1, 'unit': 'pg/mL'},
        'interferon_gamma': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'leptin': {'min': 0.5, 'max': 15.2, 'unit': 'ng/mL'},
        'adiponectin': {'min': 3, 'max': 30, 'unit': 'μg/mL'},
        'resistin': {'min': 4, 'max': 12, 'unit': 'ng/mL'},
        'ghrelin': {'min': 500, 'max': 1500, 'unit': 'pg/mL'},
        'glucagon': {'min': 50, 'max': 150, 'unit': 'pg/mL'},
        'gip': {'min': 8, 'max': 42, 'unit': 'pmol/L'},
        'glp_1': {'min': 5, 'max': 20, 'unit': 'pmol/L'},
        'amylin': {'min': 2, 'max': 15, 'unit': 'pmol/L'},
        'somatostatin': {'min': 10, 'max': 25, 'unit': 'pg/mL'},
        'neuropeptide_y': {'min': 20, 'max': 100, 'unit': 'pg/mL'},
        'peptide_yy': {'min': 5, 'max': 30, 'unit': 'pg/mL'},
        'oxyntomodulin': {'min': 5, 'max': 25, 'unit': 'pmol/L'},
        'cholecystokinin': {'min': 0.5, 'max': 5, 'unit': 'pmol/L'},
        'motilin': {'min': 50, 'max': 300, 'unit': 'pg/mL'},
        'vasoactive_intestinal_peptide': {'min': 0, 'max': 30, 'unit': 'pg/mL'},
        'substance_p': {'min': 0, 'max': 100, 'unit': 'pg/mL'},
        'neurokinin_a': {'min': 0, 'max': 30, 'unit': 'pg/mL'},
        'bombesin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'gastrin_releasing_peptide': {'min': 0, 'max': 15, 'unit': 'pg/mL'},
        'uroguanylin': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'guanylin': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'endothelin_1': {'min': 0, 'max': 5, 'unit': 'pg/mL'},
        'calcitonin_gene_related_peptide': {'min': 0, 'max': 50, 'unit': 'pg/mL'},
        'adrenomedullin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'natriuretic_peptide_brain': {'min': 0, 'max': 50, 'unit': 'pg/mL'},
        'natriuretic_peptide_c': {'min': 0, 'max': 100, 'unit': 'pg/mL'},
        'relaxin': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'urotensin_ii': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'apelin': {'min': 0, 'max': 100, 'unit': 'pg/mL'},
        'orexin_a': {'min': 0, 'max': 30, 'unit': 'pg/mL'},
        'orexin_b': {'min': 0, 'max': 30, 'unit': 'pg/mL'},
        'nesfatin_1': {'min': 0, 'max': 10, 'unit': 'ng/mL'},
        'spexin': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'kisspeptin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'phoenixin': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'galanin': {'min': 0, 'max': 50, 'unit': 'pg/mL'},
        'melanin_concentrating_hormone': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'corticotropin_releasing_hormone': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'thyrotropin_releasing_hormone': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'growth_hormone_releasing_hormone': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'somatostatin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'gonadotropin_releasing_hormone': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'prolactin_releasing_peptide': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'prolactin_inhibiting_factor': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'melanocyte_stimulating_hormone': {'min': 0, 'max': 10, 'unit': 'pg/mL'},
        'beta_endorphin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'met_enkephalin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'leu_enkephalin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'dynorphin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'nociceptin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'orphanin_fq': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'endomorphin_1': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'endomorphin_2': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'hemopressin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'neurotensin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'neuromedin_n': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'neuromedin_b': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'neuromedin_u': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'neuromedin_s': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'neuromedin_k': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'neuromedin_c': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'secretin': {'min': 0, 'max': 20, 'unit': 'pg/mL'},
        'vasopressin': {'min': 0, 'max': 5, 'unit': 'pg/mL'},
        'oxytocin': {'min': 0, 'max': 5, 'unit': 'pg/mL'},
        'angiotensin_ii': {'min': 0, 'max': 25, 'unit': 'pg/mL'},
        'aldosterone': {'min': 0, 'max': 30, 'unit': 'pg/mL'},
        'renin': {'min': 0, 'max': 30, 'unit': 'pg/mL'},
        'erythropoietin': {'min': 0, 'max': 30, 'unit': 'mIU/mL'},
        'thrombopoietin': {'min': 0, 'max': 200, 'unit': 'pg/mL'},
        'leptin': {'min': 0.5, 'max': 15.2, 'unit': 'ng/mL'},
        'adiponectin': {'min': 3, 'max': 30, 'unit': 'μg/mL'},
        'resistin': {'min': 4, 'max': 12, 'unit': 'ng/mL'},
        'ghrelin': {'min': 500, 'max': 1500, 'unit': 'pg/mL'},
        'obestatin': {'min': 0, 'max': 100, 'unit': 'pg/mL'},
        'nesfatin_1': {'min': 0, 'max': 10, 'unit': 'ng/mL'},
        'vaspin': {'min': 0, 'max': 10, 'unit': 'ng/mL'},
        'omentin': {'min': 0, 'max': 100, 'unit': 'ng/mL'},
        'chemerin': {'min': 0, 'max': 200, 'unit': 'ng/mL'},
        'apelin': {'min': 0, 'max': 100, 'unit': 'pg/mL'},
        'visfatin': {'min': 0, 'max': 50, 'unit': 'ng/mL'},
        'progranulin': {'min': 0, 'max': 100, 'unit': 'ng/mL'},
        'fibroblast_growth_factor_21': {'min': 0, 'max': 300, 'unit': 'pg/mL'},
        'irisin': {'min': 0, 'max': 100, 'unit': 'ng/mL'},
        'metrnl': {'min': 0, 'max': 100, 'unit': 'ng/mL'},
        'asprosin': {'min': 0, 'max': 10, 'unit': 'ng/mL'},
        'lipocalin_2': {'min': 0, 'max': 100, 'unit': 'ng/mL'},
        'retinol_binding_protein_4': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'adipsin': {'min': 0, 'max': 10, 'unit': 'μg/mL'},
        'adipocyte_fatty_acid_binding_protein': {'min': 0, 'max': 100, 'unit': 'ng/mL'},
        'angiopoietin_like_protein_4': {'min': 0, 'max': 100, 'unit': 'ng/mL'},
        'angiopoietin_like_protein_8': {'min': 0, 'max': 100, 'unit': 'ng/mL'},
        'fetuin_a': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'fetuin_b': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'selenoprotein_p': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'zinc_alpha_2_glycoprotein': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'afamin': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'alpha_1_acid_glycoprotein': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'alpha_1_antitrypsin': {'min': 0, 'max': 200, 'unit': 'μg/mL'},
        'alpha_2_macroglobulin': {'min': 0, 'max': 300, 'unit': 'μg/mL'},
        'antithrombin_iii': {'min': 0, 'max': 150, 'unit': 'μg/mL'},
        'apolipoprotein_a1': {'min': 0, 'max': 200, 'unit': 'μg/mL'},
        'apolipoprotein_a2': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_a4': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_a5': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_b': {'min': 0, 'max': 150, 'unit': 'μg/mL'},
        'apolipoprotein_c1': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_c2': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_c3': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_c4': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_d': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_e': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_h': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_j': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_l1': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_m': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'apolipoprotein_o': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'beta_2_microglobulin': {'min': 0, 'max': 3, 'unit': 'μg/mL'},
        'c_reactive_protein': {'min': 0, 'max': 10, 'unit': 'μg/mL'},
        'ceruloplasmin': {'min': 0, 'max': 60, 'unit': 'μg/mL'},
        'complement_c3': {'min': 0, 'max': 180, 'unit': 'μg/mL'},
        'complement_c4': {'min': 0, 'max': 60, 'unit': 'μg/mL'},
        'corticosteroid_binding_globulin': {'min': 0, 'max': 50, 'unit': 'μg/mL'},
        'cortisol': {'min': 0, 'max': 25, 'unit': 'μg/dL'},
        'cystatin_c': {'min': 0, 'max': 1.5, 'unit': 'μg/mL'},
        'd_dimer': {'min': 0, 'max': 0.5, 'unit': 'μg/mL'},
        'erythropoietin': {'min': 0, 'max': 30, 'unit': 'mIU/mL'},
        'factor_viii': {'min': 50, 'max': 150, 'unit': '%'},
        'ferritin': {'min': 0, 'max': 300, 'unit': 'ng/mL'},
        'fibrinogen': {'min': 0, 'max': 400, 'unit': 'μg/mL'},
        'folate': {'min': 0, 'max': 20, 'unit': 'ng/mL'},
        'free_thyroxine': {'min': 0, 'max': 2, 'unit': 'ng/dL'},
        'free_triiodothyronine': {'min': 0, 'max': 0.5, 'unit': 'pg/mL'},
        'growth_hormone': {'min': 0, 'max': 10, 'unit': 'ng/mL'},
        'haptoglobin': {'min': 0, 'max': 300, 'unit': 'μg/mL'},
        'hemoglobin_a1c': {'min': 0, 'max': 6, 'unit': '%'},
        'homocysteine': {'min': 0, 'max': 15, 'unit': 'μmol/L'},
        'immunoglobulin_a': {'min': 0, 'max': 400, 'unit': 'μg/mL'},
        'immunoglobulin_g': {'min': 0, 'max': 1600, 'unit': 'μg/mL'},
        'immunoglobulin_m': {'min': 0, 'max': 300, 'unit': 'μg/mL'},
        'insulin': {'min': 0, 'max': 25, 'unit': 'μIU/mL'},
        'insulin_like_growth_factor_1': {'min': 0, 'max': 300, 'unit': 'ng/mL'},
        'insulin_like_growth_factor_binding_protein_3': {'min': 0, 'max': 5000, 'unit': 'ng/mL'},
        'interleukin_6': {'min': 0, 'max': 5, 'unit': 'pg/mL'},
        'iron': {'min': 0, 'max': 200, 'unit': 'μg/dL'},
        'lactate_dehydrogenase': {'min': 0, 'max': 250, 'unit': 'U/L'},
        'lipoprotein_a': {'min': 0, 'max': 30, 'unit': 'mg/dL'},
        'parathyroid_hormone': {'min': 0, 'max': 65, 'unit': 'pg/mL'},
        'prealbumin': {'min': 0, 'max': 40, 'unit': 'μg/mL'},
        'prolactin': {'min': 0, 'max': 25, 'unit': 'ng/mL'},
        'prostate_specific_antigen': {'min': 0, 'max': 4, 'unit': 'ng/mL'},
        'retinol_binding_protein': {'min': 0, 'max': 100, 'unit': 'μg/mL'},
        'sex_hormone_binding_globulin': {'min': 0, 'max': 100, 'unit': 'nmol/L'},
        'testosterone': {'min': 0, 'max': 10, 'unit': 'ng/mL'},
        'transferrin': {'min': 0, 'max': 400, 'unit': 'μg/mL'},
        'transferrin_saturation': {'min': 0, 'max': 50, 'unit': '%'},
        'tumor_necrosis_factor_alpha': {'min': 0, 'max': 8.1, 'unit': 'pg/mL'},
        'vitamin_b12': {'min': 0, 'max': 900, 'unit': 'pg/mL'},
        'vitamin_d': {'min': 0, 'max': 100, 'unit': 'ng/mL'},
        'von_willebrand_factor': {'min': 50, 'max': 150, 'unit': '%'}
    }


def calculate_z_score(value, mean, std):
    return (value - mean) / std

def get_pathway_perturbation_score(lab_values, metabolic_db, ref_ranges):
    scores = {}
    
    for pathway_id, pathway in metabolic_db['pathways'].items():
        perturbation_score = 0
        relevant_values = 0
        abnormal_markers = []
        
        # === METABOLIC PANEL ===
        
        # Glucose - Critical for diabetes pathways
        if 'glucose' in lab_values and lab_values['glucose'] is not None:
            if 'C00031' in pathway['compounds'] or pathway['disease'] == 'Diabetes Mellitus':
                reference = ref_ranges['glucose']
                z_score = calculate_z_score(lab_values['glucose'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 2.0
                relevant_values += 2.0
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Glucose: {'High' if z_score > 0 else 'Low'}")
        
        # Lactate - Metabolic pathways
        if 'lactate' in lab_values and lab_values['lactate'] is not None:
            if 'C00186' in pathway['compounds'] or pathway['disease'] in ['Diabetes Mellitus', 'Metabolic Disorders']:
                reference = ref_ranges['lactate']
                z_score = calculate_z_score(lab_values['lactate'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.3
                relevant_values += 1.3
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Lactate: {'High' if z_score > 0 else 'Low'}")
        
        # Pyruvate - Metabolic pathways
        if 'pyruvate' in lab_values and lab_values['pyruvate'] is not None:
            if 'C00022' in pathway['compounds'] or pathway['disease'] in ['Diabetes Mellitus', 'Metabolic Disorders']:
                reference = ref_ranges['pyruvate']
                z_score = calculate_z_score(lab_values['pyruvate'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.2
                relevant_values += 1.2
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Pyruvate: {'High' if z_score > 0 else 'Low'}")
        
        # Kidney function markers
        kidney_diseases = ['Chronic Kidney Disease']
        
        if 'urea' in lab_values and lab_values['urea'] is not None:
            if pathway['disease'] in kidney_diseases:
                reference = ref_ranges['urea']
                z_score = calculate_z_score(lab_values['urea'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.2
                relevant_values += 1.2
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Urea: {'High' if z_score > 0 else 'Low'}")
        
        if 'creatinine' in lab_values and lab_values['creatinine'] is not None:
            if pathway['disease'] in kidney_diseases:
                reference = ref_ranges['creatinine']
                z_score = calculate_z_score(lab_values['creatinine'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.3
                relevant_values += 1.3
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Creatinine: {'High' if z_score > 0 else 'Low'}")
        
        if 'bun' in lab_values and lab_values['bun'] is not None:
            if pathway['disease'] in kidney_diseases:
                reference = ref_ranges['bun']
                z_score = calculate_z_score(lab_values['bun'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.1
                relevant_values += 1.1
                if abs(z_score) > 2:
                    abnormal_markers.append(f"BUN: {'High' if z_score > 0 else 'Low'}")
        
        # Ammonia - Liver/metabolic pathways
        if 'ammonia' in lab_values and lab_values['ammonia'] is not None:
            if pathway['disease'] in ['Diabetes Mellitus', 'Metabolic Disorders']:
                reference = ref_ranges['ammonia']
                z_score = calculate_z_score(lab_values['ammonia'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.0
                relevant_values += 1.0
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Ammonia: {'High' if z_score > 0 else 'Low'}")
        
        # Liver function markers
        liver_diseases = ['Diabetes Mellitus', 'Dyslipidemia', 'Metabolic Disorders']
        
        if 'alt' in lab_values and lab_values['alt'] is not None:
            if pathway['disease'] in liver_diseases:
                reference = ref_ranges['alt']
                z_score = calculate_z_score(lab_values['alt'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 0.9
                relevant_values += 0.9
                if abs(z_score) > 2:
                    abnormal_markers.append(f"ALT: {'High' if z_score > 0 else 'Low'}")
        
        if 'ast' in lab_values and lab_values['ast'] is not None:
            if pathway['disease'] in liver_diseases:
                reference = ref_ranges['ast']
                z_score = calculate_z_score(lab_values['ast'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 0.9
                relevant_values += 0.9
                if abs(z_score) > 2:
                    abnormal_markers.append(f"AST: {'High' if z_score > 0 else 'Low'}")
        
        if 'bilirubin' in lab_values and lab_values['bilirubin'] is not None:
            if pathway['disease'] in liver_diseases:
                reference = ref_ranges['bilirubin']
                z_score = calculate_z_score(lab_values['bilirubin'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 0.8
                relevant_values += 0.8
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Bilirubin: {'High' if z_score > 0 else 'Low'}")
        
        # === ENDOCRINE PANEL ===
        
        endocrine_diseases = ['Endocrine Disorders', 'Diabetes Mellitus']
        
        # Vitamin D
        if 'vitamin_d' in lab_values and lab_values['vitamin_d'] is not None:
            if pathway['disease'] in endocrine_diseases + ['Chronic Kidney Disease']:
                reference = ref_ranges['vitamin_d']
                z_score = calculate_z_score(lab_values['vitamin_d'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.0
                relevant_values += 1.0
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Vitamin D: {'High' if z_score > 0 else 'Low'}")
        
        # PTH
        if 'pth' in lab_values and lab_values['pth'] is not None:
            if pathway['disease'] in endocrine_diseases + ['Chronic Kidney Disease']:
                reference = ref_ranges['pth']
                z_score = calculate_z_score(lab_values['pth'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.1
                relevant_values += 1.1
                if abs(z_score) > 2:
                    abnormal_markers.append(f"PTH: {'High' if z_score > 0 else 'Low'}")
        
        # TSH
        if 'tsh' in lab_values and lab_values['tsh'] is not None:
            if pathway['disease'] in endocrine_diseases:
                reference = ref_ranges['tsh']
                z_score = calculate_z_score(lab_values['tsh'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.3
                relevant_values += 1.3
                if abs(z_score) > 2:
                    abnormal_markers.append(f"TSH: {'High' if z_score > 0 else 'Low'}")
        
        # Cortisol
        if 'cortisol' in lab_values and lab_values['cortisol'] is not None:
            if pathway['disease'] in endocrine_diseases:
                reference = ref_ranges['cortisol']
                z_score = calculate_z_score(lab_values['cortisol'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.2
                relevant_values += 1.2
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Cortisol: {'High' if z_score > 0 else 'Low'}")
        
        # Insulin - Critical for diabetes
        if 'insulin' in lab_values and lab_values['insulin'] is not None:
            if pathway['disease'] == 'Diabetes Mellitus':
                reference = ref_ranges['insulin']
                z_score = calculate_z_score(lab_values['insulin'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.8
                relevant_values += 1.8
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Insulin: {'High' if z_score > 0 else 'Low'}")
        
        # C-peptide
        if 'c_peptide' in lab_values and lab_values['c_peptide'] is not None:
            if pathway['disease'] == 'Diabetes Mellitus':
                reference = ref_ranges['c_peptide']
                z_score = calculate_z_score(lab_values['c_peptide'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.6
                relevant_values += 1.6
                if abs(z_score) > 2:
                    abnormal_markers.append(f"C-peptide: {'High' if z_score > 0 else 'Low'}")
        
        # Sex hormones
        if 'estrogen' in lab_values and lab_values['estrogen'] is not None:
            if pathway['disease'] in endocrine_diseases:
                reference = ref_ranges['estrogen']
                z_score = calculate_z_score(lab_values['estrogen'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 0.8
                relevant_values += 0.8
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Estrogen: {'High' if z_score > 0 else 'Low'}")
        
        if 'progesterone' in lab_values and lab_values['progesterone'] is not None:
            if pathway['disease'] in endocrine_diseases:
                reference = ref_ranges['progesterone']
                z_score = calculate_z_score(lab_values['progesterone'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 0.8
                relevant_values += 0.8
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Progesterone: {'High' if z_score > 0 else 'Low'}")
        
        if 'testosterone' in lab_values and lab_values['testosterone'] is not None:
            if pathway['disease'] in endocrine_diseases:
                reference = ref_ranges['testosterone']
                z_score = calculate_z_score(lab_values['testosterone'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 0.9
                relevant_values += 0.9
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Testosterone: {'High' if z_score > 0 else 'Low'}")
        
        # === LIPID PANEL ===
        
        lipid_diseases = ['Dyslipidemia', 'Metabolic Syndrome']
        
        # Apo A with validation
        if 'apo_a' in lab_values and lab_values['apo_a'] is not None:
            if lab_values['apo_a'] <= 0:
                print(f"Warning: Invalid APO_A value: {lab_values['apo_a']}")
            elif pathway['disease'] in lipid_diseases:
                reference = ref_ranges['apo_a']
                z_score = calculate_z_score(lab_values['apo_a'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.3
                relevant_values += 1.3
                if abs(z_score) > 2:
                    abnormal_markers.append(f"APO_A: {'High' if z_score > 0 else 'Low'}")
        
        # Apo B
        if 'apo_b' in lab_values and lab_values['apo_b'] is not None:
            if pathway['disease'] in lipid_diseases:
                reference = ref_ranges['apo_b']
                z_score = calculate_z_score(lab_values['apo_b'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.4
                relevant_values += 1.4
                if abs(z_score) > 2:
                    abnormal_markers.append(f"APO_B: {'High' if z_score > 0 else 'Low'}")
        
        # Total Cholesterol
        if 'total_cholesterol' in lab_values and lab_values['total_cholesterol'] is not None:
            if pathway['disease'] in lipid_diseases:
                reference = ref_ranges['total_cholesterol']
                z_score = calculate_z_score(lab_values['total_cholesterol'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.6
                relevant_values += 1.6
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Total cholesterol: {'High' if z_score > 0 else 'Low'}")
        
        # LDL
        if 'ldl' in lab_values and lab_values['ldl'] is not None:
            if pathway['disease'] in lipid_diseases:
                reference = ref_ranges['ldl']
                z_score = calculate_z_score(lab_values['ldl'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.8
                relevant_values += 1.8
                if abs(z_score) > 2:
                    abnormal_markers.append(f"LDL: {'High' if z_score > 0 else 'Low'}")
        
        # HDL
        if 'hdl' in lab_values and lab_values['hdl'] is not None:
            if pathway['disease'] in lipid_diseases:
                reference = ref_ranges['hdl']
                z_score = calculate_z_score(lab_values['hdl'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.5
                relevant_values += 1.5
                if abs(z_score) > 2:
                    abnormal_markers.append(f"HDL: {'High' if z_score > 0 else 'Low'}")
        
        # Triglycerides
        if 'triglycerides' in lab_values and lab_values['triglycerides'] is not None:
            if 'C00162' in pathway['compounds'] or pathway['disease'] in lipid_diseases:
                reference = ref_ranges['triglycerides']
                z_score = calculate_z_score(lab_values['triglycerides'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.6
                relevant_values += 1.6
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Triglycerides: {'High' if z_score > 0 else 'Low'}")
        
        # Lipoprotein(a)
        if 'lipoprotein_a' in lab_values and lab_values['lipoprotein_a'] is not None:
            if pathway['disease'] in lipid_diseases:
                reference = ref_ranges['lipoprotein_a']
                z_score = calculate_z_score(lab_values['lipoprotein_a'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.1
                relevant_values += 1.1
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Lipoprotein(a): {'High' if z_score > 0 else 'Low'}")
        
        # === MINERALS & ELECTROLYTES ===
        
        # Phosphorus
        if 'phosphorus' in lab_values and lab_values['phosphorus'] is not None:
            if 'C00009' in pathway['compounds'] or pathway['disease'] in ['Chronic Kidney Disease', 'Endocrine Disorders']:
                reference = ref_ranges['phosphorus']
                z_score = calculate_z_score(lab_values['phosphorus'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.0
                relevant_values += 1.0
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Phosphorus: {'High' if z_score > 0 else 'Low'}")
        
        # Magnesium
        if 'magnesium' in lab_values and lab_values['magnesium'] is not None:
            if pathway['disease'] in ['Chronic Kidney Disease', 'Endocrine Disorders', 'Diabetes Mellitus']:
                reference = ref_ranges['magnesium']
                z_score = calculate_z_score(lab_values['magnesium'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 0.9
                relevant_values += 0.9
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Magnesium: {'High' if z_score > 0 else 'Low'}")
        
        # Uric Acid
        if 'uric_acid' in lab_values and lab_values['uric_acid'] is not None:
            if pathway['disease'] in ['Metabolic Disorders', 'Chronic Kidney Disease']:
                reference = ref_ranges['uric_acid']
                z_score = calculate_z_score(lab_values['uric_acid'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.0
                relevant_values += 1.0
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Uric Acid: {'High' if z_score > 0 else 'Low'}")
        

        
        # === IMMUNE MARKERS ===
        
        immune_diseases = ['Autoimmune Diseases', 'Inflammatory Diseases']
        
        # Autoimmune markers
        for marker in ['anti_tpo', 'ana', 'anti_ccp']:
            if marker in lab_values and lab_values[marker] is not None:
                if pathway['disease'] in immune_diseases:
                    reference = ref_ranges[marker]
                    z_score = calculate_z_score(lab_values[marker],
                                              (reference['min'] + reference['max']) / 2,
                                              (reference['max'] - reference['min']) / 4)
                    perturbation_score += abs(z_score) * 0.7  # Lower weight in metabolic context
                    relevant_values += 0.7
                    if abs(z_score) > 2:
                        abnormal_markers.append(f"{marker.replace('_', '-').upper()}: {'High' if z_score > 0 else 'Low'}")
        
        # CRP
        if 'crp' in lab_values and lab_values['crp'] is not None:
            if pathway['disease'] in immune_diseases:
                reference = ref_ranges['crp']
                z_score = calculate_z_score(lab_values['crp'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 0.8  # Moderate weight
                relevant_values += 0.8
                if abs(z_score) > 2:
                    abnormal_markers.append(f"CRP: {'High' if z_score > 0 else 'Low'}")
        
        # FGF-23
        if 'fgf_23' in lab_values and lab_values['fgf_23'] is not None:
            if pathway['disease'] in ['Chronic Kidney Disease', 'Endocrine Disorders']:
                reference = ref_ranges['fgf_23']
                z_score = calculate_z_score(lab_values['fgf_23'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.0
                relevant_values += 1.0
                if abs(z_score) > 2:
                    abnormal_markers.append(f"FGF-23: {'High' if z_score > 0 else 'Low'}")
        
        # === OTHER SPECIALIZED TESTS ===
        
        # Homocysteine
        if 'homocysteine' in lab_values and lab_values['homocysteine'] is not None:
            if pathway['disease'] in ['Metabolic Disorders', 'Dyslipidemia']:
                reference = ref_ranges['homocysteine']
                z_score = calculate_z_score(lab_values['homocysteine'],
                                          (reference['min'] + reference['max']) / 2,
                                          (reference['max'] - reference['min']) / 4)
                perturbation_score += abs(z_score) * 1.1
                relevant_values += 1.1
                if abs(z_score) > 2:
                    abnormal_markers.append(f"Homocysteine: {'High' if z_score > 0 else 'Low'}")
        
        # Amino acids
        amino_acids = ['glutamine', 'glutamate']
        for aa in amino_acids:
            if aa in lab_values and lab_values[aa] is not None:
                if pathway['disease'] in ['Metabolic Disorders']:
                    reference = ref_ranges[aa]
                    z_score = calculate_z_score(lab_values[aa],
                                              (reference['min'] + reference['max']) / 2,
                                              (reference['max'] - reference['min']) / 4)
                    perturbation_score += abs(z_score) * 0.9
                    relevant_values += 0.9
                    if abs(z_score) > 2:
                        abnormal_markers.append(f"{aa.title()}: {'High' if z_score > 0 else 'Low'}")
        
        # Oxidative stress markers
        oxidative_markers = ['gsh', 'mda']
        for marker in oxidative_markers:
            if marker in lab_values and lab_values[marker] is not None:
                if pathway['disease'] in ['Oxidative Stress', 'Metabolic Disorders']:
                    reference = ref_ranges[marker]
                    z_score = calculate_z_score(lab_values[marker],
                                              (reference['min'] + reference['max']) / 2,
                                              (reference['max'] - reference['min']) / 4)
                    perturbation_score += abs(z_score) * 1.0
                    relevant_values += 1.0
                    if abs(z_score) > 2:
                        abnormal_markers.append(f"{marker.upper()}: {'High' if z_score > 0 else 'Low'}")
        
        # Calculate final scores
        if relevant_values > 0:
            normalized_perturbation = perturbation_score / relevant_values
            
            # Clinical relevance multiplier
            clinical_multiplier = 1.0
            if pathway['disease'] in ['Diabetes Mellitus', 'Dyslipidemia']:
                clinical_multiplier = 1.5
            elif pathway['disease'] in ['Metabolic Syndrome', 'Metabolic Disorders']:
                clinical_multiplier = 1.3
            elif pathway['disease'] in ['Chronic Kidney Disease', 'Endocrine Disorders']:
                clinical_multiplier = 1.1
            elif pathway['disease'] in ['Autoimmune Diseases', 'Inflammatory Diseases']:
                clinical_multiplier = 0.7
            
            # Thermodynamic score
            thermodynamic_score = min(abs(pathway.get('deltaG', 0)) / 100, 5.0)
            
            total_score = (normalized_perturbation * clinical_multiplier) + thermodynamic_score
            
            scores[pathway_id] = {
                'perturbation_score': normalized_perturbation,
                'clinical_multiplier': clinical_multiplier,
                'thermodynamic_score': thermodynamic_score,
                'total_score': total_score,
                'pathway': pathway,
                'relevant_markers': relevant_values,
                'abnormal_markers': abnormal_markers
            }
    
    return scores

def calculate_z_score(value, mean, std):
    """Calculate z-score with proper error handling"""
    if std == 0:
        return 0
    return (value - mean) / std
def calculate_z_score(value, mean, std):
    """Calculate z-score with proper error handling"""
    if std == 0:
        return 0
    return (value - mean) / std

def get_affected_pathways_analysis(lab_values, metabolic_db, ref_ranges):
    """Analyze which pathways are directly affected and at risk"""
    directly_affected = []
    at_risk = []
    
    pathway_scores = get_pathway_perturbation_score(lab_values, metabolic_db, ref_ranges)
    
    for pathway_id, data in pathway_scores.items():
        if data['perturbation_score'] > 1.5:  # High perturbation
            directly_affected.append({
                'pathway_id': pathway_id,
                'name': data['pathway']['name'],
                'disease': data['pathway']['disease'],
                'score': data['total_score'],
                'affected_biomarkers': get_affected_biomarkers(lab_values, data['pathway'], ref_ranges)
            })
        elif data['perturbation_score'] > 0.5:  # Moderate perturbation
            at_risk.append({
                'pathway_id': pathway_id,
                'name': data['pathway']['name'],
                'disease': data['pathway']['disease'],
                'score': data['total_score'],
                'risk_factors': get_risk_factors(lab_values, data['pathway'], ref_ranges)
            })
    
    return directly_affected, at_risk

def get_affected_biomarkers(lab_values, pathway, ref_ranges):
    """Get biomarkers that are abnormal and related to the pathway"""
    affected = []
    for biomarker in pathway['biomarkers']:
        biomarker_key = biomarker.lower().replace('-', '_').replace(' ', '_')
        if biomarker_key in lab_values and biomarker_key in ref_ranges:
            value = lab_values[biomarker_key]
            ref = ref_ranges[biomarker_key]
            if value < ref['min'] or value > ref['max']:
                status = 'Low' if value < ref['min'] else 'High'
                affected.append(f"{biomarker}: {status}")
    return affected

def get_risk_factors(lab_values, pathway, ref_ranges):
    """Get risk factors for pathways at risk"""
    risk_factors = []
    for biomarker in pathway['biomarkers']:
        biomarker_key = biomarker.lower().replace('-', '_').replace(' ', '_')
        if biomarker_key in lab_values and biomarker_key in ref_ranges:
            value = lab_values[biomarker_key]
            ref = ref_ranges[biomarker_key]
            # Check if value is in upper 75% of normal range (potential risk)
            upper_75 = ref['min'] + 0.75 * (ref['max'] - ref['min'])
            if ref['min'] <= value <= ref['max'] and value > upper_75:
                risk_factors.append(f"{biomarker}: Upper normal range")
    return risk_factors

def create_overview_chart(sorted_pathways):
    pathway_names = [data[1]['pathway']['name'][:20] + '...' if len(data[1]['pathway']['name']) > 20 else data[1]['pathway']['name'] for data in sorted_pathways[:5]]
    pathway_scores = [data[1]['total_score'] for data in sorted_pathways[:5]]
    
    colors = ['#EF4444', '#F56565', '#FB923C', '#22C55E', '#3B82F6']
    
    fig = go.Figure(data=[
        go.Bar(
            y=pathway_names,
            x=pathway_scores,
            orientation='h',
            marker_color=colors,
            text=[f"{score:.2f}" for score in pathway_scores],
            textposition='outside'
        )
    ])
    
    fig.update_layout(
        title="Top 5 Affected Metabolic Pathways",
        xaxis_title="Perturbation Score",
        height=400,
        margin=dict(l=200)
    )
    
    return fig

def create_thermodynamics_chart(sorted_pathways):
    pathway_names = [data[1]['pathway']['name'][:15] + '...' if len(data[1]['pathway']['name']) > 15 else data[1]['pathway']['name'] for data in sorted_pathways]
    delta_g = [data[1]['pathway']['deltaG'] for data in sorted_pathways]
    
    colors = ['#48BB78' if dg < -50 else '#F56565' if dg > -20 else '#718096' for dg in delta_g]
    
    fig = go.Figure(data=[
        go.Bar(
            x=pathway_names,
            y=delta_g,
            marker_color=colors,
            text=[f"{dg:.1f}" for dg in delta_g],
            textposition='outside'
        )
    ])
    
    fig.update_layout(
        title="Thermodynamic Favorability of Metabolic Pathways",
        yaxis_title="Gibbs Free Energy (kJ/mol)",
        height=400
    )
    
    return fig

def create_flux_chart(sorted_pathways):
    pathway_names = [data[1]['pathway']['name'][:12] + '...' if len(data[1]['pathway']['name']) > 12 else data[1]['pathway']['name'] for data in sorted_pathways[:5]]
    
    normal_flux = [100] * 5
    perturbed_flux = [100 * (1 + (data[1]['perturbation_score'] * 0.3)) for data in sorted_pathways[:5]]
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatterpolar(
        r=normal_flux,
        theta=pathway_names,
        fill='toself',
        name='Normal Flux',
        line_color='rgba(72, 187, 120, 1)'
    ))
    
    fig.add_trace(go.Scatterpolar(
        r=perturbed_flux,
        theta=pathway_names,
        fill='toself',
        name='Perturbed Flux',
        line_color='rgba(245, 101, 101, 1)'
    ))
    
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 150]
            )),
        showlegend=True,
        title="Metabolic Flux Analysis - Normal vs Perturbed States",
        height=500
    )
    
    return fig

def create_kegg_style_pathway_chart(pathway_data, selected_pathway_id=None):
    """Create a KEGG-style metabolic pathway visualization with enhanced features"""
    
    if selected_pathway_id and selected_pathway_id in pathway_data:
        pathway = pathway_data[selected_pathway_id]['pathway']
        
        fig = go.Figure()
        
        # Define positions for compounds, enzymes, and reactions
        compounds = pathway['compounds'][:4]
        compound_names = pathway['compound_names'][:4]
        enzymes = pathway['enzymes'][:3]
        enzyme_names = pathway['enzyme_names'][:3]
        reactions = pathway['reactions'][:3]
        reaction_names = pathway['reaction_names'][:3]
        reaction_types = pathway['reaction_types'][:3]
        
        # Positions for KEGG-like layout
        compound_positions = [(0, 2), (2, 3), (4, 2), (6, 1)]
        enzyme_positions = [(1, 1), (3, 1), (5, 1)]
        reaction_positions = [(1, 2.5), (3, 2.5), (5, 2.5)]
        
        # Add compounds (circles with KEGG links)
        for i, (comp_id, comp_name, pos) in enumerate(zip(compounds, compound_names, compound_positions)):
            fig.add_trace(go.Scatter(
                x=[pos[0]], y=[pos[1]],
                mode='markers+text',
                marker=dict(
                    size=40, 
                    color='lightblue', 
                    line=dict(width=2, color='darkblue'),
                    symbol='circle'
                ),
                text=comp_name.split()[0],
                textposition="middle center",
                name=f"Compound",
                hovertext=f"{comp_name} ({comp_id})<br>KEGG: https://www.kegg.jp/dbget-bin/www_bget?cpd:{comp_id}",
                showlegend=False
            ))
        
        # Add enzymes (rectangles with KEGG links)
        for i, (enz_id, enz_name, pos) in enumerate(zip(enzymes, enzyme_names, enzyme_positions)):
            fig.add_trace(go.Scatter(
                x=[pos[0]], y=[pos[1]],
                mode='markers+text',
                marker=dict(
                    size=35, 
                    color='lightgreen', 
                    line=dict(width=2, color='darkgreen'),
                    symbol='square'
                ),
                text=f"E{i+1}",
                textposition="middle center",
                name=f"Enzyme",
                hovertext=f"{enz_name} ({enz_id})<br>KEGG: https://www.kegg.jp/dbget-bin/www_bget?ec:{enz_id.split(':')[1]}",
                showlegend=False
            ))
        
        # Add reaction arrows with reversible/irreversible indicators
        for i, (pos1, pos2, rxn_type) in enumerate(zip(compound_positions[:-1], compound_positions[1:], reaction_types)):
            enz_pos = enzyme_positions[i] if i < len(enzyme_positions) else enzyme_positions[-1]
            
            # Line from compound to enzyme
            fig.add_trace(go.Scatter(
                x=[pos1[0], enz_pos[0]], y=[pos1[1], enz_pos[1]],
                mode='lines',
                line=dict(width=2, color='red'),
                showlegend=False,
                hoverinfo='skip'
            ))
            
            # Line from enzyme to next compound
            if i+1 < len(compound_positions):
                arrow_symbol = 'triangle-right' if rxn_type == 'irreversible' else 'circle'
                fig.add_trace(go.Scatter(
                    x=[enz_pos[0], pos2[0]], y=[enz_pos[1], pos2[1]],
                    mode='lines+markers',
                    line=dict(width=2, color='red'),
                    marker=dict(
                        symbol=arrow_symbol,
                        size=10,
                        color='red',
                        angleref='previous'
                    ),
                    showlegend=False,
                    hovertext=f"Reaction: {reaction_names[i]} ({reactions[i]})<br>Type: {rxn_type}",
                    hoverinfo='text'
                ))
        
        # Add pathway title and info
        fig.update_layout(
            title=f"KEGG-Style View: {pathway['name']}",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-0.5, 6.5]),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[0, 4]),
            height=500,
            annotations=[
                dict(text=f"Disease: {pathway['disease']}", x=0, y=3.5, showarrow=False, font=dict(size=12)),
                dict(text=f"ΔG: {pathway['deltaG']} kJ/mol", x=0, y=3.2, showarrow=False, font=dict(size=12)),
                dict(text="🔵 Compounds (click for KEGG)  🟩 Enzymes (click for KEGG)", x=3, y=0.2, showarrow=False, font=dict(size=10))
            ]
        )
        
    else:
        # Overview network showing pathway connections
        fig = create_pathway_overview_network(pathway_data)
    
    return fig

def create_pathway_overview_network(pathway_data):
    """Create an overview network of pathway relationships"""
    import networkx as nx
    
    # Group pathways by disease
    disease_groups = {}
    for pathway_id, data in pathway_data.items():
        disease = data['pathway']['disease']
        if disease not in disease_groups:
            disease_groups[disease] = []
        disease_groups[disease].append((pathway_id, data))
    
    fig = go.Figure()
    
    # Color mapping for diseases
    colors = {
        'Diabetes Mellitus': '#FF6B6B',
        'Dyslipidemia': '#4ECDC4', 
        'Chronic Kidney Disease': '#45B7D1',
        'Metabolic Disorders': '#A78BFA',
        'Endocrine Disorders': '#F472B6',
        'Inflammatory Diseases': '#F59E0B',
        'Autoimmune Diseases': '#10B981',
        'Oxidative Stress': '#F97316',
        'Metabolic Syndrome': '#8B5CF6',
        'Urea Cycle Disorders': '#EC4899',
        'Phenylketonuria': '#14B8A6',
        'Porphyria': '#F43F5E',
        'Gout': '#0EA5E9'
    }
    
    # Position diseases in clusters
    disease_positions = {
        'Diabetes Mellitus': (0, 0),
        'Dyslipidemia': (3, 0),
        'Chronic Kidney Disease': (1.5, 2.5),
        'Metabolic Disorders': (-2, 1),
        'Endocrine Disorders': (5, 1),
        'Inflammatory Diseases': (1.5, -2),
        'Autoimmune Diseases': (4, -2),
        'Oxidative Stress': (-1, -1),
        'Metabolic Syndrome': (2, -1),
        'Urea Cycle Disorders': (-1, 2),
        'Phenylketonuria': (4, 1.5),
        'Porphyria': (-2, -1),
        'Gout': (3, -1.5)
    }
    
    # Add disease clusters
    for disease, pathways in disease_groups.items():
        # Use default position if not specified
        center_pos = disease_positions.get(disease, (0, 0))
        
        # Add disease label
        fig.add_trace(go.Scatter(
            x=[center_pos[0]], y=[center_pos[1] + 0.8],
            mode='text',
            text=disease,
            textfont=dict(size=16, color=colors.get(disease, '#333')),
            showlegend=False,
            hoverinfo='skip'
        ))
        
        # Add pathways in cluster
        for i, (pathway_id, data) in enumerate(pathways):
            angle = i * (2 * 3.14159 / len(pathways))
            radius = 0.5
            x = center_pos[0] + radius * math.cos(angle)
            y = center_pos[1] + radius * math.sin(angle)
            
            fig.add_trace(go.Scatter(
                x=[x], y=[y],
                mode='markers+text',
                marker=dict(
                    size=max(20, min(40, data['total_score'] * 10)),
                    color=colors.get(disease, '#333'),
                    opacity=0.7,
                    line=dict(width=2, color='white')
                ),
                text=data['pathway']['name'][:8] + '...' if len(data['pathway']['name']) > 8 else data['pathway']['name'],
                textposition="bottom center",
                textfont=dict(size=8),
                name=disease,
                hovertext=f"{data['pathway']['name']}<br>Score: {data['total_score']:.2f}<br>Click for detailed view",
                showlegend=i == 0
            ))
    
    # Add connections between related pathways
    for disease1, pathways1 in disease_groups.items():
        for disease2, pathways2 in disease_groups.items():
            if disease1 != disease2:
                # Get positions, default to (0,0) if not found
                pos1 = disease_positions.get(disease1, (0, 0))
                pos2 = disease_positions.get(disease2, (0, 0))
                
                # Add connection line between disease centers
                fig.add_trace(go.Scatter(
                    x=[pos1[0], pos2[0]], y=[pos1[1], pos2[1]],
                    mode='lines',
                    line=dict(width=1, color='lightgray', dash='dash'),
                    showlegend=False,
                    hoverinfo='skip'
                ))
    
    fig.update_layout(
        title="KEGG-Style Metabolic Pathway Network Overview",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-3, 6]),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-3, 4]),
        height=600,
        hovermode='closest'
    )
    
    return fig

def display_pathway_details(pathway_data, pathway_id):
    """Display detailed information about a specific pathway"""
    data = pathway_data[pathway_id]
    pathway = data['pathway']
    
    st.subheader(f"🔬 {pathway['name']}")
    st.write(f"**Associated Disease:** {pathway['disease']}")
    st.write(f"**Clinical Significance:** {pathway['clinical_significance']}")
    st.write(f"**KEGG Map:** [View Pathway]({pathway['kegg_map_url']})")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**Key Compounds (4):**")
        for i, (compound_id, compound_name) in enumerate(zip(pathway['compounds'][:4], pathway['compound_names'][:4])):
            st.write(f"• [{compound_name} ({compound_id})](https://www.kegg.jp/dbget-bin/www_bget?cpd:{compound_id})")
        
        st.write("**Enzymes (3):**")
        for i, (enzyme_id, enzyme_name) in enumerate(zip(pathway['enzymes'][:3], pathway['enzyme_names'][:3])):
            st.write(f"• [{enzyme_name} ({enzyme_id})](https://www.kegg.jp/dbget-bin/www_bget?ec:{enzyme_id.split(':')[1]})")
    
    with col2:
        st.write("**Reactions (3):**")
        for i, (reaction_id, reaction_name, rxn_type) in enumerate(zip(pathway['reactions'][:3], pathway['reaction_names'][:3], pathway['reaction_types'][:3])):
            st.write(f"• {reaction_name} ({reaction_id}) - {'⬅️➡️' if rxn_type == 'reversible' else '➡️'}")
        
        st.write("**Affected Organs:**")
        for organ in pathway['affected_organs']:
            st.write(f"• {organ}")
    
    st.write("**Relevant Biomarkers:**")
    biomarker_text = ", ".join(pathway['biomarkers'])
    st.write(biomarker_text)
    
    # Thermodynamic information
    favorability = "Highly Favorable" if pathway['deltaG'] < -50 else "Moderately Favorable" if pathway['deltaG'] < -20 else "Less Favorable"
    st.write(f"**Thermodynamic Status:** ΔG = {pathway['deltaG']:.1f} kJ/mol ({favorability})")

def generate_enhanced_report_content(lab_values, sorted_pathways, patient_info, ref_ranges):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Analyze directly affected and at-risk pathways
    directly_affected, at_risk = get_affected_pathways_analysis(lab_values, get_metabolic_database(), ref_ranges)
    
    lab_values_text = ""
    for param, value in lab_values.items():
        ref = ref_ranges.get(param, {})
        status = "Normal"
        if ref:
            if value < ref['min']:
                status = "Low ⬇️"
            elif value > ref['max']:
                status = "High ⬆️"
        
        lab_values_text += f"{param.upper()}: {value} {ref.get('unit', '')} (Ref: {ref.get('min', 'N/A')}-{ref.get('max', 'N/A')}) - {status}\n"
    
    # Directly affected pathways analysis
    directly_affected_text = ""
    if directly_affected:
        directly_affected_text = "DIRECTLY AFFECTED PATHWAYS:\n"
        for pathway in directly_affected:
            directly_affected_text += f"• {pathway['name']} ({pathway['disease']}) - Score: {pathway['score']:.2f}\n"
            if pathway['affected_biomarkers']:
                directly_affected_text += f"  Abnormal biomarkers: {', '.join(pathway['affected_biomarkers'])}\n"
    else:
        directly_affected_text = "DIRECTLY AFFECTED PATHWAYS: None detected\n"
    
    # At-risk pathways analysis
    at_risk_text = ""
    if at_risk:
        at_risk_text = "\nPATHWAYS AT RISK:\n"
        for pathway in at_risk:
            at_risk_text += f"• {pathway['name']} ({pathway['disease']}) - Score: {pathway['score']:.2f}\n"
            if pathway['risk_factors']:
                at_risk_text += f"  Risk factors: {', '.join(pathway['risk_factors'])}\n"
    else:
        at_risk_text = "\nPATHWAYS AT RISK: None detected\n"
    
    # Top pathway details
    top_pathway = sorted_pathways[0][1]['pathway']
    pathway_details = f"""
TOP PRIORITY PATHWAY DETAILS:
Name: {top_pathway['name']}
Disease Association: {top_pathway['disease']}
Key Compounds: {', '.join(top_pathway['compound_names'][:4])}
Key Enzymes: {', '.join(top_pathway['enzyme_names'][:3])}
Key Reactions: {', '.join(top_pathway['reaction_names'][:3])}
Thermodynamic Status: ΔG = {top_pathway['deltaG']:.1f} kJ/mol
KEGG Map: {top_pathway['kegg_map_url']}
"""
    
    abnormal_count = sum(1 for param, value in lab_values.items() 
                        if param in ref_ranges and 
                        (value < ref_ranges[param]['min'] or value > ref_ranges[param]['max']))
    
    report_content = f"""AIKEGG METABOLIC PATHWAY ANALYSIS REPORT
Generated: {timestamp}

PATIENT INFORMATION
Patient ID: {patient_info.get('patient_id', 'Not specified')}
Age: {patient_info.get('age', 'Not specified')}
Gender: {patient_info.get('gender', 'Not specified')}

LABORATORY VALUES
{lab_values_text}

PATHWAY PERTURBATION ANALYSIS
{directly_affected_text}
{at_risk_text}

{pathway_details}

CLINICAL RECOMMENDATIONS
• Immediate attention required for: {len(directly_affected)} pathway(s)
• Monitoring recommended for: {len(at_risk)} pathway(s)
• Consider follow-up testing for biomarkers in affected pathways

SUMMARY
Total Parameters Analyzed: {len(lab_values)}
Abnormal Values Detected: {abnormal_count}
Pathways Directly Affected: {len(directly_affected)}
Pathways At Risk: {len(at_risk)}
"""
    
    return report_content

def generate_pdf_report(lab_values, sorted_pathways, patient_info, ref_ranges):
    """Generate PDF report using ReportLab"""
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4)
    styles = getSampleStyleSheet()
    story = []
    
    # Title
    title_style = ParagraphStyle('CustomTitle', parent=styles['Heading1'], fontSize=18, spaceAfter=30)
    story.append(Paragraph("AIKEGG METABOLIC PATHWAY ANALYSIS REPORT", title_style))
    story.append(Spacer(1, 12))
    
    # Patient Info
    story.append(Paragraph("PATIENT INFORMATION", styles['Heading2']))
    patient_data = [
        ['Patient ID:', patient_info.get('patient_id', 'Not specified')],
        ['Age:', str(patient_info.get('age', 'Not specified'))],
        ['Gender:', patient_info.get('gender', 'Not specified')],
        ['Report Date:', datetime.now().strftime("%Y-%m-%d %H:%M:%S")]
    ]
    patient_table = Table(patient_data, colWidths=[2*inch, 3*inch])
    patient_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), colors.lightgrey),
        ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 0), (-1, -1), 10),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 12),
    ]))
    story.append(patient_table)
    story.append(Spacer(1, 20))
    
    # Lab Values
    story.append(Paragraph("LABORATORY VALUES", styles['Heading2']))
    lab_data = [['Parameter', 'Value', 'Reference Range', 'Status']]
    
    for param, value in lab_values.items():
        ref = ref_ranges.get(param, {})
        status = "Normal"
        if ref:
            if value < ref['min']:
                status = "Low"
            elif value > ref['max']:
                status = "High"
        
        lab_data.append([
            param.upper(),
            f"{value} {ref.get('unit', '')}",
            f"{ref.get('min', 'N/A')}-{ref.get('max', 'N/A')} {ref.get('unit', '')}",
            status
        ])
    
    lab_table = Table(lab_data, colWidths=[1.5*inch, 1.5*inch, 1.5*inch, 1*inch])
    lab_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black)
    ]))
    story.append(lab_table)
    story.append(Spacer(1, 20))
    
    # Top Pathways
    story.append(Paragraph("TOP AFFECTED PATHWAYS", styles['Heading2']))
    pathway_data = [['Rank', 'Pathway Name', 'Disease', 'Score']]
    
    for i, (pathway_id, data) in enumerate(sorted_pathways[:5], 1):
        pathway_data.append([
            str(i),
            data['pathway']['name'],
            data['pathway']['disease'],
            f"{data['total_score']:.2f}"
        ])
    
    pathway_table = Table(pathway_data, colWidths=[0.5*inch, 3*inch, 1.5*inch, 1*inch])
    pathway_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black)
    ]))
    story.append(pathway_table)
    
    doc.build(story)
    buffer.seek(0)
    return buffer

def create_pathway_table(sorted_pathways, directly_affected, at_risk):
    """Create an Excel-like table for pathway analysis"""
    
    # Prepare data for the table
    table_data = []
    
    for i, (pathway_id, data) in enumerate(sorted_pathways, 1):
        # Determine status
        status = "Normal"
        status_color = "🟢"
        
        # Check if pathway is directly affected
        for affected in directly_affected:
            if affected['pathway_id'] == pathway_id:
                status = "Directly Affected"
                status_color = "🔴"
                break
        
        # Check if pathway is at risk
        if status == "Normal":
            for risk in at_risk:
                if risk['pathway_id'] == pathway_id:
                    status = "At Risk"
                    status_color = "🟡"
                    break
        
        table_data.append({
            'Rank': i,
            'Status': f"{status_color} {status}",
            'Pathway Name': data['pathway']['name'],
            'Disease Association': data['pathway']['disease'],
            'Perturbation Score': f"{data['perturbation_score']:.2f}",
            'Thermodynamic Score': f"{data['thermodynamic_score']:.2f}",
            'Total Score': f"{data['total_score']:.2f}",
            'Affected Organs': ', '.join(data['pathway']['affected_organs'][:3]),
            'Key Biomarkers': ', '.join(data['pathway']['biomarkers'][:3]),
            'ΔG (kJ/mol)': f"{data['pathway']['deltaG']:.1f}"
        })
    
    return table_data

def main():
    st.title("🧬 AIKEGG - AI-Enhanced Metabolic Pathway Analysis")
    st.markdown("Advanced AI platform that maps patient biochemical data to KEGG metabolic pathways for precision diagnostics")
    
    metabolic_db = get_metabolic_database()
    ref_ranges = get_reference_ranges()
    
    col1, col2 = st.columns([1, 2])
    
    with col1:
        st.header("Patient Information")
        patient_id = st.text_input("Patient ID", placeholder="Enter patient identifier")
        age = st.number_input("Age", min_value=1, max_value=120, value=None)
        gender = st.selectbox("Gender", ["", "Male", "Female", "Other"])
        
        st.header("Disease Modules")
        st.info("Comprehensive metabolic pathway analysis across multiple disease domains")
        
        # Expanded panels with all available biomarkers
        with st.expander("🩸 Metabolic Panel", expanded=True):
            glucose = st.number_input("Glucose (mg/dL)", min_value=0.0, value=None, step=0.1)
            lactate = st.number_input("Lactate (mmol/L)", min_value=0.0, value=None, step=0.01)
            pyruvate = st.number_input("Pyruvate (mmol/L)", min_value=0.0, value=None, step=0.01)
            urea = st.number_input("Urea (mg/dL)", min_value=0.0, value=None, step=0.1)
            creatinine = st.number_input("Creatinine (mg/dL)", min_value=0.0, value=None, step=0.01)
            bun = st.number_input("BUN (mg/dL)", min_value=0.0, value=None, step=0.1)
            ammonia = st.number_input("Ammonia (μg/dL)", min_value=0.0, value=None, step=0.1)
            alt = st.number_input("ALT (U/L)", min_value=0.0, value=None, step=0.1)
            ast = st.number_input("AST (U/L)", min_value=0.0, value=None, step=0.1)
            bilirubin = st.number_input("Bilirubin (mg/dL)", min_value=0.0, value=None, step=0.01)
        
        with st.expander("🧪 Endocrine Panel"):
            vitamin_d = st.number_input("Vitamin D (ng/mL)", min_value=0.0, value=None, step=0.1)
            pth = st.number_input("PTH (pg/mL)", min_value=0.0, value=None, step=0.1)
            tsh = st.number_input("TSH (mIU/L)", min_value=0.0, value=None, step=0.01)
            cortisol = st.number_input("Cortisol (μg/dL)", min_value=0.0, value=None, step=0.1)
            insulin = st.number_input("Insulin (μIU/mL)", min_value=0.0, value=None, step=0.1)
            c_peptide = st.number_input("C-peptide (ng/mL)", min_value=0.0, value=None, step=0.01)
            estrogen = st.number_input("Estrogen (pg/mL)", min_value=0.0, value=None, step=1.0)
            progesterone = st.number_input("Progesterone (ng/mL)", min_value=0.0, value=None, step=0.1)
            testosterone = st.number_input("Testosterone (ng/mL)", min_value=0.0, value=None, step=0.1)
        
        with st.expander("🧴 Lipid Panel"):
            apo_a = st.number_input("Apo A (mg/dL)", min_value=0.0, value=None, step=0.1)
            apo_b = st.number_input("Apo B (mg/dL)", min_value=0.0, value=None, step=0.1)
            total_cholesterol = st.number_input("Total Cholesterol (mg/dL)", min_value=0.0, value=None, step=0.1)
            ldl = st.number_input("LDL (mg/dL)", min_value=0.0, value=None, step=0.1)
            hdl = st.number_input("HDL (mg/dL)", min_value=0.0, value=None, step=0.1)
            triglycerides = st.number_input("Triglycerides (mg/dL)", min_value=0.0, value=None, step=0.1)
            lipoprotein_a = st.number_input("Lipoprotein(a) (mg/dL)", min_value=0.0, value=None, step=0.1)
        
        with st.expander("⚡ Minerals & Electrolytes"):
            phosphorus = st.number_input("Phosphorus (mg/dL)", min_value=0.0, value=None, step=0.1)
            magnesium = st.number_input("Magnesium (mg/dL)", min_value=0.0, value=None, step=0.01)
            uric_acid = st.number_input("Uric Acid (mg/dL)", min_value=0.0, value=None, step=0.1)
        
        with st.expander("🛡️ Immune Markers"):
            anti_tpo = st.number_input("Anti-TPO (IU/mL)", min_value=0.0, value=None, step=0.1)
            ana = st.number_input("ANA (titer)", min_value=0.0, value=None, step=1.0)
            anti_ccp = st.number_input("Anti-CCP (U/mL)", min_value=0.0, value=None, step=1.0)
            crp = st.number_input("CRP (mg/L)", min_value=0.0, value=None, step=0.1)
            fgf_23 = st.number_input("FGF-23 (RU/mL)", min_value=0.0, value=None, step=1.0)
        
        with st.expander("🧬 Other Specialized Tests"):
            homocysteine = st.number_input("Homocysteine (μmol/L)", min_value=0.0, value=None, step=0.1)
            glutamine = st.number_input("Glutamine (μmol/L)", min_value=0.0, value=None, step=1.0)
            glutamate = st.number_input("Glutamate (μmol/L)", min_value=0.0, value=None, step=0.1)
            gsh = st.number_input("GSH (μmol/L)", min_value=0.0, value=None, step=0.1)
            mda = st.number_input("MDA (μmol/L)", min_value=0.0, value=None, step=0.01)
        
        if st.button("🔬 Perform AIKEGG Analysis", type="primary"):
            lab_values = {}
            
            # Collect all non-None values from all panels
            values_map = {
                'glucose': glucose, 'lactate': lactate, 'pyruvate': pyruvate,
                'urea': urea, 'creatinine': creatinine, 'bun': bun,
                'ammonia': ammonia, 'alt': alt, 'ast': ast, 'bilirubin': bilirubin,
                'vitamin_d': vitamin_d, 'pth': pth, 'tsh': tsh, 'cortisol': cortisol,
                'insulin': insulin, 'c_peptide': c_peptide, 'estrogen': estrogen,
                'progesterone': progesterone, 'testosterone': testosterone,
                'apo_a': apo_a, 'apo_b': apo_b, 'total_cholesterol': total_cholesterol,
                'ldl': ldl, 'hdl': hdl, 'triglycerides': triglycerides,
                'phosphorus': phosphorus, 'magnesium': magnesium,
                'uric_acid': uric_acid, 
                'anti_tpo': anti_tpo, 'ana': ana, 'anti_ccp': anti_ccp,
                'crp': crp, 'fgf_23': fgf_23, 'homocysteine': homocysteine,
                'glutamine': glutamine, 'glutamate': glutamate, 'gsh': gsh,
                'mda': mda
            }
            
            for key, value in values_map.items():
                if value is not None:
                    lab_values[key] = value
            
            if not lab_values:
                st.error("Please enter at least one lab value to analyze.")
                return
            
            pathway_scores = get_pathway_perturbation_score(lab_values, metabolic_db, ref_ranges)
            sorted_pathways = sorted(pathway_scores.items(), key=lambda x: x[1]['total_score'], reverse=True)
            
            # Get affected pathways analysis
            directly_affected, at_risk = get_affected_pathways_analysis(lab_values, metabolic_db, ref_ranges)
            
            st.session_state['analysis_results'] = {
                'lab_values': lab_values,
                'sorted_pathways': sorted_pathways,
                'directly_affected': directly_affected,
                'at_risk': at_risk,
                'patient_info': {
                    'patient_id': patient_id,
                    'age': age,
                    'gender': gender
                }
            }
    
    with col2:
        if 'analysis_results' in st.session_state:
            results = st.session_state['analysis_results']
            lab_values = results['lab_values']
            sorted_pathways = results['sorted_pathways']
            directly_affected = results['directly_affected']
            at_risk = results['at_risk']
            patient_info = results['patient_info']
            
            tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
                "Overview", "Affected Pathways", "Pathway Details", 
                "Network View", "Thermodynamics & Flux", "Enhanced Report",
                "Patient Database"
            ])
            
            with tab1:
                st.subheader("🎯 Analysis Overview")
                
                # Key metrics
                col_a, col_b, col_c, col_d = st.columns(4)
                with col_a:
                    st.metric("Parameters", len(lab_values))
                with col_b:
                    st.metric("Directly Affected", len(directly_affected), delta=f"{len(directly_affected)} pathways")
                with col_c:
                    st.metric("At Risk", len(at_risk), delta=f"{len(at_risk)} pathways")
                with col_d:
                    st.metric("Top Score", f"{sorted_pathways[0][1]['total_score']:.2f}")
                
                # Overview chart
                st.plotly_chart(create_overview_chart(sorted_pathways), use_container_width=True)
                
                # Quick summary
                st.subheader("🚨 Immediate Attention Required")
                if directly_affected:
                    for pathway in directly_affected[:3]:
                        st.error(f"**{pathway['name']}** ({pathway['disease']}) - Score: {pathway['score']:.2f}")
                else:
                    st.success("No pathways require immediate attention")
            
            with tab2:
                st.subheader("📊 Pathway Perturbation Analysis")
                # Create Excel-like table
                table_data = create_pathway_table(sorted_pathways, directly_affected, at_risk)
                
                df = pd.DataFrame(table_data)

                def highlight_status(val):
                    if '🔴' in val:
                        return 'background-color: #ffebee; color: #c62828; font-weight: bold'
                    elif '🟡' in val:
                        return 'background-color: #fff3e0; color: #ef6c00; font-weight: bold'
                    elif '🟢' in val:
                        return 'background-color: #e8f5e8; color: #2e7d32; font-weight: bold'
                    return ''
                
                def highlight_scores(val):
                    try:
                        score = float(val)
                        if score > 2.0:
                            return 'background-color: #ffcdd2; font-weight: bold'
                        elif score > 1.0:
                            return 'background-color: #ffe0b2; font-weight: bold'
                        return ''
                    except:
                        return ''
                
                styled_df = df.style.applymap(highlight_status, subset=['Status']) \
                                    .applymap(highlight_scores, subset=['Total Score']) \
                                    .set_properties(**{
                                        'font-size': '12px',
                                        'border': '1px solid #ddd'
                                    })
                st.dataframe(styled_df, use_container_width=True, height=400)
                
                # Summary statistics
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Total Pathways", len(sorted_pathways))
                with col2:
                    st.metric("Directly Affected", len(directly_affected))
                with col3:
                    st.metric("At Risk", len(at_risk))
                with col4:
                    avg_score = sum(data[1]['total_score'] for data in sorted_pathways) / len(sorted_pathways)
                    st.metric("Average Score", f"{avg_score:.2f}")    

                if directly_affected:
                    st.error("🔴 **DIRECTLY AFFECTED PATHWAYS**")
                    for pathway in directly_affected:
                        with st.expander(f"{pathway['name']} ({pathway['disease']}) - Score: {pathway['score']:.2f}"):
                            st.write("**Abnormal Biomarkers:**")
                            if pathway['affected_biomarkers']:
                                for biomarker in pathway['affected_biomarkers']:
                                    st.write(f"• {biomarker}")
                            else:
                                st.write("• Based on clinical correlation and pathway relevance")

                if at_risk:
                    st.warning("🟡 **PATHWAYS AT RISK**")
                    for pathway in at_risk:
                        with st.expander(f"{pathway['name']} ({pathway['disease']}) - Score: {pathway['score']:.2f}"):
                            st.write("**Risk Factors:**")
                            if pathway['risk_factors']:
                                for factor in pathway['risk_factors']:
                                    st.write(f"• {factor}")
                            else:
                                st.write("• Moderate perturbation detected")                                        
            
            with tab3:
                st.subheader("🔬 Detailed Pathway Information")
                
                # Select pathway for detailed view
                pathway_options = {f"{data[1]['pathway']['name']} ({data[1]['pathway']['disease']})": data[0] 
                                 for data in sorted_pathways}
                selected_pathway_name = st.selectbox("Select pathway for detailed analysis:", 
                                                   list(pathway_options.keys()))
                
                if selected_pathway_name:
                    selected_pathway_id = pathway_options[selected_pathway_name]
                    pathway_data = {pid: data for pid, data in sorted_pathways}
                    display_pathway_details(pathway_data, selected_pathway_id)
            
            with tab4:
                st.subheader("🕸️ KEGG-Style Metabolic Pathway Network")
                
                # Add pathway selector for detailed view
                col_select, col_button = st.columns([3, 1])
                
                with col_select:
                    pathway_options = ["Overview"] + [f"{data[1]['pathway']['name']} ({data[1]['pathway']['disease']})" 
                                                for data in sorted_pathways]
                    selected_view = st.selectbox("Select view:", pathway_options)
                
                with col_button:
                    st.write("")  # Spacing
                    show_overview = st.button("🔄 Show Overview")
                
                pathway_data = {pid: data for pid, data in sorted_pathways}
                
                if show_overview or selected_view == "Overview":
                    network_fig = create_kegg_style_pathway_chart(pathway_data)
                    st.plotly_chart(network_fig, use_container_width=True)
                    st.info("🔵 Blue circles = Compounds (click for KEGG) | 🟩 Green squares = Enzymes (click for KEGG) | ➡️ Red arrows = Reactions")
                else:
                    # Find selected pathway ID
                    selected_pathway_id = None
                    for pid, data in sorted_pathways:
                        pathway_name = f"{data['pathway']['name']} ({data['pathway']['disease']})"
                        if pathway_name == selected_view:
                            selected_pathway_id = pid
                            break
                    
                    if selected_pathway_id:
                        detailed_fig = create_kegg_style_pathway_chart(pathway_data, selected_pathway_id)
                        st.plotly_chart(detailed_fig, use_container_width=True)
                        
                        # Show additional pathway details
                        pathway_info = pathway_data[selected_pathway_id]['pathway']
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.write("**Key Metabolites:**")
                            for comp_name in pathway_info['compound_names'][:4]:
                                st.write(f"• {comp_name}")
                        
                        with col2:
                            st.write("**Key Enzymes:**")
                            for enz_name in pathway_info['enzyme_names'][:3]:
                                st.write(f"• {enz_name}")
                        
                        st.write("**Biochemical Reactions:**")
                        for reaction_name, rxn_type in zip(pathway_info['reaction_names'][:3], pathway_info['reaction_types'][:3]):
                            st.write(f"➡️ {reaction_name} ({'⬅️➡️' if rxn_type == 'reversible' else '➡️'})")
            
            with tab5:
                st.info("This analysis shows the thermodynamic favorability and metabolic flux perturbations")
                
                # Thermodynamic Analysis
                st.subheader("🔥 Thermodynamic Favorability")
                col_thermo, col_info = st.columns([2, 1])
                
                with col_thermo:
                    st.plotly_chart(create_thermodynamics_chart(sorted_pathways), use_container_width=True)
                
                with col_info:
                    st.write("**Interpretation:**")
                    st.write("• ΔG < -50: Highly Favorable")
                    st.write("• ΔG -50 to -20: Moderately Favorable")
                    st.write("• ΔG > -20: Less Favorable")
                    
                    # Show thermodynamic summary table
                    thermo_data = []
                    for pathway_id, data in sorted_pathways[:5]:
                        favorability = "Highly Favorable" if data['pathway']['deltaG'] < -50 else "Moderately Favorable" if data['pathway']['deltaG'] < -20 else "Less Favorable"
                        thermo_data.append({
                            'Pathway': data['pathway']['name'][:20] + '...' if len(data['pathway']['name']) > 20 else data['pathway']['name'],
                            'ΔG (kJ/mol)': data['pathway']['deltaG'],
                            'Favorability': favorability
                        })
                    
                    st.write("**Top 5 Pathways:**")
                    st.dataframe(pd.DataFrame(thermo_data), use_container_width=True)
                
                st.divider()
                
                # Flux Analysis
                st.subheader("📊 Metabolic Flux Perturbation")
                col_flux, col_flux_info = st.columns([2, 1])
                
                with col_flux:
                    st.plotly_chart(create_flux_chart(sorted_pathways), use_container_width=True)
                
                with col_flux_info:
                    st.write("**Flux Analysis:**")
                    st.write("• Red area: Perturbed flux state")
                    st.write("• Green area: Normal flux state")
                    st.write("• Larger deviation = Higher perturbation")
                    
                    # Flux perturbation summary
                    flux_data = []
                    for pathway_id, data in sorted_pathways[:5]:
                        perturbation_percent = data['perturbation_score'] * 30  # Convert to percentage
                        flux_data.append({
                            'Pathway': data['pathway']['name'][:20] + '...' if len(data['pathway']['name']) > 20 else data['pathway']['name'],
                            'Perturbation %': f"{perturbation_percent:.1f}%",
                            'Status': "High" if perturbation_percent > 45 else "Moderate" if perturbation_percent > 15 else "Low"
                        })
                    
                    st.write("**Flux Perturbations:**")
                    st.dataframe(pd.DataFrame(flux_data), use_container_width=True)
            
            with tab6:
                st.subheader("📋 Enhanced Clinical Report")
                
                enhanced_report = generate_enhanced_report_content(lab_values, sorted_pathways, patient_info, ref_ranges)
                
                st.text_area("Clinical Report", enhanced_report, height=600)
                
                col_download1, col_download2, col_actions = st.columns(3)
                
                with col_download1:
                    st.download_button(
                        label="📄 Download TXT Report",
                        data=enhanced_report,
                        file_name=f"AIKEGG_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                        mime="text/plain"
                    )
                
                with col_download2:
                    try:
                        pdf_buffer = generate_pdf_report(lab_values, sorted_pathways, patient_info, ref_ranges)
                        st.download_button(
                            label="📑 Download PDF Report",
                            data=pdf_buffer,
                            file_name=f"AIKEGG_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                            mime="application/pdf"
                        )
                    except Exception as e:
                        st.error(f"PDF generation failed. Install reportlab: pip install reportlab")
                
                with col_actions:
                    # Save to CSV button
                    if st.button("💾 Save to Database", help="Save this patient's data to the permanent database"):
                        save_to_csv(csv_file, patient_info, lab_values, results)
                        st.success("Patient data saved to database!")
            
            with tab7:
                st.subheader("📊 Patient Database")
                st.info("View and search all historical patient records")
                
                # Load and display CSV data
                try:
                    df = load_csv_data(csv_file)
                    
                    if df.empty:
                        st.warning("No patient data found in the database.")
                    else:
                        # Search functionality
                        col_search1, col_search2 = st.columns(2)
                        with col_search1:
                            search_term = st.text_input("Search by Patient ID", "")
                        with col_search2:
                            date_filter = st.date_input("Filter by date range", [])
                        
                        # Apply filters
                        if search_term:
                            df = df[df['patient_id'].str.contains(search_term, case=False, na=False)]
                        
                        if date_filter:
                            if isinstance(date_filter, list) and len(date_filter) == 2:
                                df['timestamp'] = pd.to_datetime(df['timestamp'])
                                df = df[(df['timestamp'].dt.date >= date_filter[0]) & 
                                       (df['timestamp'].dt.date <= date_filter[1])]
                        
                        # Display data
                        st.dataframe(df, use_container_width=True, height=600)
                        
                        # Download options
                        csv_data = df.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            label="📥 Download CSV",
                            data=csv_data,
                            file_name="AIKEGG_patient_database.csv",
                            mime="text/csv"
                        )
                        
                        # Show basic stats
                        st.write(f"Total records: {len(df)}")
                        if not df.empty:
                            st.write(f"Last entry: {df['timestamp'].iloc[0].strftime('%Y-%m-%d %H:%M:%S')}")
                
                except Exception as e:
                    st.error(f"Error loading database: {str(e)}")

if __name__ == "__main__":
    main()
