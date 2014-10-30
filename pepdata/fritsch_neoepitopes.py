import pandas as pd
from os.path import join
from static_data import DATA_DIR


"""
Data from 'HLA-binding properties of tumor neoepitopes in humans' 
by Edward F Fritsch, Mohini Rajasagi, Patrick A Ott, et al. 
Cancer Immunology Research
"""

def load_dataframe(
        hla_type = None,
        exclude_hla_type = None):
    path = join(DATA_DIR, 'fritsch2014_neoepitopes.csv')
    df = pd.read_csv(path, skipinitialspace=True)
    hla = df['HLA Allele']
    if hla_type:
        df = df[hla.str.contains(hla_type, na=False).astype('bool')]
    if exclude_hla_type:
        df = df[~(hla.str.contains(exclude_hla_type, na=True).astype('bool'))]
    return df