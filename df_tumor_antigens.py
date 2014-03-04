#!/usr/bin/python


import pandas as pd


if __name__ == '__main__':
    data = pd.read_html("http://bio.dfci.harvard.edu/DFRMLI/datasets/tumor_epitopes.htm", header=0,infer_types=False)[0]
    data.to_csv('danafarber_verified_antigens.txt', sep='\t', index=False, encoding='utf-8')
