from download import fetch_data

"""
T-cell epitopes from Dana-Farber Repository for Machine Learning in Immunology
http://bio.dfci.harvard.edu/DFRMLI/HTML/TCellEpitopes.php
"""


def load_df_tumor():
    """
    Tumor antigens.
    This data set is a list of 718 T cell epitopes, 8-31 amino acids in length,
    derived from human tumor antigens. They are restricted by multiple HLA class
    I and class II alleles. The epitopes were experimentally validated based on
    immune recognition of HLA molecules and T cell stimulation. It could be used
    as a validation dataset for computational models to predict T cell epitopes.
    """
    path = fetch_data(
        filename = "tumor_epitopes.csv",
        download_url = \
            "http://bio.dfci.harvard.edu/DFRMLI/datasets/tumor_epitopes.htm")

def load_df_virus():
    """
    Virus antigens.
    This data set is a list of 44 HLA-A2 restricted T cell epitopes,
    9 amino acids in length, derived from human medically important viruses.
    The epitopes were experimentally characterized based on immune recognition
    of HLA molecule and T cell stimulation. It could be used as a validation
    dataset for computational models to predict T cell epitopes.
    """
    path = fetch_data(
        filename = "virus_epitopes_A2.csv",
        download_url = \
            "http://bio.dfci.harvard.edu/DFRMLI/datasets/virus_epitopes_A2.htm")

def load_df_cef():
    """
    This dataset is a list of 32 T cell epitopes, 8-12 amino acids in length,
    with sequences derived from the human Cytomegalovirus, Epstein-Barr Virus
    and Influenza Virus.
    """
    path = fetch_data(
        filename = "CEF.csv",
        download_url = "http://bio.dfci.harvard.edu/DFRMLI/datasets/CEF.htm")
