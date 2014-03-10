epitopes
=======

Python interface to IEDB and other immunology datasets (MHC binding and T-cell response assays)

**Data Sources** 

- `iedb`: [Immune Epitope Database](http://www.iedb.org), large collection of epitope assay results for MHC binding as well as T-cell/B-cell responses
- `imma2`: IMMA2 epitope immunogenic vs. non-immunogenic data set used by Tung et al. for evaluating the [POPISK](http://www.biomedcentral.com/1471-2105/12/446) immunogenicity predictor 
- `calis`: Two datasets used in Calis et al.'s [Properties of MHC Class I Presented Peptides That Enhance Immunogenicity](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003266#pcbi.1003266.s005) 
- `tantigen`: [Tumor T-cell Antigen Database](http://cvc.dfci.harvard.edu/tadb/)
- `bcipep`: [B-cell epitopes](http://www.imtech.res.in/raghava/bcipep/data.html) 
- `reference`: Peptide substrings of a reference exome
- `cancer_immunity`: Tumor associated peptides from [Cancer Immunity](http://cancerimmunity.org/peptide/mutations/)
- `tcga`: Variant peptide substrings extracted from TCGA mutations 
- `hpv`: [Human Papillomavirus T cell Antigen Database](http://cvc.dfci.harvard.edu/cvccgi/hpv/)
- `toxin`: Toxic protein sequences from [Animal Toxin Databse](http://protchem.hunnu.edu.cn/toxin/)
- 
**Installing**

Eventually we'll have a proper process for downloading data to a user directory via Python. For now, install using `python setup.py develop` and then run `get_data.sh`. 
