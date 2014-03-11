epitopes
=======

Python interface to IEDB and other immunology datasets (MHC binding and T-cell response assays)

**API**

When a dataset consists only of an unlabeled list of epitopes, then it only needs two functions:
- `load_wuzzle`: Returns set of amino acid strings 
- `load_wuzzle_ngrams`: Array whose rows are amino acids transformed into n-gram vector space. 

If, however, the dataset consists of assay results, then the following functions are available: 
- `load_wuzzle`: Load all available data from the "wuzzle" dataset (filtered by options such as `mhc_class`). 
- `load_wuzzle_values`: Group the dataset by epitope string and associate each epitope with the percentage of positive results. 
- `load_wuzzle_classes`: Split the epitopes into positive and negative classes, return a set of strings for each. 
- `load_wuzzle_ngrams`: Transform the amino acid string representation (or some reduced alphabet) into vectors of n-gram frequencies, return a sklearn-compatible `(samples, labels)` pair of arrays.   

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
- `df`: Dana Farber Repository for Machine Learning in Immunology
