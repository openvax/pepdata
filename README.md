epitopes
=======

Computational immunology seeks to investigate and model the properties of *peptides* (short strings of amino acids). Peptides can be substrings of some larger protein, naturally occurring small proteins, or synthesized for therapeutic purposes. We can investigate whether a peptide is presented on the cell surface for immune recognition by [MHC molecules](http://en.wikipedia.org/wiki/Major_histocompatibility_complex), in which case we call it an *epitope*. If an epitope elicits a response from [white blood cells](http://en.wikipedia.org/wiki/Lymphocyte) it can be called  *immunogenic*. To make useful predictions (i.e. "what peptides should go in this vaccine?") we need to partition the combinatorial space of peptides into classes such as epitopes vs. non-epitope and immunogenic vs. non-immunogenic 
Computational immunology. One approach to modeling such distinctions is to collect a large volume of data about peptides to model their properties statistically. This library seeks to aid you in that process by providing simple Python/NumPy/Pandas interfaces to commonly used immunology and bioninformatics datasets. 

**Data Sources** 

- `iedb`: [Immune Epitope Database](http://www.iedb.org), large collection of epitope assay results for MHC binding as well as T-cell/B-cell responses
- `tcga`: Variant peptide substrings extracted from [TCGA](http://en.wikipedia.org/wiki/The_Cancer_Genome_Atlas) mutations across all cancer types
- `reference`: Peptide substrings from the [human reference protein sequence](ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/pep/)
- `imma2`: IMMA2 epitope immunogenic vs. non-immunogenic data set used by Tung et al. for evaluating the [POPISK](http://www.biomedcentral.com/1471-2105/12/446) immunogenicity predictor 
- `calis`: Two datasets used in Calis et al.'s [Properties of MHC Class I Presented Peptides That Enhance Immunogenicity](http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003266#pcbi.1003266.s005) 
- `hpv`: [Human Papillomavirus T cell Antigen Database](http://cvc.dfci.harvard.edu/cvccgi/hpv/)
- `toxin`: Toxic protein sequences from [Animal Toxin Databse](http://protchem.hunnu.edu.cn/toxin/)
- `danafarber`: [Dana Farber Repository for Machine Learning in Immunology](http://bio.dfci.harvard.edu/DFRMLI/)
- `cri_tumor_antigens`: Tumor associated peptides from [Cancer Immunity](http://cancerimmunity.org/peptide/mutations/)

Planned:

- `tantigen`: [Tumor T-cell Antigen Database](http://cvc.dfci.harvard.edu/tadb/)
- `bcipep`: [B-cell epitopes](http://www.imtech.res.in/raghava/bcipep/data.html) 


**API**

When a dataset consists only of an unlabeled list of epitopes, then it only needs two functions:
- `load_wuzzle`: Returns set of amino acid strings 
- `load_wuzzle_ngrams`: Array whose rows are amino acids transformed into n-gram vector space. 

If the dataset contains additional data about the epitopes (such as HLA type u or source protein):
- `load_wuzzle`: Returns data frame with epitope strings and additional properties
- `load_wuzzle_set`: Set of epitope amino acid strings 
- `load_wuzzle_ngrams`: Array whose rows are amino acids transformed into n-gram vector space. 

If the dataset is labeled (contains positive and negative assay results), then the following functions should be available: 
- `load_wuzzle`: Load all available data from the "wuzzle" dataset (filtered by options such as `mhc_class`). 
- `load_wuzzle_values`: Group the dataset by epitope string and associate each epitope with the percentage of positive results. 
- `load_wuzzle_classes`: Split the epitopes into positive and negative classes, return a set of strings for each. 
- `load_wuzzle_ngrams`: Transform the amino acid string representation (or some reduced alphabet) into vectors of n-gram frequencies, return a sklearn-compatible `(samples, labels)` pair of arrays.   
