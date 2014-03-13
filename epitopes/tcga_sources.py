# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Bladder Urothelial Carcinoma
BLCA_URL = \
"https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/blca/gsc/broad.mit.edu/illuminaga_dnaseq/mutations/broad.mit.edu_BLCA.IlluminaGA_DNASeq.Level_2.1.3.0/BLCA-28-original.aggregated.tcga.somatic.maf"

# Pancreatic Adenocarcinoma
PAAD_URL = \
"https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/paad/gsc/broad.mit.edu/illuminaga_dnaseq/mutations/broad.mit.edu_PAAD.IlluminaGA_DNASeq.Level_2.0.3.0/PR_TCGA_PAAD_PAIR_Capture_All_Pairs_QCPASS_v3.aggregated.capture.tcga.uuid.somatic.maf"

# Head and Neck Squamous Carcinoma
HNSC_URL = \
"https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/hnsc/gsc/broad.mit.edu/illuminaga_dnaseq/mutations/broad.mit.edu_HNSC.IlluminaGA_DNASeq.Level_2.1.0.0/PR_TCGA_HNSC_PAIR_Capture_TP-NT_TP-NB.aggregated.capture.tcga.uuid.somatic.maf"

# Glioblastoma
GBM_URL = \
"https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/gbm/gsc/broad.mit.edu/illuminaga_dnaseq/mutations/broad.mit.edu_GBM.IlluminaGA_DNASeq.Level_2.100.1.0/step4_gbm_liftover.aggregated.capture.tcga.uuid.maf2.4.migrated.somatic.maf"

# map cancer types to MAF download URLs
TCGA_SOURCES = {
    'blca' : BLCA_URL,
    'paad' : PAAD_URL,
    'hnsc' : HNSC_URL,
    'gbm' : GBM_URL
}

REFSEQ_PROTEIN_URL = \
'ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz'
