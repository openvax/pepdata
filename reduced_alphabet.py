
"""
Amino acid groupings from 
'Reduced amino acid alphabets improve the sensitivity...' by 
Peterson, Kondev, et al. 
http://www.rpgroup.caltech.edu/publications/Peterson2008.pdf
"""


def dict_from_list(groups):
  result = {}
  for i, group in enumerate(groups):
    for c in group:
      result[c.upper()] = i
      result[c.lower()] = i
  return result 
  

gbmr4 = dict_from_list(["ADKERNTSQ", "YFLIVMCWH", "G", "P"])

sdm12 = dict_from_list(
  ["A", "D", "KER", "N",  "TSQ", "YF", "LIVM", "C", "W", "H", "G", "P"]
)

hsdm17 = dict_from_list(
  ["A", "D", "KE", "R", "N", "T", "S", "Q", "Y", "F", "LIV", "M", "C", "W", "H", "G", "P"])

"""
Other alphabets from 
http://bio.math-inf.uni-greifswald.de/viscose/html/alphabets.html
"""

#hydrophilic vs. hydrophobic
hp2 = dict_from_list(["AGTSNQDEHRKP", "CMFILVWY"])

murphy10 = dict_from_list(
  ["LVIM", "C", "A", "G", "ST", "P", "FYW", "EDNQ", "KR", "H"]
)


alex6 = dict_from_list(["C", "G", "P", "FYW", "AVILM", "STNQRHKDE"])

aromatic2 = dict_from_list(["FHWY", "ADKERNTSQLIVMCGP"])