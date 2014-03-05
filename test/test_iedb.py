import iedb 

def test_hla_a24():
  """
  Test that HLA restriction actually decreases number of results and 
  that regular expression patterns are being used correctly
  """
  df_all = iedb.load_tcell_prct()
  df_a24_1 = iedb.load_tcell_prct(hla_type='HLA-A24')
  df_a24_2 = iedb.load_tcell_prct(hla_type='HLA-A\*24')
  df_a24_combined = iedb.load_tcell_prct(hla_type = 'HLA-A24|HLA-A\*24')
  assert len(df_a24_1) < len(df_all)
  assert len(df_a24_2) < len(df_all)
  assert len(df_combined) == len(df_a24_1) + len(df_a24_2)

