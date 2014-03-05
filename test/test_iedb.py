from epitope import iedb

def test_tcell_hla_a24():
  """
  Test that HLA restriction actually decreases number of results and
  that regular expression patterns are being used correctly
  """
  df_all = iedb.load_tcell()
  df_a24_1 = iedb.load_tcell(hla_type='HLA-A24')
  df_a24_2 = iedb.load_tcell(hla_type='HLA-A\*24')
  df_a24_combined = iedb.load_tcell(hla_type = 'HLA-A24|HLA-A\*24')
  assert len(df_a24_1) < len(df_all)
  assert len(df_a24_2) < len(df_all)
  assert len(df_a24_combined) <= len(df_a24_1) + len(df_a24_2), \
   "Expected %d <= %d + %d" % (len(df_a24_combined), len(df_a24_1), len(df_a24_2))


def test_mhc_hla_a2():
  """
  Test that HLA restriction actually decreases number of results and
  that regular expression patterns are being used correctly
  """
  df_all = iedb.load_mhc()
  df_a2_1 = iedb.load_mhc(hla_type='HLA-A2')
  df_a2_2 = iedb.load_mhc(hla_type='HLA-A\*02')
  df_a2_combined = iedb.load_mhc(hla_type = 'HLA-A2|HLA-A\*02')
  assert len(df_a2_1) < len(df_all)
  assert len(df_a2_2) < len(df_all)
  assert len(df_a2_combined) <= len(df_a2_1) + len(df_a2_2), \
    "Expected %d <= %d + %d" % (len(df_a2_combined), len(df_a2_1), len(df_a2_2))
