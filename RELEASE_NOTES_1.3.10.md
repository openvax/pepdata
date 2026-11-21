# PepData 1.3.10

Maintenance and hygiene release:

- Fix `pepdata.chou_fasman` (`alpha_helix_score`, `beta_sheet_score`, `turn_score`) to key residue dictionaries by one-letter amino acid codes (`"A"`–`"V"`) instead of integer alphabet indices, matching every other table in `pepdata`.
- Upgrade remaining `http://` IEDB download URLs in `pepdata.iedb.alleles` and `pepdata.iedb.tcell` to `https://`.
- Remove unused `mhc_class_series` lookup in `pepdata.iedb.tcell.load_dataframe`, fix the `df: pd.DataFrame` annotation on `pepdata.iedb.columns.get_mhc_assay`, replace the lambda preprocessor assignment in `pepdata.peptide_vectorizer`, and tighten `[tool.ruff.lint]` ignores to catch `F821`, `F841`, `E731`, `E722`, and `E721`.
