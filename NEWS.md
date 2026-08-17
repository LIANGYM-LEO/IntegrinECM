# IntegrinECM 0.1.2

- Added file-path and MD5 input auditing plus an explicit completed run status.
- Added `combine_integrin_ecm_results()` with sample, metadata, parameter, and
  duplicate-module compatibility checks.
- Added `write_reference_workbook()` for the 284-record interaction resource
  and all fixed structure catalogs and component stoichiometry.
- Expanded `write_score_workbook()` with individually selectable technical
  sheets and the `intermediates`, `references`, and `all` aliases.
- Added export of gene coverage, composition matrices, geometric-mean `u` and
  `q` intermediates, selected catalog components, and detailed input audits.

# IntegrinECM 0.1.1

- Made workbook export portable across macOS, Windows, and Linux by writing
  XLSX files directly with `writexl` and removing the platform-dependent font
  unpacking and re-compression step.
- Retained `font_name` as a deprecated compatibility argument; supplied values
  now produce a warning and are not applied.
- Added strictly matched expression and metadata input files for a bundled
  10-sample TCGA-BRCA Normal workflow.
- Added a directly runnable bundled-example script and end-to-end tests for
  input alignment, default scoring, and workbook creation.
