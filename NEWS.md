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
