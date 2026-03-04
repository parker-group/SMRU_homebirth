# Data Dictionary: Homebirth Analysis Dataset

This document describes the variables used in the formatted dataset file:

- [`data/homebirth_analysis.csv`](../data/homebirth_analysis.csv)

Each row represents **one participant / pregnancy observation**.

| Variable | Description | Units / Coding |
|----------|-------------|----------------|
| `ANC_code` | Pseudonymized participant identifier (unique; one row per participant) | ID code |
| `OR_site` | Pseudonymized outreach site / clinic identifier (used as the clustering unit in the mixed-effects model) | ID code |
| `home_birthYN` | Birth location indicator | `1` = home birth, `0` = facility birth |
| `hi` | Health insurance status code |  `1` = MFUND; `2` = Thai health insurance; `3` = no insurance  |
| `Age` | Maternal age | Years |
| `yrproject` | Study year / project phase identifier | Categorical |
| `MEthnicCAT` | Participant-reported ethnicity category | Categorical |
| `CanRead` | Literacy indicator | `1` = can read, `0` = cannot read |
| `Parity` | Number of prior births | Count |
| `Caesarean` | History of prior cesarean delivery | 'Yes' or 'No' |
| `hx_home_sum` | Number of prior home births | Count |
| `HowLongYear` | Number of years participant has lived in the study area | Years |
| `EnterDate` | Date of study entry | Date (MM/DD/YYYY) |
| `enrolEGA` | Gestational age at enrollment | Weeks |
| `EGA_birth` | Gestational age at birth | Weeks |
| `NOC_MX_max` | Total number of antenatal care (ANC) visits recorded for the pregnancy | Count |
| `prematureY_N` | Preterm birth indicator |  `1` = preterm (<37 weeks gestation); `2` = term birth  |

Notes  
- `ANC_code` and `OR_site` are pseudonymized identifiers created for the public dataset.  
- Additional variables used in the statistical models (e.g., `hi2`, `Age_bin25p`, parity-gated terms) are constructed within the analysis pipeline and are not included as columns in the public dataset.  
- Derived variables are created in `scripts/analysis_main.R`.
