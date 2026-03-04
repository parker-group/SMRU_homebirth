# Data Dictionary: Homebirth Analysis Dataset

This document describes the variables used in the formatted dataset file:

  * `data/homebirth_analysis.csv`

Variable  Description  
`ANC_code` Pseudonymized participant ID (unique; one row per participant)  
`OR_site` Pseudonymized outreach site / clinic ID (used for clustering / random intercept)  
`home_birthYN` Home birth indicator (1 = home birth, 0 = not home birth)  
`hi` Health insurance status code (raw categorical coding as provided in the dataset)  
`Age` Maternal age in years  
`yrproject` Project year / study period identifier  
`MEthnicCAT` Ethnicity category (categorical)  
`CanRead` Literacy indicator (1 = can read, 0 = cannot read)  
`Parity` Number of prior births (parity)  
`Caesarean` Prior cesarean history (raw categorical coding as provided in the dataset)  
`hx_home_sum` Number of prior home births  
`HowLongYear` Years lived in current area  
`EnterDate` Study entry date (MM/DD/YYYY as stored in the dataset)  
`enrolEGA` Gestational age at enrollment (weeks)  
`EGA_birth` Gestational age at birth (weeks)  
`NOC_MX_max` Number of antenatal care (ANC) visits recorded for the pregnancy  
`prematureY_N` Preterm birth indicator (raw categorical coding as provided in the dataset)  

Notes  
`ANC_code` and `OR_site` are pseudonymized identifiers created for public release. Derived variables used in modeling (e.g., `hi2`, `Age_bin25p`, parity-gated terms) are created within `scripts/analysis_main.R` and are not included as columns in the public dataset.
