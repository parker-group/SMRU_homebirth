# Data Dictionary: Homebirth Analysis Dataset

This document describes the variables used in the formatted dataset file:

- [`data/homebirth_analysis.csv`](../data/homebirth_analysis.csv)

Each row represents **one participant / pregnancy observation**.

---

## Identifiers

| Variable | Description |
|----------|-------------|
| `ANC_code` | Pseudonymized participant identifier (unique; one row per participant) |
| `OR_site` | Pseudonymized outreach site / clinic identifier (used as the clustering unit in the mixed-effects model) |

---

## Birth Outcome

| Variable | Description |
|----------|-------------|
| `home_birthYN` | Birth location indicator (`1` = home birth, `0` = facility birth) |
| `prematureY_N` | Preterm birth indicator (raw categorical coding as provided in the dataset) |

---

## Maternal Demographics

| Variable | Description |
|----------|-------------|
| `Age` | Maternal age in years |
| `MEthnicCAT` | Participant-reported ethnicity category |
| `CanRead` | Literacy indicator (`1` = can read, `0` = cannot read) |
| `hi` | Health insurance status code (raw categorical coding as provided in the dataset) |

---

## Reproductive History

| Variable | Description |
|----------|-------------|
| `Parity` | Number of prior births |
| `Caesarean` | History of prior cesarean delivery (raw categorical coding as provided in the dataset) |
| `hx_home_sum` | Number of prior home births |

---

## Pregnancy Care

| Variable | Description |
|----------|-------------|
| `NOC_MX_max` | Total number of antenatal care (ANC) visits recorded for the pregnancy |
| `enrolEGA` | Gestational age at enrollment (weeks) |
| `EGA_birth` | Gestational age at birth (weeks) |

---

## Residence and Time Variables

| Variable | Description |
|----------|-------------|
| `HowLongYear` | Number of years the participant has lived in the study area |
| `EnterDate` | Date of study entry |
| `yrproject` | Study year / project phase identifier |

---

## Notes

- `ANC_code` and `OR_site` are **pseudonymized identifiers** created for the public dataset.
- Additional variables used in the statistical models (e.g., `hi2`, `Age_bin25p`, parity-gated terms) are **constructed within the analysis pipeline** and are not included as columns in the public dataset.
- Derived variables are created in:

```
scripts/analysis_main.R
```
