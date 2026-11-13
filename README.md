# Seeing Behind the Meter

Team 01 - Binghamton University Power Electronics Research Group
2026 Avangrid Energy Innovation Hackathon

MATLAB tool that reads hourly AMI loads and weather, estimates which homes have air conditioning, electric heating, heat pumps, EVs, water heaters, pools, hot tubs and solar, segments the customers, and ranks homes for demand response and a behind-the-meter battery.

## Methodology

1. Features (`step2`): 82 numbers per home: size, weather response, charging blocks and spikes, zero readings, daily and monthly load shape.
2. Weather model: daily kWh = b0 + bH*HDD60 + bC*CDD65, least squares, one fit per home. Air conditioning and electric heating are thresholds on bC and bH.
3. Threshold rules (`step3`): water heater, pool, hot tub, solar, solar + battery and other signatures. Each rule gives a 0-1 confidence from a sigmoid around its threshold.
4. Trained models (`step3`): L1-regularised logistic regression for EV and heat pump, trained on Avangrid's customer lists with 5-fold cross-validation. The lists hold about 25% of the real customers, so the top N = labels/0.25 homes are flagged and precision is estimated as (labelled homes in top k) / (0.25 k).
5. Segmentation (`step4`): k-means (k = 6) on the normalised daily and monthly load shape, plus rule-based business segments.
6. Demand flexibility and battery (`step5`): flexible kW per home in the 100 system peak hours; one-year hourly simulation of the 20.16 kWh / 7.68 kW battery with a controller that uses past data only; score, tiers, sizing, benefit and cost.

## Setup

- MATLAB R2024a with the Statistics and Machine Learning Toolbox (`requirements.txt`).
- Put Avangrid's data folders, the weather file and the customer list file in one folder and set `root` and `src` at the top of `step1_load_data.m` and `step1b_add_folder.m`.
- Set `wd` (and `nd` in `matlab_v2`) to the folder where results are saved.

## Execution

Run the scripts in `matlab` in this order.

| Script | What it does | Output | Time |
| --- | --- | --- | --- |
| `step1_load_data.m` | Reads csv files, weather and customer lists | `ami.mat` | 20-30 min |
| `step1b_add_folder.m` | Adds the corrected csv files | `ami.mat` | 3 min |
| `step2_features.m` | Features per home | `features.mat` | 30 s |
| `step3_classify.m` | Technology flags and confidence scores | `class.mat` | 15 s |
| `step4_segments.m` | Clusters and business segments | `seg.mat` | 10 s |
| `step5_bess.m` | Demand flexibility, battery simulation, tiers, economics | `bess.mat`, `Avangrid_results.xlsx` | 2 min |
| `step6_figures.m`, `step7_evidence.m` | Figures | `figs\` | 30 s |
| `step8_model_summary.m` | Summary tables of the weather model | printed | 5 s |

`Avangrid_results.xlsx` has one row per customer: flags, confidence scores, cluster, segments, flexible kW, battery tier, score, units, benefit and cost.

The scripts in `matlab_v2` are tests that do not change the results: `tod1_fingerprint.m` (weather model per 3-hour block), `ev_hp_confidence.m` (charging sessions, heat-pump checks), `model_compare.m` (other model types, feature families, weather).

## Assumptions

Battery $1,000/kWh installed, avoided service upgrade $8,000 (new load only), peak relief $100/kW-yr, energy 28 c/kWh on-peak and 10 c/kWh off-peak, 15 years, 7% discount rate.
