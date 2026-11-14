| run | tasks x runs | reward (95 % CI) | tumour cleared | gap to reference (95 % CI) | tags (a run can carry several) |
|---|---|---|---|---|---|
| Reference scientist | 10 x 3 | 0.64 (0.59 to 0.69) | 30/30 | - | - |
| Claude Fable 5.2 | 10 x 3 | 0.59 (0.52 to 0.66) | 30/30 | 0.06 (0.01 to 0.11) | avoidable_harm 9 |
| Claude Opus 5.5 | 10 x 3 | 0.56 (0.48 to 0.63) | 30/30 | 0.08 (0.03 to 0.16) | avoidable_harm 4, too_toxic 2 |
| GPT-6.1 Sol | 10 x 3 | 0.15 (0.05 to 0.27) | 24/30 | 0.49 (0.39 to 0.59) | too_toxic 18, avoidable_harm 5, clone_not_covered 3, wrong_drugs 2, missed_clone 1 |
| Claude Haiku 5.5 | 10 x 2 | 0.28 (0.16 to 0.42) | 18/20 | 0.36 (0.24 to 0.47) | too_toxic 11, avoidable_harm 5, wrong_drugs 1, missed_clone 1 |
| Claude Haiku 5.5 with its notes | 10 x 2 | 0.30 (0.17 to 0.43) | 16/20 | 0.35 (0.23 to 0.46) | too_toxic 7, avoidable_harm 4, missed_clone 2, clone_not_covered 2, suboptimal 1 |
| Claude Haiku 5.5 with its second notes | 10 x 2 | 0.21 (0.09 to 0.33) | 14/20 | 0.44 (0.31 to 0.56) | too_toxic 6, avoidable_harm 6, clone_not_covered 5, missed_clone 1 |
| Hedge | 10 x 1 | 0.05 (0.00 to 0.12) | 4/10 | 0.59 (0.52 to 0.66) | missed_clone 5, too_toxic 5, avoidable_harm 1, wrong_drugs 1 |
| Blind guess | 10 x 1 | 0.04 (0.00 to 0.12) | 4/10 | 0.61 (0.53 to 0.68) | too_toxic 8, missed_clone 6, avoidable_harm 1 |
| Guideline regimen | 10 x 1 | 0.00 (0.00 to 0.01) | 0/10 | 0.64 (0.59 to 0.69) | missed_clone 10, too_toxic 3 |
| Three-day screen | 10 x 1 | 0.00 (0.00 to 0.01) | 1/10 | 0.64 (0.59 to 0.69) | missed_clone 8, too_toxic 2, wrong_drugs 1 |
| Stop early | 10 x 1 | 0.00 (0.00 to 0.00) | 0/10 | 0.64 (0.59 to 0.69) | no_regimen 10 |
| Do nothing | 10 x 1 | 0.00 (0.00 to 0.00) | 0/10 | 0.64 (0.59 to 0.69) | no_regimen 10 |
