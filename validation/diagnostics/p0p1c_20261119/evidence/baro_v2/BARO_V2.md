# BAROTROPIC OSC refinement, criteria v2

| run | settled v2 | period s | semi-amplitude Pa | mean Pa | window s |
|---|---|---:|---:|---:|---|
| N51_dt0.03125 | True | 23.344137 | 15354.38 | 266775.6 | 638.4-871.8 |
| N51_dt0.0625 | True | 23.265321 | 13764.28 | 266573.6 | 636.5-869.1 |
| N51_dt0.125 | True | 23.108397 | 9905.70 | 266174.6 | 355.6-586.7 |
| N101_dt0.0625 | True | 23.422593 | 16002.43 | 266957.2 | 640.5-874.7 |
| N101_dt0.125 | True | 23.266063 | 12762.72 | 266541.4 | 636.7-869.3 |

| axis | coarse | fine | period | amplitude | mean | v2 converged | status |
|---|---|---|---:|---:|---:|---|---|
| space | N51_dt0.125 | N101_dt0.125 | 0.0068 | 0.2239 | 0.0014 | False | FAIL |
| time | N51_dt0.125 | N51_dt0.0625 | 0.0067 | 0.2803 | 0.0015 | False | FAIL |
| time | N51_dt0.0625 | N51_dt0.03125 | 0.0034 | 0.1036 | 0.0008 | False | FAIL |
| combined | N51_dt0.125 | N101_dt0.0625 | 0.0134 | 0.3810 | 0.0029 | False | FAIL |
