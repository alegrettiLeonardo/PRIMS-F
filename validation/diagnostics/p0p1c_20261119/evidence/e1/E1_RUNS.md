# Roe vs E1: dynamic replays (frozen operator, same deck)

| reference | E1 run | t ref (s) | t E1 (s) | reason E1 | acc/rej ref | acc/rej E1 | dt_min E1 | cpu ref | cpu E1 | max/RMS dp_in | max dMl, dMg (kg) | A*(0.25) ref | A*(0.25) E1 |
|---|---|---:|---:|---|---|---|---:|---:|---:|---|---|---:|---:|
| SS1_CENTER | SS1_E1_CENTER | 78.577 | 80.861 | NEWTON_GLOBALIZATION_FAILURE | 706/29 | 678/23 | 2.46e-06 | 212.7 | 176.7 | 0.0200/0.0063 | 1.525, 0.0016 | - | - |
| SS1_CENTER | SS1_E1_STEADY | 78.577 | 79.886 | NEWTON_GLOBALIZATION_FAILURE | 706/29 | 669/23 | 1.03e-06 | 212.7 | 754.2 | 0.0124/0.0039 | 0.945, 0.0010 | - | - |
| SS1_STEADY | SS1_E1_STEADY | 78.838 | 79.886 | NEWTON_GLOBALIZATION_FAILURE | 661/23 | 669/23 | 1.03e-06 | 628.7 | 754.2 | 0.0096/0.0030 | 0.731, 0.0008 | - | - |
| B_ISO_MASS | ISO_E1_CENTER | 223.645 | 223.883 | NEWTON_GLOBALIZATION_FAILURE | 1829/25 | 1825/24 | 2.21e-07 | 311.4 | 329.1 | 0.0102/0.0021 | 0.731, 0.0015 | 1.159e+04 | 1.159e+04 |
| B_ISO_MASS | ISO_E1_STEADY | 223.645 | 223.842 | NEWTON_GLOBALIZATION_FAILURE | 1829/25 | 1823/24 | 3.69e-07 | 311.4 | 1581 | 0.0054/0.0011 | 0.386, 0.0008 | 1.159e+04 | 1.158e+04 |
| B_ISO_MASS_N101 | ISO_E1_CENTER_N101 | 249.400 | 226.307 | OPERATOR_TERMINATED_COST_STAGNATION | 2037/24 | 2423/-1 | 4.56e-05 | 811.3 | - | 0.0037/0.0008 | 0.258, 0.0005 | 1.133e+04 | 1.13e+04 |
| B_ISO_MASS_N101 | ISO_E1_STEADY_N101 | 249.400 | 226.884 | OPERATOR_TERMINATED_COST_STAGNATION | 2037/24 | 2160/-1 | 3.16e-05 | 811.3 | - | 0.0019/0.0004 | 0.137, 0.0003 | 1.133e+04 | 1.132e+04 |

Last accepted junction state:

| run | p_in Pa | mdot_g junction kg/s | alpha pipe end | alpha riser base | holdup base |
|---|---:|---:|---:|---:|---:|
| SS1_CENTER | 275139 | 5.447e-11 | 0.5600 | 0.0188 | 0.9812 |
| SS1_E1_CENTER | 278302 | 1.887e-11 | 0.5586 | 0.0152 | 0.9848 |
| SS1_E1_STEADY | 277261 | 3.409e-12 | 0.5593 | 0.0164 | 0.9836 |
| SS1_STEADY | 275678 | 4.664e-11 | 0.5600 | 0.0181 | 0.9819 |
| B_ISO_MASS | 265513 | 1.727e-11 | 0.5505 | 0.0057 | 0.9943 |
| ISO_E1_CENTER | 264853 | 9.376e-12 | 0.5489 | 0.0060 | 0.9940 |
| ISO_E1_STEADY | 265429 | 1.686e-12 | 0.5496 | 0.0057 | 0.9943 |
| B_ISO_MASS_N101 | 256165 | 2.451e-12 | 0.5510 | 0.0049 | 0.9951 |
| ISO_E1_CENTER_N101 | 273433 | 8.043e-06 | 0.5494 | 0.0024 | 0.9976 |
| ISO_E1_STEADY_N101 | 275839 | 4.636e-06 | 0.5498 | 0.0012 | 0.9988 |
