# NP-3989 Table A-1: 1986 and 1992 transcriptions (steam-water L, declared properties)

| case | P psia | D ft | jf, jg | table C0 / Vgj / alpha | identity residual (bound) | V1986 C0 / Vgj / alpha | V1992 C0 / Vgj / alpha | version | cause | status |
|---:|---:|---:|---|---|---|---|---|---|---|---|
| 1 | 14.7 | 0.05 | 5, 10 | 1.2037 / 0.5979 / 0.5361 | +8.8e-06 (1.4e-04) | 1.2038 / 0.5984 / 0.5361 | 1.2038 / 0.5984 / 0.5361 | both (upflow form identical) | reproduced | PASS |
| 2 | 1000 | 0.05 | 5, 10 | 1.1116 / 0.1410 / 0.5947 | -1.2e-05 (1.3e-04) | 1.1107 / 0.1373 / 0.5953 | 1.1107 / 0.1373 / 0.5953 | both | unexplained: not a version difference (upflow form identical); see diagnostic fit | FAIL |
| 3 | 1000 | 1 | 5, 10 | 1.1119 / 0.2234 / 0.5914 | -4.2e-04 (1.3e-04) | 1.1109 / 0.2263 / 0.5921 | 1.1109 / 0.2263 / 0.5921 | both | printed row violates alpha(C0 j+Vgj)=jg beyond its rounding bound (defective reference row) | FAIL (reference inconsistent) |
| 4 | 14.7 | 1 | 5, 10 | 1.1922 / 0.9054 / 0.5323 | +1.1e-04 (1.4e-04) | 1.1922 / 0.9063 / 0.5322 | 1.1922 / 0.9063 / 0.5322 | both (upflow form identical) | reproduced | PASS |
| 5 | 14.7 | 1 | -5, -10 | 1.0400 / 5.1172 / 0.9538 | +1.5e-04 (1.3e-04) | 1.0399 / 5.1162 / 0.9539 | 1.5869 / 11.3909 / 0.8056 | 1986 | downflow revised in 1992 (C3=2(C10/2)^B2, C9, Eq.3-4b); 1986 equations reproduce the row | PASS (version difference demonstrated) |
| 6 | 1000 | 1 | -5, -10 | 1.0498 / 4.8444 / 0.9171 | +1.2e-04 (1.3e-04) | 1.0514 / 4.8318 / 0.9142 | 1.0921 / 4.6374 / 0.8515 | both | unexplained: not a version difference (upflow form identical); see diagnostic fit | FAIL |
| 7 | 1000 | 0.05 | -5, -10 | 1.2948 / 1.8775 / 0.5700 | -3.6e-05 (1.3e-04) | 1.2927 / 1.8363 / 0.5697 | 1.2934 / 1.8193 / 0.5688 | both | unexplained: not a version difference (upflow form identical); see diagnostic fit | FAIL |
| 8 | 14.7 | 0.05 | -5, -10 | 1.3036 / 4.6590 / 0.6714 | -5.0e-05 (1.3e-04) | 1.3026 / 4.6665 / 0.6724 | 1.3037 / 4.6584 / 0.6713 | both | reproduced by both versions | PASS |

Case 2 diagnostic (not a correction): an effective Re_f = 81956 in B1/C3 (printed 184219) gives max rel. error 5.7e-03 (C0 1.1107, Vgj 0.1402, alpha 0.5952).

Case 3 diagnostic (not a correction): an effective Re_f = 83603 in B1/C3 (printed 3684867) gives max rel. error 1.3e-02 (C0 1.1109, Vgj 0.2263, alpha 0.5921).
