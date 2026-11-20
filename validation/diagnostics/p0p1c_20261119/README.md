# Diagnósticos da continuação P0/P1 (19–20/11/2026)

Ferramentas isoladas; nenhuma altera o solver. Relatório:
[QUALIFICATION_P0P1_CONTINUATION_20261120](../../../docs/QUALIFICATION_P0P1_CONTINUATION_20261120.md).

| Arquivo | Uso |
|---|---|
| `SETTLING_CRITERIA_v2.json` | critérios pré-registrados (`813da75`) |
| `settling_criteria.py` | implementação pura dos critérios v2; usada por `../p0p1_20261119/analyze_refinement.py` |
| `ray_map.f90` | raio de Newton no estado terminal: λ* por bisseção, resíduo admissível, mapa de raízes cocorrentes Φ(u_g) ou traço de junção |
| `junction_bounds.f90` | fluxos de fase dos dois traços vs Roe e HLL; decomposição característica da dissipação gasosa; parcela do Δp centro‑a‑centro |
| `steady_face_check.f90` | face da junção com `junction_trace` `CENTER` e `STEADY` no mesmo estado |
| `collect_evidence.py` | índice de execuções: entradas, binários, status, passos, custo, hashes e seis balanços pelo auditor do ledger |
| `evidence/` | logs compactos dos replays, índice de execuções e resumos v2 |

Compilação dos drivers (objetos debug do mesmo checkout; `offline_api.f90` vem de
`../qualification2_20261119/globalization_flux/`):

```sh
objs=$(find tmp/check -name '*.o' ! -name main.o)
FL="-std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fopenmp -Itmp/check"
gfortran $FL -c validation/diagnostics/qualification2_20261119/globalization_flux/offline_api.f90 -o offline_api.o
gfortran $FL -I. validation/diagnostics/p0p1c_20261119/ray_map.f90 offline_api.o $objs -llapack -lblas -o ray_map
./ray_map PASTA_COM_last_rejected_state.dat_E_probe_parameters.dat CG_SI
```

`CG_SI` é a velocidade do som da EOS **da execução** (343 m/s BAROTROPIC;
√(287·293,15) m/s ISOTHERMAL). `offline_api` fixa malha 51+51 e Bhagwat no riser.

Análise de acomodação/convergência:

```sh
python3 validation/diagnostics/p0p1_20261119/analyze_refinement.py --runs PASTA_COM_N*_dt* --output NOVA
python3 validation/diagnostics/p0p1c_20261119/collect_evidence.py RESULTADOS SAIDA.json
```
