# Chexal–Lellouche em todo o riser: ferramentas e evidências (21/11/2026)

Relatório: [CHEXAL_LELLOUCHE_RISER_QUALIFICATION](../../../docs/CHEXAL_LELLOUCHE_RISER_QUALIFICATION.md).
Nenhuma ferramenta entra no solver.

| Arquivo | Uso |
|---|---|
| `chexal92_ref.py` | Transcrição Python independente das equações de 1992 (ar–água; cocorrente em qualquer ângulo e contracorrente vertical HIGH). Contracorrente inclinada, horizontal negativa e gás ↓ com líquido ↑ levantam `DomainError` |
| `np3989_sample_check.py` | EPRI NP-3989 Rev. 1, Tabela A-1 (vapor–água vertical) contra uma transcrição independente de 1992 com L de vapor–água; banda de propriedades declarada |
| `core_crosscheck.f90/.py` | Núcleo Fortran (`chexal_coefficients_si`) contra `chexal92_ref.py` em 945 estados |
| `global_jacobian_check.f90`, `prepare_state.py`, `run_global_jacobian.sh` | Jacobiano global montado (51+51, MATCHED, Roe, LAGGED) contra diferenças centradas do resíduo completo em estados aceitos; modo `JUNCTION` para o traço da junção (dissipação no estado aceito, como o solver) |
| `fields_closure_check.py` | Resíduo Chexal recalculado em todos os estados aceitos de riser gravados (`fields.csv`) |
| `domain_boundary_probe.py` | Decodifica o status terminal e recupera j_l no estado aceito, localizando a fronteira de líquido estagnado |
| `run_campaign.py` | Campanha P1 (decks primários com override de `closure_riser`, `tol_res=1e-10`, opção `--grid N:dt`) sob `junction_constitutive/run_guarded.py` |
| `summarize_campaign.py` | Classificação, seis balanços (`audit_ledger`), eventos, inventários e ciclos |

Compilação das ferramentas Fortran (objetos de `make check_all`):

```sh
objs=$(find tmp/check -name '*.o' ! -name main.o); mkdir -p tmp/cx
for t in global_jacobian_check core_crosscheck; do
  gfortran -std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fopenmp \
    -Itmp/check -Jtmp/cx validation/diagnostics/chexal_riser_20261121/$t.f90 $objs -llapack -lblas -o tmp/cx/$t
done
```

`global_jacobian_check` fixa a geometria D (Lr = 12,546253993197080 m), ρ0 = 997,05, P0 = 101325
e as propriedades dos decks D; ω_u = max(jl, jg) e c_g vêm do deck via `prepare_state.py`.

Evidências em `evidence/`: `NP3989_TABLE_A1.*`, `CORE_CROSSCHECK.json`, `global_jacobian/`,
`FIELDS_CLOSURE_CHECK.json`, `P1_CAMPAIGN.*`, `P1_REFINE.*`, `DOMAIN_PROBE.json`, `REGRESS_BHAGWAT.json`.
