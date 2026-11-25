# Junção: limites constitutivos, referência independente e E2 (20/11/2026)

Ferramentas isoladas; nenhuma entra no solver. Relatório:
[PRIMSF_JUNCTION_CONSTITUTIVE_QUALIFICATION](../../../docs/PRIMSF_JUNCTION_CONSTITUTIVE_QUALIFICATION.md).

| Arquivo | Uso |
|---|---|
| `junction_limits.f90` | P0.1: células e face da junção (α, ρ_g, u_g, j_g, U_g, p, j_l, C0, Vgj, autovalores; F_conv, F_diss, F_total por operador) e limites LIM1 (u_g=0, α>0), LIM2 (α→0) e LIM3 (j_g→0±) com as rotinas de fechamento do solver; linhas `FACE_E2` com `junction_riser_limit='ALPHA0'` |
| `summarize_limits.py` | tabelas P0.1 (`evidence/P01_LIMITS.*`) |
| `junction_riemann_ref.f90` | P0.2: referência independente (problema de Riemann com camada estacionária; T1 conservativa, T2 pressão contínua; curvas integrais; Newton). Modo `TESTS`: uniforme, partição, ondas pequenas |
| `summarize_ref.py` | tabela P0.2 (`evidence/P02_REFERENCE.*`) |
| `e2_newton_probe.f90` | Gate P0.3: reproduz a direção de Newton do passo recusado de um replay E2 e o mapeamento do traço do riser |
| `run_guarded.py` | execução com limites de custo pré-declarados (COMPLETED, FAIL, TIMEOUT, COST_STAGNATION) |
| `p11_time_convergence.py` | P1.1: janela pré-registrada (ciclos com início em [15, 100] s), ajuste modal, ordens observadas |
| `ss3_countercurrent_local.py` | P1.2: relação de deriva assinada na célula crítica de SS3 |

Compilação (objetos debug do mesmo checkout, `make check_all`):

```sh
objs=$(find tmp/check -name '*.o' ! -name main.o)
FL="-std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fopenmp -Itmp/check"
mkdir -p tmp/jc
gfortran $FL -Jtmp/jc -c validation/diagnostics/qualification2_20261119/globalization_flux/offline_api.f90 -o tmp/jc/offline_api.o
for t in junction_limits junction_riemann_ref e2_newton_probe; do
  gfortran $FL -Itmp/jc -Jtmp/jc validation/diagnostics/junction_constitutive/$t.f90 tmp/jc/offline_api.o $objs -llapack -lblas -o tmp/jc/$t
done
./tmp/jc/junction_limits FIXTURE C_GAS JUNCTION_FLUX JUNCTION_TRACE [RISER_LIMIT]
./tmp/jc/junction_riemann_ref FIXTURE C_GAS JUNCTION|TESTS
./tmp/jc/e2_newton_probe FIXTURE C_GAS JUNCTION_FLUX JUNCTION_TRACE RISER_LIMIT
```

`C_GAS` é a velocidade do som da EOS da execução (343 m/s BAROTROPIC;
290,0587009555135 m/s ISOTHERMAL). `offline_api` fixa a malha 51+51 e o Bhagwat no riser.
