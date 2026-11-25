# t1_seal_jac_20261122 — ferramentas e evidência da rodada P0/P1/P2

Relatório: [docs/PRIMSF_T1_SEAL_JACOBIAN_PERFORMANCE_QUALIFICATION.md](../../../docs/PRIMSF_T1_SEAL_JACOBIAN_PERFORMANCE_QUALIFICATION.md).

| Arquivo | Uso |
|---|---|
| `build.sh PROG.f90 OUT` | Liga um diagnóstico offline (`offline_api.f90`) aos objetos release do worktree |
| `seal_jac_check.f90` | Face 53 \| 54 em estados congelados. Modos: OPS (operadores A/B/C), JAC (implícito × Richardson da face completa × FD), SWITCH (troca H4/H7, J₊/J₋), NOISE (ruído do fluxo), RH (Rankine–Hugoniot da frente H7) |
| `newton_floor.f90` | Newton denso offline de um passo a partir de um estado aceito; localiza o piso do resíduo por célula e componente |
| `ckpt_fixture.py CKPT TEMPLATE OUT` | Fixture offline a partir de `checkpoint.dat` |
| `bench_summary.py RUNS OUT.json OUT.md [a:b …]` | Métricas P0.4 e diferenças de solução |
| `momentum_balance.py RUN [OUT]` | Balanço integral de momento por duto, riser e sistema (`momentum_ledger.csv`) |
| `cycle_metrics.py RUN [OUT]` | Eventos SEAL, PEN e BLOW; métricas por ciclo |
| `evidence/` | Logs e tabelas leves; a evidência bruta fica em `~/Documentos/PRIMS-F-t1seal-jac-results-20261122/` |

Ambiente das ferramentas offline:

```
PRIMS_PROBE_CLOSURE_RISER=5 PRIMS_JUNCTION_TRANSMISSION=T1 PRIMS_JUNCTION_SEAL=LIQUID_BRIDGE PRIMS_SEAL_STEPS=16
```
