# pipe_riser (Fortran Moderno)

Este projeto e um **esqueleto compilavel** em Fortran (F2008) que segue a modularizacao pedida:

- `module_simulate_pipe_riser`  (simulatePipeRiser.py)
- `module_jacobian`            (jacobianMatrix.py)
- `module_boundary_conditions` (boundaryConditions.py)
- `module_catenary`            (catenary.py)

## Estado da traducao

- `module_jacobian`: portado com fidelidade para o Jacobiano do fluxo.
- `module_boundary_conditions`: portado com NSCBC e **ghost cells por extrapolacao linear**.
- `module_simulate_pipe_riser`: portado no nivel de estrutura/loops; inclui OpenMP no loop de interfaces.

### Pontos que ainda precisam do port 100% fiel

1. `module_find`: hoje a recuperacao de `ul, ug` e um **placeholder** (o Python resolve o fechamento drift-flux). Para fidelidade total, porte:
   - `closureLaw.py` (Swananda + derivadas, padroes)
   - `find_ulug_from_uj_genflw_simp` e rotinas auxiliares (brentq etc.)
   - `criteriaStabFlow.py` e `voidFraction.py` conforme necessario.

2. `module_roe_solver`: MUSCL + |A| (Roe) esta simplificado. Para fidelidade total, porte `roeSolver.py`.

3. IO: o Python grava HDF5 e tem plots; aqui deixei somente `print`.

## Compilacao

No diretorio raiz:

```bash
make
make run
```

Dependencias:
- gfortran
- BLAS/LAPACK (`-llapack -lblas`)
- OpenMP (opcional, ja habilitado no Makefile)

## Proximo passo sugerido

1) Substituir o placeholder de `compute_primitive_variables` em `module_find` pelo port fiel do Python.
2) Portar `roeSolver.py` para `module_roe_solver`.
3) Validar uma simulacao simples (caso benchmark) comparando com o Python.
