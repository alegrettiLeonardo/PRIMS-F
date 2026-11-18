# Estado para continuação no Codex

## Base e objetivo

Repositório `alegrettiLeonardo/PRIMS-F`. O [PR #1](https://github.com/alegrettiLeonardo/PRIMS-F/pull/1)
foi incorporado ao `main` no commit `4572f95ee61795c852b8d87c0f142505dcd11378`.
A continuação com os diagnósticos e este handoff está na branch
`audit/chexal-continuation-20261118`, baseada nesse merge e entregue em novo PR draft.
O primeiro commit de implementação foi `c31521ed5bff1b4df90a184b44fb88765f3b8dd1`;
a base original era `16c9879cc8c134591bf0653e47d70e399377d8a9`.
Verifique a ponta das branches antes de trabalhar e preserve alterações locais.

Último pedido: avaliar Chexal para contrafluxo e continuar no WSL em
`/home/leonardoab/Documentos/PRIMS-F`. **Chexal ainda não está implementado.**

## O que já existe

- Fortran 2018; operadores em SI, derivadas, Fanning, gravidade e escala corrigidos;
  primitivas, Jacobianos, catenária, estacionário e contornos revistos.
- IDs de fechamento: 1 estratificado, 2 Bendiksen, 3 Bhagwat conforme transcrição
  Azevedo–Baliño–Burr 2015, 4 variante Montgomery. Bhagwat atual exige cocorrente.
- FV de primeira ordem, Euler implícito, Newton e LAPACK. Campanha usa `ROE`,
  `MATCHED`, `LAGGED`. MUSCL não está qualificado no novo caminho.
- **2185 verificações Fortran em 6 programas e 16 testes Python aprovados.**
  Isso não constitui validação física. Comandos e dependências: README/CI.
- Catenária da dissertação: Lp 57,4 m, declive 2°, X 6,435 m, Z 9,886 m,
  D 0,05254 m; Lr 12,546254 m, topo 76,473695° desde a horizontal. Perfil principal
  usa `rho_g=p/343^2`; `_ISO` é sensibilidade separada. Gás de entrada convertido
  nas condições padrão 101300 Pa / 293,15 K. Não misturar instalações/EOS.

## Achados que orientam a próxima alteração

1. Controle estável completa 2400 s. OSC barotrópico oscila com passos menores,
   mas amortece com passos maiores; não atribuir estabilidade apenas ao fechamento.
   Em 51 + 51 células, dt máximo0,125 s: 64 intervalos na janela 900–2400 s, período
   próximo 23 s. Dois refinamentos completos apresentam dez intervalos finais
   com amplitude quase constante: 51 + 51 / dt máximo 0,0625 s e 101 + 101 /
   dt máximo 0,125 s, ambos com período próximo de 23,27 s. As amplitudes ainda
   diferem entre discretizações. Leia `VALIDATION_RESULTS.md` e os envelopes;
   não reclassifique essas oscilações como validação física.
2. SS1/SS3 não completam os transientes. Há causas distintas: limite celular
   `ul=0`; traço gasoso `MATCHED` negativo com células admissíveis; caminho de Roe
   fora do domínio ao usar inclinação média. Veja `TRANSIENT_DOMAIN_LIMITS.md`
   e `validation/diagnostics/final_domain_summary.json`. Fluxo numérico negativo
   na junção não prova contrafluxo físico nas células.
3. A busca estacionária atual só varre `|logit(alpha)|<=28`: pode retornar falso
   `NO_ROOT` para alfa positivo abaixo de aproximadamente 6,9144e-13. O candidato
   em `validation/diagnostics/alpha_tail_candidate/` foi testado separadamente,
   **não aplicado**. Corrige apenas essa busca; não libera contracorrente.
4. Chexal 1992 p.262 / PDF 8 restringe contracorrente ao tubo vertical, gás subindo,
   líquido descendo; p.264 / PDF 10 só cocorrente horizontal. P.268 / PDF 14 exige
   seleção entre duas raízes abaixo da CCFL. Não extrapolar isso para a catenária.
   Referências e cuidados: `CHEXAL_1992_SCOPE.md`.
5. Baliño 2010 tabela 3 p.653 / PDF 11, caso 11: SS3 experimental **76 s**, confirmado
   visualmente no novo PDF; a dissertação traz 26 s. O erro percentual impresso
   na tabela também é inconsistente com 72,2/76. OSC caso 13: faixa 19,3–33,8 s;
   a dissertação informa 33,9 s. Não ajustar coeficientes para resolver divergências.

## Retomada das execuções interrompidas

A campanha entregue reúne **28 configurações: 13 completas, 13 falhas do solver
e duas interrupções da sessão de execução**, além de duas tentativas anteriores
arquivadas por timeout de 1800 s de relógio: 30 execuções no total. As interrupções
não são falhas físicas. Não há processo que esta entrega prometa continuar.

As duas repetições interrompidas pertencem a D_OSC / BAROTROPIC / ID3 / MATCHED /
LAGGED / ROE: 201 + 201 células, dt máximo 0,125 s, último registro em
1351,5568528 s; e 51 + 51 células, dt máximo 0,03125 s, último registro em
1329,8109900 s. Os dados das primeiras tentativas chegaram mais longe e permanecem
em `validation/runtime_limited`; não unir séries reiniciadas como se fossem uma.

Para fechar a campanha anterior, repita essas configurações com os fontes
congelados da entrega, limite de relógio suficiente e rótulo novo. O limite de
relógio não altera o tempo físico solicitado. O solver não oferece retomada por
checkpoint validada; iniciar novamente não deve ser chamado de continuação do
mesmo estado. Preserve a comparação por versão se aplicar o candidato de alfa.

## Mapa mínimo

| Trabalho | Arquivos |
|---|---|
| C0, Vgj, fonte, derivadas | `scr/module_closure_law.f90`, `scr/mod_void_helpers.f90` |
| Primitivas/Jacobiano/Roe | `scr/module_find.f90`, `scr/module_jacobian.f90`, `scr/module_matrix_roe.f90` |
| Ramo, estacionário, contornos | `scr/module_flow_regime.f90`, `scr/module_steady_state.f90`, `scr/module_boundary_conditions.f90` |
| Fluxos, junção, Newton | `scr/module_simulate_pipe_riser.f90`, `scr/module_roe_solver.f90` |
| Entradas/testes | `scr/module_config.f90`, `input/benchmarks/`, `tests/` |
| Reprodução | `scripts/run_validation.py`, `scripts/run_study.py`, `scripts/analyze_*.py`, `scripts/compare_validation.py` |
| Fontes e qualificação | `BENCHMARK_REFERENCE.md`, `benchmark_sources.json`, `IMPLEMENTATION_PHYSICS.md`, `VALIDATION_RESULTS.md` |

Os PDFs da conversa não são distribuídos pelo repositório. Localize-os no WSL
antes de transcrever fórmulas; não suponha que caminhos do ambiente anterior
existam no computador do usuário. Os manifests preservam caminhos históricos
como proveniência; scripts novos devem funcionar a partir da raiz do clone.

Use rótulos novos para novas execuções, preserve resultados anteriores e registre
hashes. Diferencie execução incompleta, timeout de infraestrutura, verificação
numérica, convergência e concordância experimental.

`make test` executa as suítes Fortran. Execute também
`python3 -m unittest discover -s tests -p 'test_validation_evidence.py' -v`
para os 16 testes Python. O pós-processamento requer NumPy, SciPy e Matplotlib;
mpmath é necessário apenas para regenerar as referências de alta precisão.

O runner `run_study --resume` não é um retry automático: pode reutilizar uma
tentativa encerrada com falha e recusa a pasta `INTERRUPTED`. O analisador
`analyze_osc_envelopes.py` ainda seleciona somente os padrões
`final_D_OSC_C3_BAROTROPIC_MATCHED_LAGGED_ROE_*` e
`refined_D_OSC_C3_BAROTROPIC_MATCHED_LAGGED_ROE_*`. Antes de comparar Chexal ou
rótulos novos, generalize essa seleção com opções explícitas e testes, sem
misturar fisicamente os conjuntos antigos. Para um smoke de um único caso,
use `run_validation.py`; reduzir o tempo em `run_study` ainda lança toda a matriz.
