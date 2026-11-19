# Verificação da seleção no oleoduto

A [implementação](PIPE_REGIME_BURR15.md) usa a Eq. (15), p. 3, conferida no
PDF de Burr. A transição reversa e a seleção automática de produção
ficaram **BLOCKED**. O detector calcula perda de estabilidade da solução
estratificada com a escrita de diagnóstico fora do operador. As
[fontes](BURR15_SOURCES.md) registram a diferença do fator de altura
em relação à expressão inicial, a inconsistência entre as densidades
nas Eqs. (15) e (53)–(54) e a ausência de retorno aplicável.

## Versão e regressões

- Operador `8a85e37319995b37795d2c798d6f865b863b1834`, descendente direto de
  `ee5810e`. Trabalho local em `audit/burr15-regime-20261119`.
- Build release com GNU Fortran 13.3, Fortran 2018 e `-O2`, LAPACK/BLAS.
  Binário da campanha
  `6b3164d2a18e08e173c1e65e5058bfc024d648a94b64fa10c7bc8ab8a24d37e5`.
  Build debug com `make check_all`, checks, backtrace e traps.
- `make test FC=gfortran LIBS='-llapack -lblas'` **PASS**, 3.573 verificações
  em 13 programas, incluindo 126 para o detector e os modos fixos.
  A suíte da CI Python passou 19 casos. Os verificadores independentes e a
  análise de campos passaram 20 casos adicionais. A
  [regressão de CLI](../validation/diagnostics/burr15_config_cli_20261119/)
  passou 18 verificações, incluindo deriva fixa horizontal, admissibilidade
  Chexal e consistência da saída final.
- [Regressão real do legado](../validation/diagnostics/burr15_legacy_regression_20261119/)
  com D_Stable e D_OSC, 51+51 células, 1 s, deu `monitor.csv` e `fields.csv`
  byte a byte iguais a `ee5810e`. Em `results_*.dat`, a recuperação das
  velocidades do riser passou a usar o ID simulado. Em D_Stable/D_OSC, os
  maiores deltas de `ug` na saída final foram 0,140665/0,290124 m/s.

## Série inicial

Oito execuções BAROTROPIC/Bhagwat no riser, `MATCHED/LAGGED/ROE`, com
geometria D e vazões dos decks. Controles idênticos em cada par,
51+51 células, 2400 s solicitados e teto de passo 0,125 s.

| Caso | Tempo atingido nos dois modos (s) | Passos aceitos/rejeitados por modo | Resultado | Identidade FIXED × diagnóstico |
|---|---:|---:|---|---|
| D_Stable | 2400 | 19222/0 | **PASS** em término/estacionariedade | CSVs físicos idênticos |
| D_OSC | 2400 | 19222/0 | **PASS** em término, **FAIL** em acomodação dessa janela | CSVs físicos idênticos |
| D_SS3 | 97,131852 | 809/24 | **FAIL** numérico 6000 | CSVs e snapshot final de rejeição idênticos |
| D_SS1 | 78,571788 | 662/24 | **FAIL** numérico 6000 | CSVs e snapshot final de rejeição idênticos |

O controle estável teve pressão média de base 501752,340 Pa na janela
900–2400 s e excursão de aproximadamente 6e-11 Pa. Os máximos erros
integrais líquido/gás foram 5,0783e-4/6,5445e-7 kg, equivalentes a
6,1911e-6/2,2332e-6 dos inventários iniciais. Pela massa processada
líquida de cada fase, as razões são 1,6479e-7/1,5139e-7. O critério
de 1e-6 por inventário inicial falha em ambos os modos. A identidade
entre as trajetórias torna esse desvio comum e o teste de regressão passa.

Em OSC, os onze picos cercados por mínimos delimitaram dez ciclos
de 2154,820 a 2385,850 s, com período de 23,1030 s e semi-amplitude
média de 9705,75 Pa. A semi-amplitude mudou 2,398% entre o primeiro
e o último pico, acima do limite de 1%. A leitura preliminar por
cruzamentos 20/80 estimou 23,0939 s e variação de 2,149%. Ambos os
métodos indicam uma janela ainda em acomodação.

Os máximos erros integrais OSC foram 9,5815e-10 kg para líquido e
1,8010e-10 kg para gás. As razões por inventário inicial são
1,3331e-11/1,0177e-9 e pela massa processada líquida são
3,3756e-13/2,4772e-11. O maior erro por ciclo é
1,7138e-11/6,6322e-13 kg (líquido/gás).
A diferença entre os dois modos em qualquer métrica dos CSVs é zero
em toda a cobertura comum. A comparação entre regimes comutados permanece
**BLOCKED** pelo modelo de transição.

| Caso diagnóstico | Estados do oleoduto avaliados | Perdas apontadas pela Eq. (15) | Margem mínima e máxima (m²/s²) |
|---|---:|---:|---:|
| D_Stable | 481644 | 481644 | 0,10395 a 1,40013 |
| D_OSC | 416058 | 416058 | 0,43828 a 5,10369 |
| D_SS3 | 19686 | 17307 | −0,07426 a 3,53644 |
| D_SS1 | 15249 | 2202 | −0,04302 a 1,74327 |

A reavaliação Python independente encontrou zero divergências de
status ou classificação. A maior discrepância de limiar/margem
quadrática foi 2,887e-15 m²/s². Em D_Stable e D_OSC, o indicador já
estava positivo em `t=0` e continuou positivo em todos os estados
amostrados. Isso limita o uso da classificação local como evidência
de uma transição de regime ou de uma classe global.

Nos replays dos dois snapshots terminais, o resíduo salvo foi reproduzido
exatamente. Em SS1, a junção `MATCHED` recusou as 18 propostas de Newton
entre `lambda=1` e `2**(-17)` com `jg*` de −6,055e-9 a −3,240e-14 m/s
e código 3012. Em SS3, a primeira recusa foi `NO_ROOT` na célula 69
com código 1014 e gás positivo na junção. O código 6000 registra o
esgotamento da busca linear. O gás negativo de SS1 pertence às propostas
recusadas. A duração atingida em ambos os casos é insuficiente para
qualificar período, o mecanismo físico ou a divergência experimental
de 26 s (dissertação, Tabela 6.7, p. 136) e 76 s (Baliño et al. 2010,
Tabela 3, caso 11, p. 653).

Nos estados amostrados de OSC, nenhuma velocidade de fase foi negativa
e o menor fluxo gasoso bruto de base amostrado foi
1,7701e-3 kg/s. As frações de vazio variaram de 0,2665 a 0,5512.
Os cruzamentos espaciais dos níveis 0,4/0,5/0,6 foram
7584/871/0, excluído o salto da junção. São observações de campo
que acompanham a oscilação numérica com gás atravessando a base.
A janela de 900 a 2400 s contém 66 episódios compatíveis com
acumulação/compressão, 49 excursões de produção e 65 episódios
compatíveis com blowdown. O fluxo gasoso amostrado permaneceu acima do
limiar inferior do detector de interrupção/retomada. O protocolo
identifica esses episódios pelas séries de pressão, inventários e
fluxos, com os limiares da análise registrados.
A identificação experimental do mecanismo e de frentes materiais
permanece **NOT_VERIFIED**.

## Refinamentos do controle estável

Dois eixos, dois níveis em cada. Foram mantidos o binário e as condições
físicas. Nenhum critério de ordem assintótica foi aplicado.

| Malha por trecho | Teto / passo mediano amostrado (s) | Aceitos/rejeitados | Pressão média na base (Pa), 900–2400 s | Máx. erro de massa L/G (kg) | Máx. erro / inventário inicial L/G |
|---|---:|---:|---:|---:|---:|
| 51+51 | 0,125 / 0,125 | 19222/0 | 501752,340 | 5,0783e-4 / 6,5445e-7 | 6,1911e-6 / 2,2332e-6 |
| 101+101 | 0,125 / 0,125 | 19222/0 | 501786,781 | 1,2509e-4 / 7,2690e-7 | 1,5250e-6 / 2,4804e-6 |
| 51+51 | 0,0625 / 0,0625 | 38418/0 | 501752,350 | 1,5479e-4 / 1,2851e-6 | 1,8871e-6 / 4,3853e-6 |

Os dois refinamentos completaram 2400 s e preservaram a estacionariedade.
O maior erro por massa processada líquida, entre os três níveis, foi
1,6479e-7 em líquido e 2,9729e-7 em gás. O teste de erro normalizado
pelo inventário inicial excede 1e-6 em todas as configurações.

## OSC após refinamento e extensão

As rodadas adicionais usam o modo fixo no mesmo binário. A série inicial
estabeleceu a identidade com o detector em toda a trajetória observada
para cada caso. Os onze picos cercados por mínimos são a referência
nas três janelas abaixo.

| Malha por trecho | Passo mediano amostrado (s) | Tempo atingido (s) | Passos aceitos/rejeitados | Período médio (s) | Semi-amplitude média (Pa) | Variação entre picos extremos |
|---|---:|---:|---:|---:|---:|---:|
| 101+101 | 0,125 | 2400 | 19222/0 | 23,266397 | 12773,70 | +0,02128% |
| 51+51 | 0,0625 | 2400 | 38418/0 | 23,265639 | 13772,54 | +0,00052% |
| 51+51 | 0,125 | 3600 | 28822/0 | 23,108640 | 9915,43 | +0,00641% |

Para os dez ciclos analisados, as médias e excursões não nulas dos 13
canais também ficaram dentro de 1%. A maior mudança relativa de excursão
entre as primeiras e últimas cinco janelas foi de 0,0183% no refino
espacial, 0,00119% no temporal e 0,0118% na extensão. Nesta última, os
limites dos dez intervalos são 3356,386–3587,473 s.

Na janela comum 900–2400 s, a excursão pico a pico da pressão de base
aumentou 30,30% no refinamento espacial e 40,47% no temporal comparados
ao caso grosso. Entre as janelas acomodadas, com a extensão como
referência, a semi-amplitude do caso refinado em malha é 28,83% maior,
e a do refinado em passo é 38,90% maior. As médias de base na janela
comum foram 265946,1 Pa no caso grosso, 266388,6 Pa no refino espacial
e 266510,8 Pa no temporal. Os períodos ficaram próximos de 23 s nas
três configurações. Há sensibilidade à discretização. Com dois níveis
por eixo, convergência assintótica de amplitude, forma de onda e
defasagens permanece **NOT_VERIFIED**.

| OSC | Máx. erro integral L/G (kg) | Máx. erro / inventário inicial L/G | Máx. erro por ciclo L/G (kg) |
|---|---:|---:|---:|
| Refino espacial | 9,2814e-10 / 1,2338e-10 | 1,2913e-11 / 6,9726e-10 | 1,7138e-11 / 3,8858e-13 |
| Refino temporal | 1,5286e-9 / 8,0766e-10 | 2,1268e-11 / 4,5642e-9 | 7,6184e-11 / 3,1125e-13 |
| Extensão | 1,5521e-9 / 1,8010e-10 | 2,1594e-11 / 1,0177e-9 | 7,4735e-11 / 6,6916e-13 |

Na extensão, a normalização por massa processada líquida foi
3,6448e-13/1,6517e-11 no intervalo integral e chegou a
2,7351e-12/9,5629e-12 nos ciclos (líquido/gás). As séries incluem
inventários separados e fluxos de entrada, base e saída. O balanço
integral independente por trecho permanece **NOT_VERIFIED** porque o
solver não acumula fluxos aceitos na base. Os contadores de passo são
completos e a série de `dt_s` é amostrada. Nas séries de 2400 s
sem rejeições, foram observados 19207/19222 passos positivos para
teto 0,125 s e 38403/38418 para teto 0,0625 s. O menor passo positivo
observado foi 0,006192 s, com passos iniciais menores previstos no deck.

Os tempos de parede dos refinamentos OSC foram cerca de 1082 s em
malha e 988 s em passo. A extensão custou 815 s. Considerando os
quatro planos, houve 13 execuções, nove completas e quatro falhas
numéricas. Nenhuma terminou por timeout.

## Status por requisito

| Requisito | Status | Evidência e limite |
|---|---|---|
| Legado, modos fixos, operadores e contornos | **PASS** | Regressão contra `ee5810e` em 1 s, smoke de deriva fixa horizontal/descendente e testes de Jacobianos no mesmo fechamento. |
| Eq. (15) literal, unidades, sinais e limites | **PASS** | Referências independentes, condições monofásicas explícitas e comparação SI nos CSVs. A qualificação física do próprio critério permanece aberta. |
| Transição reversa e comutação em Newton | **BLOCKED** | Falta regra reversível rastreável às fontes. As derivadas correspondem ao fechamento fixo efetivamente usado. |
| Regressão de admissibilidade e junção | **PASS** | Fluxos conservativos nos testes e igualdade do histórico fixo/diagnóstico, inclusive nos snapshots terminais. SS1 continua falhando no traço MATCHED. |
| Balanço integral global | **PASS** em OSC, **FAIL** no controle pelo inventário | O erro do controle é comum aos dois modos. Sem alteração de tolerância. |
| Acomodação numérica de OSC | **PASS** na extensão e nos dois refinamentos | Oscilação com gás positivo na base amostrada. Estudo assintótico **NOT_VERIFIED**. |
| SS1 e SS3 a 2400 s | **FAIL** | Esgotamento da busca linear nos dois casos, com mecanismos de recusa distintos. Período e regime **NOT_VERIFIED**. |
| Validação física da seleção e de severe slugging | **NOT_VERIFIED** | Falta modelo de retorno, convergência e confronto experimental compatível. |

## Artefatos e reprodução

O [índice dos manifestos](../validation/diagnostics/burr15_campaign_20261119/summary.json)
registra SHA, compilador/dependências, entradas, comandos, status,
passos, tempo de parede e hashes. Dados brutos, decks, planos selados,
binários locais e medições estão em
`/home/leonardoab/Documentos/PRIMS-F-burr15-results-20261119`.
As pastas são `v1`, `v2_stable`, `v2_osc` e `v2_osc3600`. Em cada
uma, `runs` guarda manifestos, logs, snapshots e CSVs. `analysis_spatial*`
guarda candidatos de frente, fluxos locais assinados e métricas de fase.
`v1/independent_burr15/v2` inclui conferência da Eq. (15) e
`v1/comparison.json` contém as comparações de identidade.

Execute os comandos de build e testes da [nota de uso](PIPE_REGIME_BURR15.md)
a partir da árvore congelada. Para reproduzir um caso, copie seu deck e
aponte `out_dir` para um destino novo. A linha de comando e o ambiente
originais constam de cada manifesto. Para pós-processar as saídas:

```sh
python3 scripts/analyze_burr15_campaign.py --output /caminho/analise-nova RUN_FIXED
python3 scripts/verify_burr15_observations.py --fields RUN_AUTO/fields.csv \
  --burr15 RUN_AUTO/burr15.csv --diameter 0.05254 --gravity 9.81 \
  --expected-pipe-cells 51 --max-gap 0.5 --output /caminho/burr15-verificado.json
python3 tests/test_pipe_regime_cli.py --binary bin/prims_f --output tmp/b15r
```

O [replay](../validation/diagnostics/burr15_rejection_replay_20261119/)
tem comando, fontes e resumo por proposta. O PR #2 permanece draft e
estas alterações estão somente na branch local.
