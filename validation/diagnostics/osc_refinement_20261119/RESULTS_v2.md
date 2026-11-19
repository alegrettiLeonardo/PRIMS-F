# OSC: refinamentos independentes até N401

Operador `69813ee70551fa4a3f075919aeafc45de9693f11`; binário `7442291fd6cdc311627e79500dc5a4804813b42c7bf4daecf88b4b8b2322288a`. Os manifestos, controles efetivos, fontes, bibliotecas, testes e caminhos dos dados externos estão vinculados por hash em `terminal_index.json`.

| N por trecho | dt teto (s) | Estado | Período 20/80 (s) | Semi-amplitude média (kPa) | CPU do solver (s) | Maior erro/inventário | Balanços e 10 ciclos |
|---:|---:|---|---:|---:|---:|---:|---|
| 51 | 0.125 | COMPLETED | 23.108617598 | 9.909017 | 821.134 | 1.029598e-07 | PASS |
| 101 | 0.125 | COMPLETED | 23.266445936 | 12.770664 | 1793.007 | 3.275331e-08 | PASS |
| 201 | 0.125 | COMPLETED | 23.345494597 | 14.011107 | 3437.280 | 2.235372e-08 | PASS |
| 401 | 0.125 | COMPLETED | 23.385039694 | 14.596904 | 6929.323 | 1.308336e-08 | PASS |
| 51 | 0.0625 | COMPLETED | 23.265643212 | 13.763825 | 1547.012 | 9.703594e-07 | PASS |
| 51 | 0.03125 | COMPLETED | 23.344197706 | 15.347389 | 3022.722 | 2.723117e-07 | PASS |

A última coluna exige conclusão, contrato aritmético P2, seis balanços estritos de inventário/massa processada, sessenta balanços por ciclo e acomodação dos treze canais. Ela não declara convergência nem validação experimental. Os custos de parede usados no índice vêm do solver; o coordenador teve uma reserva de fila documentada.

| Par adjacente | Eixo | Diferença relativa do período | Diferença relativa da semi-amplitude | Critério da semi-amplitude |
|---|---|---:|---:|---|
| N51/dt=0.125 → N101/dt=0.125 | h | 0.678352% | 22.407971% | FAIL |
| N101/dt=0.125 → N201/dt=0.125 | h | 0.338603% | 8.853285% | FAIL |
| N201/dt=0.125 → N401/dt=0.125 | h | 0.169104% | 4.013161% | FAIL |
| N51/dt=0.125 → N51/dt=0.0625 | dt | 0.674925% | 28.006805% | FAIL |
| N51/dt=0.0625 → N51/dt=0.03125 | dt | 0.336505% | 10.318133% | FAIL |

As diferenças da tabela são `abs(fino-grosso)/max(abs(fino),abs(grosso))`. Os testes usam 1% mais os pisos pré-declarados; os valores completos por canal, média, mínimo, máximo e excursão permanecem em `adjacent_orders_v1.json`.

| Três níveis | Eixo | p do período | p da semi-amplitude |
|---|---|---|---|
| N51/dt=0.125 / N101/dt=0.125 / N201/dt=0.125 | h | 1.019137 (POSITIVE_OBSERVED_ORDER) | 1.230025 (POSITIVE_OBSERVED_ORDER) |
| N101/dt=0.125 / N201/dt=0.125 / N401/dt=0.125 | h | 1.010053 (POSITIVE_OBSERVED_ORDER) | 1.093679 (POSITIVE_OBSERVED_ORDER) |
| N51/dt=0.125 / N51/dt=0.0625 / N51/dt=0.03125 | dt | 0.999234 (POSITIVE_OBSERVED_ORDER) | 1.283484 (POSITIVE_OBSERVED_ORDER) |

As razões espaciais usam h=Lp/N verificado nos campos. Ordens são somente observadas; o outro eixo não foi demonstrado resolvido e não há GCI ou prova assintótica. Sinais alternantes, diferenças abaixo do piso e passos efetivos fora do teto permanecem explícitos.

A pressão marcadora e os fluxos de face ficam em locais físicos fixos. Os holdups legados “base” e “saída” vêm de centros móveis com N; suas ordens espaciais não são de um mesmo ponto. O complemento `fixed_trace_comparison.json` compara `1-alpha_riser_trace` na junção fixa, nas mesmas dez bordas, sem alterar os treze critérios originais.

Formas de onda por picos de pressão usam 256 pontos de fase (`adjacent_phase_index.json`); o complemento do traço usa 256 pontos entre os marcadores 20/80 e fica rotulado separadamente. A cadência observada está em `monitor_cadence_final.json`: 256 bins interpolados não equivalem a 256 observações independentes.

Os seis balanços globais e por ciclo, cotas aritméticas, contadores, passos efetivos, dados de cada trecho e controles preservados constam do índice. Nenhum resíduo de Newton foi subtraído do erro físico ou adicionado à tolerância física de 1e-6. Bloqueio total de gás não foi usado como requisito universal de OSC/SS3.

Os resultados são uma qualificação numérica desta coorte. Não demonstram por si só o regime experimental nem autorizam selecionar parâmetros para produzir slugging. Qualquer extensão de duração tem de permanecer em plano e índice separados; esta tabela conserva as integrações originais de 3600 s.
