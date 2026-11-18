# Resultados da implementação e campanha de validação — 18/11/2026

**Operadores verificados; captura de ciclos completos de severe slugging não aprovada.**

Foram executadas 24 tentativas de 2400 s: 11 chegaram ao tempo solicitado. A rodada final de testes isolados/integrados passou em 6/6 programas, com 2185 verificações. A campanha mantém fontes e executável identificados por SHA-256 em cada manifesto.

O perfil principal é a catenária da dissertação com EOS gasosa p/343², entrada de gás referida a condições padrão, Bhagwat conforme Azevedo–Baliño–Burr2015 no riser e fechamento estratificado no oleoduto. MATCHED, Roe com dissipação LAGGED, FV1 e Euler implícito são as opções numéricas efetivas. Os pares ISOTHERMAL alteram somente a lei de densidade gasosa.

## Matriz efetivamente executada

| Caso | Malha | dt máximo (s) | EOS | Tempo atingido (s) | Resultado | Erro massa L relativo | Erro massa G relativo |
|---|---:|---:|---|---:|---|---:|---:|
| D_Stable | 51+51 | 0.5 | BAROTROPIC | 2400,000 | completa | 1,23e-08 | 4,00e-08 |
| D_Stable | 101+101 | 0.5 | BAROTROPIC | 2400,000 | completa | 1,69e-08 | 4,74e-08 |
| D_Stable | 201+201 | 0.5 | BAROTROPIC | 2400,000 | completa | 8,02e-09 | 3,39e-08 |
| D_Stable | 51+51 | 0.25 | BAROTROPIC | 2400,000 | completa | 7,24e-08 | 5,50e-08 |
| D_Stable | 51+51 | 0.125 | BAROTROPIC | 2400,000 | completa | 1,65e-07 | 1,51e-07 |
| D_Stable | 51+51 | 0.5 | ISOTHERMAL | 2400,000 | completa | 1,10e-08 | 4,54e-08 |
| D_OSC | 51+51 | 0.5 | BAROTROPIC | 2400,000 | completa | 9,35e-08 | 1,84e-08 |
| D_OSC | 101+101 | 0.5 | BAROTROPIC | 2400,000 | completa | 8,80e-08 | 1,87e-08 |
| D_OSC | 201+201 | 0.5 | BAROTROPIC | 2400,000 | completa | 8,82e-08 | 3,49e-09 |
| D_OSC | 51+51 | 0.25 | BAROTROPIC | 2400,000 | completa | 3,88e-13 | 1,23e-10 |
| D_OSC | 51+51 | 0.125 | BAROTROPIC | 2400,000 | completa | 3,45e-13 | 2,48e-11 |
| D_OSC | 51+51 | 0.5 | ISOTHERMAL | 382,980 | incompleta | 7,03e-15 | 1,09e-12 |
| D_SS3 | 51+51 | 0.5 | BAROTROPIC | 97,966 | incompleta | 1,20e-15 | 1,32e-13 |
| D_SS3 | 101+101 | 0.5 | BAROTROPIC | 104,874 | incompleta | 2,24e-15 | 3,81e-14 |
| D_SS3 | 201+201 | 0.5 | BAROTROPIC | 109,508 | incompleta | 1,18e-15 | 1,91e-14 |
| D_SS3 | 51+51 | 0.25 | BAROTROPIC | 97,193 | incompleta | 1,33e-15 | 3,16e-14 |
| D_SS3 | 51+51 | 0.125 | BAROTROPIC | 97,132 | incompleta | 1,09e-10 | 4,35e-11 |
| D_SS3 | 51+51 | 0.5 | ISOTHERMAL | 72,632 | incompleta | 9,02e-16 | 1,05e-14 |
| D_SS1 | 51+51 | 0.5 | BAROTROPIC | 77,814 | incompleta | 2,39e-15 | 2,89e-15 |
| D_SS1 | 101+101 | 0.5 | BAROTROPIC | 80,718 | incompleta | 3,23e-15 | 1,96e-15 |
| D_SS1 | 201+201 | 0.5 | BAROTROPIC | 83,835 | incompleta | 2,37e-15 | 1,33e-15 |
| D_SS1 | 51+51 | 0.25 | BAROTROPIC | 78,317 | incompleta | 4,91e-15 | 1,14e-15 |
| D_SS1 | 51+51 | 0.125 | BAROTROPIC | 78,572 | incompleta | 6,24e-15 | 3,74e-15 |
| D_SS1 | 51+51 | 0.5 | ISOTHERMAL | 50,281 | incompleta | 1,23e-15 | 2,71e-15 |

O erro relativo usa o maior valor entre inventário inicial e massa acumulada de entrada. O critério é erro absoluto ≤ 1e-8 kg + 1e-6 vezes essa escala; o pacote conserva também o erro absoluto. Conservação até a parada não transforma uma execução incompleta em validação.

## Interpretação

D_Stable permanece globalmente amortecido. D_OSC barotrópico amortece com dt máximo0,5s nas três malhas e com dt máximo0,25s em51+51, mas apresenta64intervalos oscilatórios com dt máximo0,125s em51+51: período mediano23,0s, amplitude pico a pico19597,4Pa e amplitude ainda variável entre janelas. Portanto o solver obtém uma solução numérica oscilatória com amplitude crescente, mas a dinâmica ainda não está convergida em passo. Os casos SS1 e SS3 não chegam a dez ciclos úteis. A sensibilidade isotérmica deve ser lida separadamente, sem substituir silenciosamente a EOS declarada na dissertação.

Em SS3, nas malhas 51/101 com dt máximo0,5 s, a continuação local examinada alcança ul=0 e a direção de Newton exige ul<0 com gás ascendente, fora do domínio qualificado do fechamento ID3. A malha201 atinge antes a limitação da junção; o caso51/dt0,25 falha em um ponto da quadratura de Roe avaliado na inclinação média, com células ainda admissíveis. Em SS1 o traço de gás calculado pelo acoplamento MATCHED/Roe pode tornar-se inadmissível com as células ainda válidas; isso não prova retorno físico do gás. Ver `TRANSIENT_DOMAIN_LIMITS.md` e os probes reproduzíveis para os estados e versões examinados.

O refinamento quantifica médias, a forte sensibilidade da dinâmica OSC e o instante de interrupção dos casos incompletos. Ele não demonstra convergência de amplitude/período. OSC0,125s mantém gás atravessando a base; a ausência de bloqueio no detector não excluiOSC/SS3. A identificação de slug severo exige interpretar essa sequência com o regime e confirmar convergência e concordância experimental.

## Reprodução

```sh
make -j4 FC=gfortran LIBS='-llapack -lblas'
make test FC=gfortran LIBS='-llapack -lblas'
python3 scripts/run_study.py --workers 4 --time 2400 --label final
python3 scripts/run_study.py --study osc-refinement --workers 4 --time 2400 --label refined
python3 scripts/analyze_validation.py --pattern '[fr]*' --output validation/analysis_final --plots
```

As tabelas/JSON de comparação e sequência estão em `validation/analysis_final`. Os dados brutos, decks, hashes e logs são fornecidos no pacote de evidências da implementação. A especificação das fontes está em `BENCHMARK_REFERENCE.md` e `benchmark_sources.json`.
