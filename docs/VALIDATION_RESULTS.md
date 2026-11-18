# Resultados da implementação e campanha de validação — 18/11/2026

**Operadores verificados; OSC apresenta ciclos numéricos localmente acomodados em duas resoluções, mas a captura física e a convergência global de severe slugging ainda não estão aprovadas.**

A campanha contém **28 configurações** com duração solicitada de 2400 s: **13 completas, 13 com falha do solver e 2 interrompidas pela perda da sessão de execução**. Duas tentativas anteriores dessas últimas configurações atingiram o limite de 1800 s de tempo de execução e foram arquivadas separadamente: são **30 processos/tentativas**, não 30 configurações físicas. Os retries tinham limite de 7200 s, mas a sessão foi perdida antes da conclusão. Não se atribui essa interrupção ao solver nem se fabrica um run_status de sucesso ou falha física.

A produção permanece na versão com **6/6 programas e 2185 verificações Fortran**, além dos **16 testes do pós-processamento**. A correção de busca da cauda de alfa é um candidato separado, revisado localmente e **não aplicado** ao código ou às execuções desta campanha. Fontes e executável são identificados por SHA-256 em cada manifesto.

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
| D_OSC | 51+51 | 0.5 | ISOTHERMAL | 382,980 | falha do solver | 7,03e-15 | 1,09e-12 |
| D_SS3 | 51+51 | 0.5 | BAROTROPIC | 97,966 | falha do solver | 1,20e-15 | 1,32e-13 |
| D_SS3 | 101+101 | 0.5 | BAROTROPIC | 104,874 | falha do solver | 2,24e-15 | 3,81e-14 |
| D_SS3 | 201+201 | 0.5 | BAROTROPIC | 109,508 | falha do solver | 1,18e-15 | 1,91e-14 |
| D_SS3 | 51+51 | 0.25 | BAROTROPIC | 97,193 | falha do solver | 1,33e-15 | 3,16e-14 |
| D_SS3 | 51+51 | 0.125 | BAROTROPIC | 97,132 | falha do solver | 1,09e-10 | 4,35e-11 |
| D_SS3 | 51+51 | 0.5 | ISOTHERMAL | 72,632 | falha do solver | 9,02e-16 | 1,05e-14 |
| D_SS1 | 51+51 | 0.5 | BAROTROPIC | 77,814 | falha do solver | 2,39e-15 | 2,89e-15 |
| D_SS1 | 101+101 | 0.5 | BAROTROPIC | 80,718 | falha do solver | 3,23e-15 | 1,96e-15 |
| D_SS1 | 201+201 | 0.5 | BAROTROPIC | 83,835 | falha do solver | 2,37e-15 | 1,33e-15 |
| D_SS1 | 51+51 | 0.25 | BAROTROPIC | 78,317 | falha do solver | 4,91e-15 | 1,14e-15 |
| D_SS1 | 51+51 | 0.125 | BAROTROPIC | 78,572 | falha do solver | 6,24e-15 | 3,74e-15 |
| D_SS1 | 51+51 | 0.5 | ISOTHERMAL | 50,281 | falha do solver | 1,23e-15 | 2,71e-15 |
| D_OSC | 51+51 | 0.0625 | BAROTROPIC | 2400,000 | completa | 5,45e-13 | 1,11e-10 |
| D_OSC | 51+51 | 0.03125 | BAROTROPIC | 1329,811 | ambiente interrompido | 1,73e-12 | 5,86e-10 |
| D_OSC | 101+101 | 0.125 | BAROTROPIC | 2400,000 | completa | 3,28e-13 | 1,70e-11 |
| D_OSC | 201+201 | 0.125 | BAROTROPIC | 1351,557 | ambiente interrompido | 1,75e-13 | 2,25e-11 |

O erro relativo usa o maior valor entre inventário inicial e massa acumulada de entrada. O critério é erro absoluto ≤ 1e-8 kg + 1e-6 vezes essa escala; o pacote conserva também o erro absoluto. Maiores resíduos absolutos: líquido 5,08e-04 kg e gás 6,54e-07 kg; maiores relativos 1,65e-07 e 1,51e-07, respectivamente. Conservação até a parada não transforma uma execução incompleta em validação.

## OSC: crescimento inicial e ciclos finais

D_Stable permanece globalmente amortecido. D_OSC barotrópico amortece com teto de 0,5 s nas três malhas e com 0,25 s em 51+51, mas oscila ao reduzir o passo. O envelope usa máximos e mínimos refinados localmente por uma parábola de três amostras; a semi-amplitude é metade da diferença entre o pico e a linha de base interpolada entre os dois mínimos vizinhos. Os ajustes abaixo usam somente 40–400 s e semi-amplitudes ≥ 1 Pa. São taxas do sinal, não autovalores demonstrados do operador.

| N por trecho | dt máximo (s) | Situação | Taxa inicial de log-amplitude (s⁻¹) | R² |
|---:|---:|---|---:|---:|
| 51 | 0.5 | completa | -0,010982774 | 1,000000000 |
| 51 | 0.25 | completa | -0,002002954 | 1,000000000 |
| 51 | 0.125 | completa | 0,002513812 | 0,999999911 |
| 51 | 0.0625 | completa | 0,004769835 | 0,999999408 |
| 51 | 0.03125 | parcial; interrupção do ambiente | 0,005894235 | 0,999998428 |
| 101 | 0.125 | completa | 0,004150573 | 0,999999776 |
| 201 | 0.125 | parcial; interrupção do ambiente | 0,004985521 | 0,999999997 |

O crescimento aumenta quando se reduz o passo. Para o problema escalar z′=(σ+iω)z, Euler implícito tem fator 1/[1−h(σ+iω)] e pode amortecer um modo contínuo instável se o passo for grande. Um ajuste descritivo desses envelopes dá mudança de sinal perto de 0,194535 s. O modo LAGGED separa a dissipação congelada dos fluxos físicos implícitos; esse ajuste escalar não demonstra o espectro do solver. A ordem observada da taxa inicial nos passos 0,125/0,0625/0,03125 s é 1,004627, válida somente para essa grandeza e janela, sem alegação de GCI global.

### Dez últimos ciclos completos delimitados pelo detector

Onze picos consecutivos, cada um com os dois mínimos vizinhos observados, definem dez intervalos. As três execuções abaixo alcançaram 2400 s; as janelas explícitas evitam confundir o último pico delimitado com o instante final exato do processo. A amplitude é semi-amplitude, não excursão pico a pico.

| N por trecho | dt (s) | Janela (s) | Período médio (s) | Semi-amplitude média (kPa) | Variação de amplitude |
|---:|---:|---|---:|---:|---:|
| 51 | 0.125 | 2154,819–2385,849 | 23,103011 | 9,699735 | 2,398601% |
| 51 | 0.0625 | 2145,556–2378,212 | 23,265630 | 13,763806 | 0,000523% |
| 101 | 0.125 | 2144,055–2376,719 | 23,266396 | 12,769678 | 0,021283% |

Os dois refinamentos completos, 51+51/0,0625 s e 101+101/0,125 s, dão evidência positiva de **acomodação local dos ciclos numéricos**: as amplitudes variam apenas 0,000523% e 0,021283% nesses dez ciclos. As semi-amplitudes médias ainda diferem entre discretizações. Comparando cada refinamento ao caso 51+51/0,125 s, muda-se respectivamente somente o teto de passo ou somente a malha. Comparar os dois novos casos diretamente muda ambos os controles, portanto não é um par independente de convergência. A estabilidade local não prova ciclo limite assintótico, convergência global ou concordância experimental.

O fluxo de gás na base permanece positivo nos instantes observados dos três casos oscilatórios completos. O detector não encontra a sequência de bloqueio forte/liberação; sua ausência não exclui OSC/SS3. OSC não foi confundido com SS1. Período, amplitude, inventários, fluxos de produção e sequência física precisam ser confrontados conjuntamente.

## Interrupções do ambiente e tentativas arquivadas

- OSC, 201+201 células, dt máximo 0.125 s: último tempo observado 1351,556852760 s, manifesto `INTERRUPTED`, sem `run_status.txt` conclusivo. Os dados disponíveis são parciais e não cobrem 900–2400 s.
- OSC, 51+51 células, dt máximo 0.03125 s: último tempo observado 1329,810990000 s, manifesto `INTERRUPTED`, sem `run_status.txt` conclusivo. Os dados disponíveis são parciais e não cobrem 900–2400 s.

Os dois processos anteriores dessas configurações estão em `validation/runtime_limited`, identificados como `TIMEOUT` com limite de 1800 s de execução. Os retries e os registros anteriores não foram concatenados para fabricar uma trajetória completa. Ambas as configurações ainda precisam ser repetidas até 2400 s em ambiente persistente.

## Limites distintos dos casos SS1, SS3 e OSC isotérmico

Em SS3, nas malhas 51/101 com dt máximo0,5 s, a continuação local examinada alcança ul=0 e a direção de Newton exige ul<0 com gás ascendente, fora do domínio qualificado do fechamento ID3. A malha201 atinge antes a limitação da junção; o caso51/dt0,25 falha em um ponto da quadratura de Roe avaliado na inclinação média, com células ainda admissíveis. Em SS1 o traço de gás calculado pelo acoplamento MATCHED/Roe pode tornar-se inadmissível com as células ainda válidas; isso não prova retorno físico do gás. Ver `TRANSIENT_DOMAIN_LIMITS.md` e os probes reproduzíveis para os estados e versões examinados.

O caso OSC isotérmico falhou em 382,980199622 s na junção, com células ainda admissíveis; é distinto dos OSC barotrópicos completos e das duas interrupções do ambiente. O status 3021 do SS1/0,25 s também tem uma limitação numérica de busca: uma raiz positiva pode existir abaixo do piso da varredura. Não foi tratado automaticamente como prova de inexistência física.

## Correção candidata da cauda de alfa: ainda não aplicada

O contraprobe sintético encontra alpha=3,61178969909e-13 para fluxo gasoso positivo, mas a varredura de produção termina em aproximadamente 6,9144e-13 e retorna NO_ROOT. O candidato em `validation/diagnostics/alpha_tail_candidate` acrescenta uma busca geométrica seguida de bisseção relativa somente após a busca antiga falhar, no domínio cocorrente ascendente ID3. A monotonicidade da cauda é verificada; não há piso substituído pela raiz nem extensão para fluxos negativos. Passou em 55 guards O0/O3 e na revisão independente de 129 contratos, preservando sucessos anteriores. Esses números não foram somados aos 2185 testes da produção. Aplicar, repetir a suíte e repetir os casos interrompidos pelo solver são tarefas futuras; o candidato não resolve automaticamente ul<0 ou a quadratura de Roe fora do domínio.

## Referências experimentais e divergências confirmadas

- A dissertação [D], tabela 6.7, p.138/PDF140, declara OSC experimental 33,9 s e simulado 3,3 s; 3,3 s não é alvo experimental. Baliño–Burr–Nemoto2010 [B], tabela 3, p.653/PDF11, caso 13, fornece faixa experimental modulada **19,3–33,8 s**. O período numérico tardio de aproximadamente 23,27 s está dentro dessa faixa, mas essa compatibilidade isolada não valida o modelo.
- Para SS3, a página primária de [B] foi agora conferida visualmente: o caso 11 imprime **Texp=76 s, Tsim=72,2 s e erro −1%**. O erro recalculado é −5%; há uma inconsistência na própria tabela. A dissertação imprime **26,0 s**. Trata-se de divergência documental real; o valor literal D não foi substituído silenciosamente.
- O caso SS1 externo próximo informa 85,6 s, mas a temperatura e a vazão líquida diferem ligeiramente do perfil D. Não foram obtidos sinais experimentais brutos de produção no topo.

## Reprodução

```sh
make -j4 FC=gfortran LIBS='-llapack -lblas'
make test FC=gfortran LIBS='-llapack -lblas'
python3 tests/test_validation_evidence.py
python3 scripts/run_study.py --workers 4 --time 2400 --label final
python3 scripts/run_study.py --study osc-refinement --workers 4 --time 2400 --label refined
python3 scripts/analyze_validation.py --pattern '[fr]*' --output validation/analysis_final --plots
python3 scripts/compare_validation.py --pattern '[fr]*' --metrics validation/analysis_final/metrics.json --output validation/analysis_final
python3 scripts/analyze_phase_sequence.py --pattern '[fr]*' --metrics validation/analysis_final/metrics.json --output validation/analysis_final
python3 scripts/analyze_osc_envelopes.py --output validation/analysis_final
```

As tabelas/JSON de comparação e sequência estão em `validation/analysis_final`. Os dados brutos, decks, hashes e logs são fornecidos no pacote de evidências da implementação. A especificação das fontes está em `BENCHMARK_REFERENCE.md` e `benchmark_sources.json`.
