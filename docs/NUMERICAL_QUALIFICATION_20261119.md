# PRIMS-F — qualificação numérica e comparação experimental

Continuação local de `9f2ed21668269bc60c8dbeaa5853b5700a081ac2` em
`audit/numerics-20261119`. Operador congelado para a campanha
`69813ee70551fa4a3f075919aeafc45de9693f11`, GNU Fortran 13.3.0, Fortran
2018, release `-O2`, LAPACK/BLAS. SHA-256 do binário
`7442291fd6cdc311627e79500dc5a4804813b42c7bf4daecf88b4b8b2322288a`.

## Defeitos, correções e verificações

| Problema | Causa demonstrada | Correção e teste | Status |
|---|---|---|---|
| Pressão inválida na saída Chexal | A recuperação interior deixava `status=0` antes de validar a pressão prescrita. Pressão inválida ainda podia acionar trap em debug. | Status separado e guardas de finitude antes das comparações. Testes com pressão zero/NaN, pressão interior abaixo do domínio EOS, fluxo nulo e fase ausente. | **PASS** |
| Roe sem histórico | A ausência de velocidades era convertida em histórico artificial de duas velocidades nulas. | Primeiro nó recupera sem hints. O teste contracorrente vertical distingue AUTO ambíguo de HIGH/LOW e de histórico específico. | **PASS** |
| Contornos nos snapshots SS1/SS3 | Espectros reais no interior e nos estados impostos. Contagem de duas características entrantes na entrada e uma na saída. | Auditoria de sinais e pressão absoluta. Testes de fluxo constante, equilíbrio hidrostático, uniformidade e Jacobiano da montagem. | **PASS** no escopo ensaiado |
| SS1, face 53 | O termo de dissipação de Roe cancela o transporte convectivo de gás e a proposta passa a ter fluxo gasoso negativo no traço. As células vizinhas têm gás presente e velocidades positivas. | Registro online de cada proposta recusada e replay do Jacobiano. O matching característico linear testado também produz massa gasosa negativa. | **BLOCKED** |
| SS3, célula 69 | A linearização chegou ao limite cocorrente da lei Bhagwat. Um limite inferior positivo exclui todas as raízes das propostas amortecidas no intervalo admissível. | `NO_PHYSICAL_ROOT` identificado separadamente de `ROOT_SEARCH_FAILURE`, `INADMISSIBLE_ROOT` e da falha de globalização. Cálculo independente com 80 dígitos e regressão em dois snapshots. | **PASS** em classificação, transiente **BLOCKED** |
| Perda acumulada de massa no controle | Contribuições não nulas dos resíduos de massa dos passos aceitos acumulam-se no balanço. | Integração assinada dos fluxos usados na montagem, com os dois lados da junção e inventários por trecho. Reexecutado com duas tolerâncias de resíduo. | Medições na seção de conservação |

[Contornos e características](../validation/diagnostics/numerics_p0p1_20261119/README.md)
registram os espectros, as restrições impostas, o confronto com as
Eqs. (39)–(42) de Baliño2010 e as hipóteses que separam o modelo FV do
modelo de características do artigo. A diferença entre vazão líquida
volumétrica prescrita no artigo e vazão mássica fixa no solver foi
mantida visível. A saída tem pressão de separador absoluta fixa.
A pressão do separador e as vazões flutuavam na instalação experimental.

A junção mantém a pressão de traço comum e o fluxo por fase com os sinais
opostos nos dois volumes, sem impor igualdade de holdup. As regressões
de fluxo/Jacobiano comparam os operadores no mesmo fechamento e na
mesma configuração de dissipação. O teste hidrostático verifica o
gradiente e os fluxos de repouso com fase ausente. A integração FV
monofásica com Jacobiano reduzido continua sem qualificação.

### Telemetria de SS1 e SS3

Foram feitas novas tentativas com 2400 s solicitados, a partir do início
e com a física e os controles efetivos da série Burr15. Ambas terminaram
por `NEWTON_GLOBALIZATION_FAILURE`. O campo `proposal_status` identifica
a recusa anterior no operador. O resíduo de uma proposta que não chega
a ser montada é marcado como indisponível.

| Caso | Último tempo aceito | Aceitos / rejeitados | Última iteração / backtrack | Primeira recusa na proposta terminal |
|---|---:|---:|---:|---|
| SS1 | 78,571788477 s | 662 / 24 | 13 / 18 | `3012`, face 53, traço do riser |
| SS3 | 97,131852339 s | 809 / 24 | 8 / 18 | `1017`, célula 69 |

`rejected_proposals.csv` guarda tempo, passo tentado, iteração,
\(\lambda\), a norma da linearização e o local de cada proposta
recusada. A última proposta de ambos usou \(\lambda=2^{-17}\).
`replay_69813ee` guarda os estados conservativos das propostas, as
primitivas/fluxos observáveis, os blocos do Jacobiano e os resíduos
por célula. Consulte o [limite de SS3](SS3_COCURRENT_BOUND_20261119.md)
para a distinção entre ponto linearizado admissível e proposta sem raiz.

Uma direção alternativa de Newton, com região de confiança ou restrição
ativa, pode conservar admissibilidade dentro do modelo. A evidência
exclui todas as propostas testadas nessa direção. A continuação física
de SS1 precisa representar a frente, o inventário bloqueado e a
compressão descritos nas Eqs. (1)–(14), pp. 645–646 de Baliño2010.
No SS3, continuar além da fronteira cocorrente requer sustentar a nova
lei de contrafluxo inclinado. A célula 69 está inclinada 0,901239 rad
desde a horizontal. O domínio contracorrente vertical de Chexal1992
não cobre essa geometria.

## Eq. (15) e dados primários

O detector conserva a forma literal da Eq. (15), com fator quadrático
\((1-h_L/D)^2\) e densidade \(\rho_l\). As Eqs. (53)–(54) levam a
\(\rho_g\) no denominador. Em dois envelopes de medições estratificadas
de Carraretto et al. (2020), a margem dimensional mínima foi positiva
em 0,346475 e 0,446875 m/s mesmo propagando a incerteza de leitura e de
vazão. O detector literal discordou do regime observado nos dois
envelopes. Isso constitui **FAIL** do teste de classificação local
nesses estados. A comutação e uma regra reversa seguem **BLOCKED**.
[Fontes, instrumentação e cálculo independente](P4_P5_PRIMARY_EVIDENCE_20261119.md)

O estado não estratificado publicado no caso 14 de Baliño não fornece
uma altura local de interface plana. A guarda de configuração recusa
AUTO com fechamento drift. A classificação de um evento de perda
efetivamente observado permanece **BLOCKED** pela ausência, nas fontes
conferidas, de holdup, velocidade local e rótulo simultâneos de transição.

## Conservação e tolerância de resíduo

O par estável usa os mesmos 51+51 volumes, passo inicial 0,001 s,
teto 0,125 s, EOS, geometria, perturbação e contornos. Além do diretório de
saída, só `tol_res` difere. Ambos completaram 2400 s com 19.222 passos
aceitos e zero rejeições. A norma mais apertada é uma condição numérica
qualificada para esse controle, e o deck histórico continua rastreável.

Nas tabelas, `E` é o maior erro integral absoluto até 2400 s.
`M0` é o inventário inicial e `Mproc` soma as
massas absolutas processadas na entrada e na saída da seção considerada.

| Fase / seção | `1e-8`: E (kg) | E/M₀ | E/Mproc | `1e-10`: E (kg) | E/M₀ | E/Mproc |
|---|---:|---:|---:|---:|---:|---:|
| Líquido, sistema | 5,0783e-4 | 6,1911e-6 | 8,2395e-8 | 7,8086e-7 | 9,5197e-9 | 1,2670e-10 |
| Líquido, oleoduto | 3,0565e-6 | 5,1383e-8 | 4,9592e-10 | 1,4620e-9 | 2,4578e-11 | 2,3721e-13 |
| Líquido, riser | 5,0477e-4 | 2,2394e-5 | 8,1900e-8 | 7,8232e-7 | 3,4707e-8 | 1,2693e-10 |
| Gás, sistema | 6,5445e-7 | 2,2332e-6 | 7,5697e-8 | 6,0646e-9 | 2,0695e-8 | 7,0146e-10 |
| Gás, oleoduto | 1,2382e-6 | 4,5048e-6 | 1,4322e-7 | 1,7416e-9 | 6,3362e-9 | 2,0144e-10 |
| Gás, riser | 1,8927e-6 | 1,0407e-4 | 2,1892e-7 | 4,3230e-9 | 2,3771e-7 | 5,0002e-10 |

O balanço pelo inventário inicial fica **FAIL** com `1e-8` e
**PASS** em todos os seis canais com `1e-10`. A soma ponderada dos
resíduos de massa aceitos acompanha a deriva do balanço. Para líquido
global com `1e-8`, a soma termina em 5,07830e-4 kg, ante erro de
5,07829e-4 kg. Os fluxos internos cancelam nos integrais, com discrepância
zero entre as duas faces da junção em ambas as fases.

| tol_res | Iterações Newton | Solves lineares | Avaliações de resíduo, incluindo Jacobiano | CPU (s) | Parede do solver (s) |
|---|---:|---:|---:|---:|---:|
| 1e-8 | 19452 | 2677 | 22129 | 258,390 | 258,946 |
| 1e-10 | 19934 | 4323 | 24257 | 269,555 | 270,143 |

O custo de CPU aumentou 4,32%. Ambos ficaram estacionários e as
estatísticas de pressão, holdup, fluxos e inventários passaram a
comparação de 1% com os pisos de resolução documentados.
O auditor aritmético inicial, com margem fixa, marcou discrepâncias
da ordem de 3,45e-14 kg após milhares de somas. A segunda revisão usa
o limite de somatório dependente de `n` e da soma das magnitudes.
A margem aritmética verifica a reprodução do acumulador, e o critério
de conservação de `1e-6` segue aplicado ao erro medido.
Os dois conjuntos de análises foram preservados.
[Instrumentação, sinais e verificador](P2_ACCEPTED_TELEMETRY.md)

Até os términos, o maior `E/M0` por seção/fase foi 1,558e-9 em
SS1 e 2,649e-9 em SS3. O resultado é **PASS** para conservação
das trajetórias parciais. A caracterização dos regimes segue **BLOCKED**.

## Convergência

O refinamento até N401 e a campanha de passo adicional foram condicionados
à estabilização de P0–P2. P1 continua **BLOCKED** em SS1/SS3. A extensão
OSC em malha 51+51 mede acomodação em uma discretização específica.
Os volumes são uniformes em comprimento de arco dentro de cada trecho,
com cerca de 1,1255 m no oleoduto e 0,2460 m no riser.
A série Burr15 anterior conserva sua proveniência e fornece evidência
histórica de sensibilidade, sem integrar a nova coorte. A região
assintótica e o alvo de 1% para amplitude e período permanecem
**NOT_VERIFIED**.

## OSC calculado e confronto com Baliño, Burr e Nemoto

OSC completou 3600 s, 28.822 passos aceitos, zero rejeições e 821,134 s
de CPU. Onze picos da pressão na primeira célula do riser, cercados por
mínimos, delimitam dez ciclos de 3356,386 a 3587,473 s. O período médio
é 23,108640 s, a semi-amplitude média é 9915,43 Pa e sua variação entre
os picos extremos é +0,00641%. O protocolo independente por cruzamentos
20/80 encontrou 23,108618 s e acomodação dos 13 canais.

O mínimo de fluxo gasoso efetivamente aceito na junção foi
1,75792e-3 kg/s. O maior balanço normalizado pelo inventário inicial
entre as seis combinações de fase e trecho foi 1,0296e-7, em gás no
riser. Os máximos globais por ciclo foram 7,4735e-11 kg no líquido e
6,6916e-13 kg no gás. O erro integral e a normalização por massa
processada estão no relatório do ledger.

Os campos variam entre frações de vazio 0,2652 e 0,5512. As isolinhas
de holdup na figura são calculadas dentro de cada seção. Os
cruzamentos de níveis de fração de vazio 0,4/0,5/0,6 são observações
geométricas, com 10448/3176/0 candidatos nos campos disponíveis.
Entre 900 e 3600 s, pressão, inventários e fluxos deram 118
episódios compatíveis com acumulação/compressão, 101 excursões de
produção líquida e 117 episódios compatíveis com blowdown.
A penetração é compatível com gás atravessando a base em todos
os passos aceitos. A correspondência de frentes materiais e
eventos com o ensaio permanece **NOT_VERIFIED**.

| Referência experimental | Observável publicado | Resultado PRIMS-F | Comparabilidade |
|---|---|---|---|
| Caso 2, SS1 | `Texp=85,6 s`, `DeltaPexp=0,51 bar` | Truncado em 78,57 s na junção | Período/amplitude **BLOCKED** |
| Caso 11, SS3 | `Texp=76 s`, `DeltaPexp=0,86 bar` | Truncado em 97,13 s na inversão cocorrente | Período/amplitude **BLOCKED** |
| Caso 13, OSC | `Texp=19,3–33,8 s`, `DeltaPexp=Max 0,69 bar` | `T=23,10864 s`, semi-amplitude 0,099154 bar | Período cai no intervalo publicado. Os escalares de amplitude têm definições distintas. |
| Caso 13, Fig. 18(a), janela 1200–1500 s | Excursão central digitalizada 0,6453 bar, banda 0,5363–0,7598 bar | Excursão 0,1982 bar na junção, janela final de cerca de 300 s | Diferença aproximada −0,4471 bar frente ao centro digitalizado. A excursão calculada fica abaixo de toda a banda gráfica. |

As janelas da sobreposição têm origens temporais independentes. A
Fig. 18(a) fornece uma curva digitalizada de pressão absoluta com
544 colunas, sem acesso aos registros originais de 10 Hz. O espectro
comparativo usa uma janela de Hann e potência relativa por série.
A curva digitalizada tem picos em bins aproximados de 0,02662 e
0,05324 Hz. O cálculo tem pico dominante em 0,04332 Hz, com
harmônico em 0,08665 Hz. Os bins da digitalização não substituem os
períodos publicados. O padrão espectral e a excursão observada
permanecem discrepantes apesar da proximidade de escalas temporais.

A pressão do traço da junção e a pressão da primeira célula aparecem
separadamente no gráfico. Os ensaios registraram flutuações dos
contornos, e o cálculo usa uma pressão de separador constante e
vazões impostas. A instalação inclui volume de ar adicional e
instrumentação sem séries sincronizadas disponíveis para esta
comparação. As diferenças de propriedades, EOS e entrada do caso 2
estão identificadas na [nota experimental](P4_P5_PRIMARY_EVIDENCE_20261119.md).
Desempenho experimental para OSC e reprodução física de severe slugging
nos três casos permanecem **NOT_VERIFIED**.

## Execução, hashes e reprodução

As cinco execuções tiveram três términos normais e duas falhas numéricas
SS1/SS3. O [índice](../validation/diagnostics/numerics_campaign_20261119/index.json)
contém SHA, flags do verificador, dependências, comandos, manifestos por
execução e hashes dos fluxos aceitos, propostas, campos, medições e
figuras. O release foi gerado no worktree detached
`PRIMS-F-numerics-run-20261119` e os resultados externos estão em
`/home/leonardoab/Documentos/PRIMS-F-numerics-results-20261119`.
Use as pastas `v1_2400`, `v1_osc3600`, `analysis`,
`osc_experimental_final`, `osc_field_analysis` e `replay_69813ee`.
As séries originais e as análises do auditor anterior foram preservadas.

O teste Fortran no SHA congelado passou **3.625 verificações em 13
programas**. Passaram também **58 testes Python**, **18 verificações
de CLI** no binário com traps e **cinco smoke tests da instrumentação**
no release congelado. Builds debug e release foram executados.

```sh
make -j4 FC=gfortran LIBS='-llapack -lblas'
make check_all FC=gfortran LIBS='-llapack -lblas'
make test FC=gfortran LIBS='-llapack -lblas'
PYTHONPATH=scripts python3 -m unittest discover -s tests -p 'test_*.py'
python3 tests/test_pipe_regime_cli.py --binary bin/prims_f_check_all --output tmp/b15n
python3 tests/test_accepted_telemetry_cli.py --binary bin/prims_f --output tmp/telemetry_nova
python3 scripts/analyze_accepted_conservation.py --run-dir RUN --output NOVO_BALANCO.json
python3 scripts/analyze_burr15_campaign.py --output NOVA_ANALISE OSC_RUN
python3 scripts/plot_osc_experimental.py --run OSC_RUN --digitization FIG18_CSV --output NOVAS_FIGURAS
```

O [replay](SS3_COCURRENT_BOUND_20261119.md) e a
[digitalização](P4_P5_PRIMARY_EVIDENCE_20261119.md) têm comandos separados.
Para reproduzir um caso, use o deck efetivo registrado e um destino novo
em `out_dir`. O comando e o ambiente originais estão no manifesto.
Os executáveis e os PDFs de referência estão fora do Git.
