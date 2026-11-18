# PRIMS-F — Fortran 2018

Solver unidimensional de escoamento gás–líquido em oleoduto descendente e riser
em catenária. Esta revisão implementa operadores, derivadas, inicialização,
contornos e testes com referência à dissertação de Leonardo Alegretti Belaroza.

**Estado da qualificação:** os operadores passam nos testes isolados e de
integração. A captura de ciclos completos de *severe slugging* ainda não está
validada. OSC apresenta uma solução numérica oscilatória quando o teto de passo
é reduzido, com fluxo gasoso positivo observado na base; a resposta dinâmica
ainda depende da discretização. Há também casos amortecidos e casos que atingem
o limite de domínio do fechamento ou da junção. Um processo terminado com
código diferente de zero é uma execução incompleta, mesmo que produza CSV.
Consulte `docs/VALIDATION_RESULTS.md` para os resultados da campanha e
`docs/TRANSIENT_DOMAIN_LIMITS.md` para os diagnósticos de SS1 e SS3.

Para continuar a auditoria no WSL, comece por `docs/CODEX_HANDOFF_STATE.md`.
O prompt pronto para o Codex está em `docs/CODEX_CONTINUATION_PROMPT_PT.md`,
e o domínio documentado do futuro fechamento em `docs/CHEXAL_1992_SCOPE.md`.
Chexal ainda não está implementado nesta entrega.

## Compilação e testes

Dependências do executável: GNU Fortran com suporte a Fortran 2018,
BLAS/LAPACK e runtime OpenMP. A implementação foi executada com GNU Fortran
14.2.0. O Makefile mantém diretórios separados para os builds release e debug.

```sh
make -j4 FC=gfortran LIBS='-llapack -lblas'
make test FC=gfortran LIBS='-llapack -lblas'
make check_all FC=gfortran LIBS='-llapack -lblas'
```

O teste recompila os módulos em um diretório próprio com `-std=f2018`,
verificação de limites, traps para operações inválidas/divisão por zero/overflow
e pilha não executável. O resumo JSON contém hashes, comandos e logs.
Não é necessário habilitar trampolins de procedimentos internos.

As seis suítes verificam geometria, fechamentos, primitivas/Jacobianos,
estacionário, contornos/álgebra temporal e operadores de volumes finitos.
Os valores de Bhagwat são confrontados com sete fixtures calculadas
independentemente em 80 dígitos; as derivadas usam verificações por perturbação
no mesmo ramo físico. `scripts/generate_bhagwat_fixtures.py` requer `mpmath`.

## Benchmark da dissertação

Use os arquivos em `input/benchmarks/`. A especificação completa, os números de
página e as divergências entre fontes estão em
`docs/BENCHMARK_REFERENCE.md` e `docs/benchmark_sources.json`.

| Caso | Pressão do separador (Pa abs.) | Temperatura (K) | jl de referência (m/s) | jg padrão (m/s) |
|---|---:|---:|---:|---:|
| D_Stable | 420000 | 308.45 | 0.594 | 0.690 |
| D_OSC | 206000 | 293.15 | 0.547 | 1.160 |
| D_SS3 | 201000 | 292.95 | 0.280 | 0.480 |
| D_SS1 | 208000 | 293.25 | 0.530 | 0.110 |

A geometria usa Lp=57.4 m, inclinação do oleoduto −2°, X=6.435 m,
Z=9.886 m e D=0.05254 m. O comprimento calculado do riser é
12.5462539932 m. `nodes` e `n_riser` são os números de células físicas de cada
trecho; células fantasmas não entram nessa contagem.

O perfil principal mantém a EOS declarada na dissertação, `rho_g=p/343^2`.
Os pares `_ISO` alteram somente a EOS gasosa para `rho_g=p/(Rg*Tgas)` e são
uma análise de sensibilidade identificada. A vazão de gás é convertida a
fluxo mássico usando 101300 Pa e 293.15 K; ela não é recalculada a partir da
pressão da base ou do separador. As entradas `mul` e `mug` são fluxos de massa
por área em kg/(m² s); `mu_l` e `mu_g` são viscosidades em Pa s.

```sh
./bin/prims_f input/benchmarks/D_Stable.nml
python3 scripts/run_validation.py --case D_OSC --cells 51 --dt 0.5 --time 2400
```

`run_validation.py` cria um deck, um diretório de resultados e um manifesto por
execução. Nomes existentes são protegidos; `--overwrite` é uma ação explícita.
O processo retorna falha se o solver não alcançar a duração solicitada.

## Estudo de malha, passo e EOS

```sh
python3 scripts/run_study.py --workers 4 --time 2400 --label final
python3 scripts/run_study.py --study osc-refinement --workers 4 --time 2400 --label refined
python3 scripts/analyze_validation.py --pattern '[fr]*' \
  --output validation/analysis_final --plots
python3 scripts/compare_validation.py --metrics validation/analysis_final/metrics.json \
  --pattern '[fr]*' --output validation/analysis_final
python3 scripts/analyze_phase_sequence.py --metrics validation/analysis_final/metrics.json \
  --pattern '[fr]*' --output validation/analysis_final
```

O estudo contém 24 tentativas: quatro casos, três malhas (51, 101 e 201
células em cada trecho) a dt máximo=0.5 s, três tetos de passo
(0.5, 0.25 e 0.125 s) na malha 51+51, e quatro sensibilidades isotérmicas.
O ponto de base é compartilhado entre os refinamentos. `--resume` reutiliza
somente tentativas finalizadas com hashes de fontes/executável e duração
compatíveis; uma tentativa que falhou continua classificada como falha.

O estudo adicional `osc-refinement` contém quatro tentativas motivadas pela
oscilação observada em OSC na malha51+51 com dt máximo0.125 s: passos0.0625 e
0.03125 s nessa malha, e malhas101+101/201+201 com dt máximo0.125 s. Essa
verificação distingue a resposta oscilatória obtida com passo menor do
amortecimento obtido com passos maiores; o resultado não é decidido pelo
rótulo do caso.

O pós-processamento de gráficos/espectro requer `numpy`, `scipy` e
`matplotlib`. A janela nominal é 900–2400 s, sujeita à verificação de
acomodação. O analisador principal atribui uma janela diagnóstica a execuções
curtas. Os scripts de comparação e sequência mantêm a janela solicitada; use
`--window-start 0` explicitamente para examinar o trecho inicial de um caso
incompleto. Essas janelas curtas não aprovam ciclos acomodados. O limiar de 100 Pa na detecção de picos é uma regra de análise
declarada; ausência de picos acima dele não prova estabilidade experimental.

## Método e consistência física

As conservativas são as massas de líquido e gás por volume e o momento da
mistura. As duas EOS e a mesma lei de fechamento determinam as primitivas;
o Jacobiano diferencia essa recuperação implicitamente. Densidades,
velocidades e derivadas respeitam as escalas de adimensionalização.

O caminho qualificado usa volumes finitos de primeira ordem e Euler implícito.
`use_muscl=.true.` retorna diagnóstico: a reconstrução legada não integra este
caminho de verificação. Os comprimentos de controle vêm das faces reais da
malha. Todas as células físicas fazem parte do sistema não linear.

O fluxo utiliza uma matriz de Roe integrada por quadratura e sua decomposição
espectral com LAPACK. O sistema em blocos é resolvido com pivotamento por
`DGBSV`. Autovalores complexos, sistemas singulares e raízes ausentes
retornam status; não são substituídos por identidades ou estados recortados.

`dissipation_mode='LAGGED'` congela a viscosidade numérica no estado aceito no
início de cada passo, tanto no residual quanto no Jacobiano. Os fluxos físicos
permanecem implícitos. Essa escolha evita usar uma dissipação descontínua
durante Newton ao atravessar ramos de derivadas do fechamento; sua precisão
temporal precisa ser avaliada por refinamento. A opção `FULLY_IMPLICIT` é
mantida como variante explícita.

No contorno de entrada são impostos os dois fluxos de massa; a pressão vem
do estado interior e a fração de gás é recuperada com continuação do ramo.
Na face de saída é imposta a pressão do separador, preservando fração e
momento interior para recuperar as velocidades no mesmo fechamento.

Na junção, `MATCHED` usa pressão e fluxos de massa comuns e dois traços de
fração de gás. Os fluxos de momento podem diferir, representando a hipótese
de acoplamento com força pontual. Essa construção é testada algebricamente;
ela não é uma solução exata de Riemann para todos os estados. `CONSERVATIVE`
adota outra hipótese discreta, com fluxo único. `BLOCKING` é uma projeção
experimental para modelos específicos: não equivale ao modelo de frente
móvel de Baliño e não é usado nos resultados principais.

## Fechamentos e domínio

| ID | Fechamento | Uso e limitação |
|---|---|---|
| 1 | Estratificado, balanço de tensões | Oleoduto; seleção explícita da raiz |
| 2 | Bendiksen conforme Baliño 2010 | Alternativa com salto em Fr=3.5; pode haver lacuna na inversão conservativa |
| 3 | Bhagwat–Ghajar conforme Azevedo, Baliño e Burr 2015 | Riser dos benchmarks; domínio cocorrente qualificado |
| 4 | Variante reproduzida por Montgomery | Identificador distinto para evitar mistura de fórmulas |

A fonte efetivamente transcrita para ID3 é a seção 6.3, equações 116–128,
p.121 de IMAACA2015. O artigo original de 2014 não foi auditado integralmente.
Veja `docs/BHAGWAT_SOURCE_STATUS.md`. O domínio qualificado permite uma fase
estagnada e exige velocidades não negativas. Uma tentativa de líquido
descendente com gás ascendente exige um fechamento contracorrente validado.

`regime` é um rótulo de saída. Ele não impõe SS1/SS3/OSC nem força ciclos.
Todos os casos recebem a mesma perturbação inicial localizada de 1% no
momento; não há perturbação periódica nem ajuste dos alvos experimentais.
As antigas chaves de histerese não ativam uma troca silenciosa de fechamento.

## Controles e saídas

Use `dt_start_s`, `dt_max_s` e `dt_min_s` para controles em segundos.
`dt_min_s>0` tem precedência sobre o piso adimensional legado `dt_min`.
As entradas legadas positivas `dt_max` e `dt_start` ainda têm precedência
sobre `dt_max_s` e `dt_start_s`; os decks do estudo não as ativam.
O estudo usa piso de falha 1e−7 s para interromper estagnação numérica.
O passo efetivo é adaptado pelas iterações de Newton e rejeições, respeitando
o teto. O piso é verificado ao rejeitar e reduzir um passo. `CFL` e algumas
chaves legadas são aceitas para compatibilidade, mas não controlam o passo
desse integrador implícito; não interpretar o teto como passo sempre fixo.

Cada diretório contém:

- `monitor.csv`: pressão, fração de gás e líquida, fluxos na base/saída,
  inventários totais/por trecho, integrais de entrada/saída e resíduos de massa.
- `fields.csv`: campos espaciais em SI e tempos efetivamente registrados.
- `run_metadata.txt`: namelist efetivo, geometria, escalas e massas impostas.
- `run_status.txt`: tempo atingido/solicitado, status e passos aceitos/rejeitados.
- `last_rejected_state.dat`: estado conservativo da última rejeição, quando há.
- `manifest.json` e `stdout.log`: proveniência e resultado do processo, quando
  a execução é feita pelo runner Python.

O intervalo de saída é um mínimo entre amostras; não cria dados entre passos.
O espectro usa interpolação explicitamente identificada. O balanço integral
emprega os fluxos efetivamente aceitos em cada passo, sem depender da
frequência de gravação do CSV. A análise da sequência física exige examinar
acumulação, compressão, produção, penetração gasosa e blowdown conforme o regime.
O detector de bloqueio forte não constitui condição necessária para OSC/SS3.
Período ou conservação isolados não aprovam um ciclo físico.
