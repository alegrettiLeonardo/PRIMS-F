# Seleção de fechamento no oleoduto

Base `ee5810e` do PR #2. Operador congelado para verificação e campanha em
`8a85e37319995b37795d2c798d6f865b863b1834`.

## Diagnóstico e implementação

Em `ee5810e`, o driver só expunha `closure_riser` e o oleoduto usava ID 1.
O índice de regime devolvido por `compute_primitive_variables` era um rótulo
do fechamento escolhido. `PRIMITIVE_UNSUPPORTED_SWITCH` bloqueava a antiga
opção de troca sem lei qualificada. Os parâmetros `stab_*` eram aceitos para
compatibilidade. A seleção efetiva do transiente usava o sinal de `theta`,
que deixava o oleoduto horizontal com operadores diferentes no estacionário,
no interior e na entrada. A saída `results_*.dat` também reconstruía o riser
com Bendiksen implicitamente, inclusive quando o transiente usava outro ID.

| `pipe_regime_mode` em `&PARAMS` | `closure_pipe` | Operador |
|---|---|---|
| `FIXED_STRATIFIED` (padrão) | 1 ou omitido | Equilíbrio estratificado local |
| `FIXED_DRIFT` | ID explícito 2–7 | 2 Bendiksen, 3 Bhagwat, 4 Montgomery, 5 Chexal1992 AUTO, 6 HIGH, 7 LOW |
| `AUTO_BURR15` | 1 ou omitido | Estratificado, com detector em `burr15.csv` |

O ID do oleoduto chega ao estacionário, à recuperação de primitivas, ao
fluxo, à fonte e respectivos Jacobianos, à entrada e aos dois traços da
junção. A seção da malha resolve a ambiguidade de ângulo zero. Os operadores
avulsos mantêm a convenção angular quando a seção não é informada. Na junção
`MATCHED`, o mapeamento usa a mesma pressão e os mesmos fluxos de massa por
fase nas duas leis, com duas frações de vazio admissíveis. O caminho Roe
usa o fechamento do oleoduto após o mapeamento. A saída final usa os mesmos IDs.

A Eq. (15), p. 3, compara a velocidade real do gás, não a superficial,
com `(1-hL/D)*sqrt((rhol-rhog)*g*cos(beta)*A*alpha/(rhol*Sl))`.
A forma quadrada tem `(1-hL/D)**2`. Nesta edição, as Eqs. (53)–(54)
implicam trocar `rhol` por `rhog` no denominador. Geometria, páginas,
proveniência, referências numéricas e a distinção entre o perímetro molhado
`Sl` e a largura interfacial estão na [nota de fontes](BURR15_SOURCES.md).

O detector avalia o candidato estratificado obtido do estado conservado
aceito. Cada diagnóstico computa `hL/D`, `ug**2`,
`ug_crit**2` e a margem em m²/s². `unstable=1` indica margem não negativa
com `status=0`. O solver usa `theta` em radianos da horizontal, com
`theta=-beta` no oleoduto descendente. O detector aceita geometrias
circulares horizontais e descendentes com `cos(theta)>0`, velocidades
não negativas, duas fases e `rhol>rhog>0`. Essa guarda delimita a avaliação
da fórmula. A faixa de validação física da Eq. (15) permanece aberta.

| Status registrado | Significado |
|---:|---|
| 0 | Avaliação válida |
| 1 | Estado ou propriedades inválidos |
| 2 | Fase ausente |
| 3 | Velocidade negativa em pelo menos uma fase |
| 4 | Ângulo fora da avaliação horizontal/descendente |
| 5 | Candidato estratificado indisponível |
| 6 | Área e diâmetro incompatíveis com seção circular ou inválidos |
| 7 | Avaliação sem representação numérica adequada |

`burr15.csv` acompanha as amostras de campo das células do oleoduto.
As velocidades assinadas, densidades e frações originais estão em
`fields.csv`. Leituras do detector com status diferente de zero não
classificam a estabilidade. A continuidade de pressão e de fluxos é
verificada no traço da junção. Uma proposta de Newton ou de Roe recusada,
incluindo gás negativo no traço, não é uma observação de contrafluxo físico.

## Limites da seleção

A fonte fornece um teste de perda de estabilidade do estratificado e
diverge na densidade usada entre suas formas dimensional e adimensional.
Falta um critério de retorno aplicável aos estados não estratificados do
transiente. Portanto a seleção automática para produção está **BLOCKED**.
Avaliar o candidato estratificado com lei fixa evita a circularidade no
detector. Usar o complemento da desigualdade como retorno requer uma lei
adicional e pode produzir comutação espúria se as duas leis recuperarem
velocidades distintas a partir do mesmo momento. Nenhuma banda de histerese
é adicionada. As correlações de drift conservam seus domínios, inclusive
a restrição de Chexal1992 em contracorrente inclinada.

`FI_FANNING`, atrito de parede, EOS, condições físicas, controladores e
tolerâncias do solver permanecem os da base para cada deck. O novo detector
usa a mesma trajetória que o modo fixo estratificado, com diagnósticos
avaliados fora de Newton e das equações de conservação.

## Reprodução e critérios da campanha

```sh
make test FC=gfortran LIBS='-llapack -lblas'
make -j4 FC=gfortran LIBS='-llapack -lblas'
make check_all FC=gfortran LIBS='-llapack -lblas'
python3 -m unittest discover -s tests -p test_validation_evidence.py -v
PYTHONPATH=scripts python3 -m unittest discover -s tests -p 'test_burr15*.py' -v
```

Usar os decks `input/benchmarks/D_*.nml` BAROTROPIC, Bhagwat no riser,
`ROE`, `MATCHED`, `LAGGED` e a perturbação publicada no deck, com
saídas diferentes para cada execução. Congelar o SHA, parâmetros, binário,
fontes e bibliotecas. A série inicial usa 51+51 células, 2400 s e
`dt_max_s=0.125`. Refino espacial altera apenas as contagens para 101+101.
Refino temporal altera apenas o teto para 0,0625 s. Conferir o passo aceito
nas amostras de `monitor.csv` e os totais de `run_status.txt`. Encerramento
normal com o tempo físico alvo atingido é necessário para qualificar ciclos.

Comparar a identidade dos CSVs físicos entre `FIXED_STRATIFIED` e o modo
diagnóstico. Para OSC, usar onze picos consecutivos cercados por mínimos
para dez ciclos completos. Exigir excursão de pelo menos 100 Pa e RMS de
pelo menos 1 Pa. Acomodação exige variação de até 1% da amplitude de ponta
a ponta e de médias e excursões não nulas nos 13 canais entre as primeiras
e as últimas cinco janelas. A janela 900–2400 s serve para comparação no
mesmo intervalo. Os critérios de acomodação não demonstram convergência.

O balanço por fase é `M(t)-M(0)-Min(t)+Mout(t)`, usando fluxos integrados
pelo solver. Reportar valores absolutos, por inventário inicial e por
massa líquida processada, além das diferenças por ciclo. Os inventários
do oleoduto e riser são separados. Sem cumulativos aceitos na base,
balanços independentes por trecho têm limitação de amostragem.

As observações espaciais incluem cruzamentos de níveis `alpha=0.4,0.5,0.6`
interpolados entre centros da mesma seção. Registrar posição e orientação
sem transformar um salto da junção em frente. São candidatos geométricos
dependentes da discretização e do nível. A interpretação de acumulação,
produção, penetração e blowdown usa a série de fluxos e inventários,
com o protocolo de `scripts/analyze_phase_sequence.py`. Pressão
oscilante e classe local do detector são insuficientes para confirmar
severe slugging. A discrepância SS3 de 26 e 76 s continua explícita no
[relatório anterior](PRIMSF_CONTINUATION_VALIDATION_20261118.md).
