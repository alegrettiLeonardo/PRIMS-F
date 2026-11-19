# PRIMS-F — continuação, correções e validação em 18/11/2026

**Estado:** Chexal 1992 ar–água está implementado como fechamento selecionável;
as verificações numéricas passaram e as duas repetições da baseline completaram
2400 s. **A validação experimental permanece aberta:** a amplitude da baseline
ainda varia com o refinamento, a catenária não admite a extensão contracorrente
vertical de Chexal e o traço gasoso negativo MATCHED é um limite separado.
A comparação transiente final foi executada com términos e falhas identificados
por versão, conforme os resultados abaixo. O PR deve permanecer **draft** até a qualificação física.

## Estado preservado e versões

Base de continuação: `origin/audit/chexal-continuation-20261118`, `05ad799`, após
o merge do PR #1 em `4572f95`. O checkout original
`/home/leonardoab/Documentos/PRIMS-F` permaneceu em `main`, preservado. O trabalho
foi realizado em `PRIMS-F-chexal`, branch `audit/chexal-implementation-20261118`;
campanhas usam worktrees separados. Não houve reset destrutivo, concatenação
de tentativas interrompidas nem aplicação retroativa de uma correção a dados antigos.

| Versão / etapa | Resultado realmente obtido | Evidência / limite |
|---|---|---|
| `05ad799`, baseline congelada | Bhagwat ID 3; sem alpha-tail/Chexal. Histórico: 28 configurações, 13 completas, 13 falhas e duas interrupções; duas tentativas adicionais tinham timeout. | [Handoff](CODEX_HANDOFF_STATE.md), [resultados históricos](VALIDATION_RESULTS.md). Interrupção de infraestrutura não é falha física. |
| `47433f8`, alpha-tail | Candidato integrado antes da mudança de física: O0/O3, 55 guards e 398 verificações estacionárias por versão; 101 registros antigos idênticos. Revisão O2: 129 contratos. | [Integração alpha-tail](../validation/diagnostics/alpha_tail_integration_20261118/README.md). Erro relativo máximo das novas raízes: 3,263e-13. Domínio Bhagwat preservado. |
| `db8d47c`, recuperação de evidência | `prior_evidence.json` recuperado do ZIP fornecido com o SHA original `f88afc…169af89`. | [Proveniência](../validation/diagnostics/alpha_tail_integration_20261118/history_recovery.json). A ausência inicial e sua resolução permanecem identificadas. |
| `ee4b7c6`, geometria no caminho Roe | Interpolação de U e inclinação; 315 verificações específicas. Snapshot SS3 antes recusado com status 14 passa no replay corrigido. | [Roe e geometria](../validation/diagnostics/roe_geometry_20261118/README.md). Não prova hiperbolicidade global nem elimina os outros limites. |
| `3096c3b` a `868eeb9`, núcleo Chexal | Correlação SI, ramos, CCFL e referências independentes com 60 dígitos; 139 verificações do núcleo na integração anterior à correção de status. | [Implementação e fontes](CHEXAL_1992_IMPLEMENTATION.md). Não contém coeficientes da variante 1997. |
| `40a7d44`, `868eeb9`, `7321e67`, integração | Residual, primitivas, Jacobianos, estacionário, contornos e junção integrados. Evidência anterior: 3331 verificações / dez programas e 19 testes Python. | [Contrato de integração](CHEXAL_1992_INTEGRATION.md), [resumo anterior](../validation/diagnostics/chexal_integration_20261118/fortran_summary.json). Campanha congelada nessa etapa não é atribuída à correção seguinte. |
| `9fec7ca`, `fdff391`, análise | Seleção explícita de coortes; recusa de mistura de fontes/física; integração trapezoidal compatível com NumPy 1/2. | [Seleção de coortes](../scripts/osc_run_selection.py), [testes de evidência](../tests/test_validation_evidence.py). As tolerâncias são escolhas de análise, não calibração física. |
| `0e84f19`, fechamento da baseline | As duas repetições de `05ad799` chegaram a 2400 s, sem rejeições, preservando todos os hashes Fortran. | [Relatório e métricas](../validation/baseline_repeat_20261118/README.md). Compilador histórico 14.2 versus atual 13.3 explicitado. |
| `b0d33cb`, derivadas próximas de velocidade zero | Stencil unilateral no interior do domínio, guarda bilateral adaptativa e mesmas assinaturas. Suíte dessa versão: **3439 verificações / 12 programas Fortran e 19 testes Python**, todos PASS; build release executado. | [Resumo da correção](../validation/diagnostics/chexal_velocity_tail_20261118/README.md), [referência independente](CHEXAL_VELOCITY_DERIVATIVE_REFERENCE.md). Os testes registram o commit pai `0e84f19` e os hashes incorporados em `b0d33cb`; [CI remota 44389057802](https://github.com/alegrettiLeonardo/PRIMS-F/actions/runs/44389057802): SUCCESS. |
| `490936a`, CCFL/status sobre `b0d33cb` | Guarda de interioridade/resíduo/máximo local em CCFL; propagação de NUMERIC/SINGULAR durante a busca de velocidades. Última suíte local: **3447 verificações / 12 programas Fortran e 19 testes Python**, todos PASS; núcleo 147 checks. | [Evidência CCFL/status](../validation/diagnostics/chexal_ccfl_status_20261118/README.md). A evidência registra o pai e os hashes testados; [CI desta versão, 44390781711](https://github.com/alegrettiLeonardo/PRIMS-F/actions/runs/44390781711): SUCCESS. |

O total atual 3447 inclui os 87 checks de cauda/guarda e 21 de referência
independente; os oito checks adicionais são regressões de CCFL/status.
Não se somam novamente os subtotais da tabela. Comandos efetivamente executados:

```sh
make test FC=gfortran LIBS='-llapack -lblas'
python3 -m unittest discover -s tests -p 'test_validation_evidence.py' -v
```

GNU Fortran 13.3.0, Fortran 2018, verificações de runtime e traps de invalid/zero/
overflow estão registrados nos resumos. O último hash do adaptador de fechamento
testado é `2c7e9b7d81651faa6580be044fcf04e8430fc54670afd8b10dd614d2e955d233`.
Esses resultados verificam operadores e regressões; não são dados experimentais.

## Pedido versus realização

| Requisito | Situação | Evidência ou trabalho restante |
|---|---|---|
| Preservar estado local, partir do handoff e isolar versões | Realizado | Checkout original preservado; worktrees de implementação, baseline e comparação separados. |
| Revisar e integrar alpha-tail em commit próprio antes de Chexal | Realizado | `47433f8`; comparação antes/depois e evidência histórica preservadas. |
| Fechamento Chexal 1992 selecionável, mantendo Bhagwat | Realizado no código; qualificação física parcial | IDs 5/AUTO, 6/HIGH e 7/LOW; Bhagwat continua ID 3. |
| Ramos/CCFL, j=0, ausência de fases, sinais, unidades e domínio | Verificado nos casos exercitados | Referências independentes, testes de escala, dois ramos, ida/volta U–primitivas e status explícitos. Enumeração finita não é prova global de unicidade. |
| Corrigir operadores, Roe variável e examinar MATCHED negativo separadamente | Realizado no escopo documentado | Testes de fluxo/Jacobianos/contornos/junção; limite do traço horizontal negativo permanece aberto. |
| Executar Fortran e CI Python | Realizado | Última evidência local: 3447/12 + 19; CI remota confirmada para `490936a` (3447/12 + 19). Testes locais não equivalem à validação experimental. |
| Repetir os dois refinamentos anteriores sem alpha-tail/Chexal | Realizado | Ambas as repetições em `05ad799` completas em 2400 s, com dez intervalos finais analisados. |
| Comparação nova de fechamentos, controles e refinamentos separados | Executada, com falhas explícitas | Coorte final490936a: seis execuções completas e sete falhas; controles estáveis e OSC Bhagwat analisados. Nenhuma falha foi convertida em sucesso. |
| Validar slugging físico e retirar o PR de draft | Não demonstrado | Exige convergência, comparação experimental coerente e resolução dos limites documentados. |

## Modelo, sinais e parâmetros preservados

O residual implementado é `R = alpha*[C0*(jl+jg)+Vgj] - jg`, inclusive quando
`jl+jg=0`, sem divisão pelo fluxo total. Densidades/velocidades entram no núcleo
em SI; os adaptadores recuperam as escalas antes de avaliar a lei. Reynolds são
assinados: `Re_l=rho_l*jl*D/mu_l` e `Re_g=rho_g*jg*D/mu_g`. Valores absolutos
aparecem somente onde a fonte os prescreve. O ângulo do solver é em radianos
desde a horizontal; Chexal usa graus desde a vertical, com transformação de eixo
e sinais quando necessário. Não houve clipping de fração de vazio, fluxos,
raízes ou espectros para permitir uma execução.

Na comparação de fechamentos, mantêm-se, para cada caso, a EOS, as vazões e a
geometria do mesmo deck. Perfil principal: `rho_g=p/343²`, com EOS líquida linear;
`_ISO` é sensibilidade distinta. Catenária D: Lp=57,4 m, declive=2°, X=6,435 m,
Z=9,886 m, D=0,05254 m; Lr≈12,546254 m e inclinação de topo≈76,473695° desde a
horizontal. Logo ela **não contém uma seção vertical**.

No deck D_OSC: rho_l0=997,05 kg/m³, Cl=1498 m/s, mu_l=0,89e-3 Pa·s,
mu_g=1,85e-5 Pa·s, sigma=0,0728 N/m e g=9,81 m/s². A tensão superficial é
herdada do repositório, não um valor alegado como especificado pela dissertação.
O gás de entrada é convertido de condições padrão 101300 Pa / 293,15 K.
Vazões superficiais padrão não são confundidas com velocidades locais. As
constantes e vazões de cada regime continuam nos [decks de referência](../input/benchmarks/),
sem ajuste para forçar slugging. MATCHED/LAGGED/ROE e a perturbação de 1% ficam
identificados em cada manifesto.

## Critérios de leitura dos resultados

- **Término:** `COMPLETED` exige término normal e tempo físico solicitado atingido.
  `FAILED`, `TIMEOUT`, `RUNNING` e interrupção permanecem distintos.
- **Acomodação numérica:** onze picos consecutivos encerrados entre mínimos
  delimitam dez intervalos. O protocolo de pressão usa excursão mínima de 100 Pa
  e RMS mínimo de 1 Pa; esses pisos não impõem um regime físico. Adota-se 1% para
  a variação da semi-amplitude entre o primeiro e o último pico. Médias e excursões
  dos sinais não nulos de fase/inventário são comparadas entre os primeiros e os
  últimos cinco intervalos, também com referência de 1%. Sinal plano não recebe
  ciclos artificiais; controle estável exige análise própria de estacionariedade.
- **Balanço integral por fase:** `E=M(t)-M(0)-Min(t)+Mout(t)`. A escala é o maior
  valor entre inventário inicial, entrada acumulada e 1e-12 kg. O critério efetivo
  é `max|E| <= 1e-8 kg + 1e-6*escala`. Fechar massa não demonstra convergência
  dinâmica nem concordância experimental.
- **Refinamento:** malha e teto de passo variam separadamente; a janela comum
  900–2400 s e os dez intervalos finais são resultados distintos. Amostras são
  ponderadas no tempo. Mudanças de fonte, compilador, EOS ou fechamento devem
  aparecer na comparação, sem misturar coortes silenciosamente.
- **Regime físico:** pressão, fração líquida, fluxos de cada fase, inventários e
  balanço são considerados conjuntamente. Acumulação/compressão, produção,
  penetração e blowdown são inferidos com limites explícitos; pressão oscilante
  isolada não valida slugging. Bloqueio completo de gás não é requisito universal
  de OSC/SS3.

## Baseline congelada: resultados finais

Os dois casos foram reiniciados do zero em `05ad799`, sem alpha-tail/Chexal,
com EOS BAROTROPIC, ID 3, MATCHED/LAGGED/ROE. Não houve retomada por checkpoint.

| Repetição | Término | Passos aceitos / rejeitados | Período nos dez intervalos finais | Semi-amplitude média | Variação entre primeiro/último pico |
|---|---:|---:|---:|---:|---:|
| 201+201, dt máximo 0,125 s | 2400 s | 19222 / 0 | 23,345470 s | 14011,547 Pa | +0,002630% |
| 51+51, dt máximo 0,03125 s | 2400 s | 76815 / 0 | 23,344185 s | 15347,394 Pa | −0,000031% |

Nos 13 canais analisados — pressão, holdups de base/saída, quatro vazões de fase
e seis inventários — a maior variação relativa das excursões médias entre os
primeiros e os últimos cinco intervalos é **0,002773%** e **0,000640%**,
respectivamente. As duas janelas passam a referência de 1%; isso não estabelece
um ciclo limite assintótico ou validação experimental.

Erros máximos de massa líquido/gás: N201, **9,044925e-10 / 9,227064e-11 kg**;
N51, **2,725528e-9 / 2,370512e-9 kg**. O maior erro relativo entre as duas
repetições e fases é **3,2619e-10**, abaixo do critério de balanço.

Comparados aos casos históricos completos, o refino 101→201 em dt=0,125 s
aumenta a semi-amplitude em **9,725%**; o refino 0,0625→0,03125 s em N51,
em **11,505%**. Os períodos mudam aproximadamente 0,34% nas duas comparações.
**A amplitude ainda não está demonstrada convergente.** As fontes e entradas
físicas coincidem, mas os históricos usam Fortran 14.2 e as repetições, 13.3:
há também mudança de compilador/binário, explicitamente preservada como limite
da reprodução. A comparação não isola experimentalmente esse efeito.

Em 900–2400 s, os proxies de acumulação/compressão, produção e blowdown são
65/62/65 em N201 e 65/64/65 em N51. Os mínimos do fluxo gasoso na base são
0,00136519 e 0,00124482 kg/s: a transferência gasosa continua positiva. Não foi
inventado um bloqueio/reinício para completar uma sequência predeterminada.
Não se mediu posição de frente nem temporização experimental das fases.

Dados, hashes, entradas e critérios completos:
[README da baseline](../validation/baseline_repeat_20261118/README.md),
[summary.json](../validation/baseline_repeat_20261118/summary.json),
[diferenças de refinamento](../validation/baseline_repeat_20261118/last_ten_refinements.csv),
[manifesto do ZIP local](../validation/baseline_repeat_20261118/artifact_manifest.json).
O arquivo `PRIMSF_Baseline_Repeat_20261118.zip` permanece fora do Git, com séries,
campos e figuras próprias; não contém PDFs bibliográficos nem executáveis.

## Comparação final Bhagwat/Chexal

A coorte homogênea final usa **490936a**, com as correções de derivada e CCFL/status.
As 13 execuções terminaram: **seis COMPLETED e sete FAILED, nenhum timeout**.
A primeira coorte `7321e67` permanece preservada em
[relatório anterior](../validation/comparisons/frozen7321e67_20261118/README.md),
com [suplemento de acomodação dos 13 canais](../validation/comparisons/frozen7321e67_accommodation13_20261118/README.md).
Não se atribuem replays ou dados anteriores ao código corrigido.

| Caso base, N51+51/dt máximo 0,125 s | Bhagwat 3 | Chexal 5 |
|---|---|---|
| D_Stable | COMPLETED, 2400,000000 s | COMPLETED, 2400,000000 s |
| D_SS1 | FAILED, 78,571788 s | FAILED, 52,795265 s |
| D_SS3 | FAILED, 97,131852 s | FAILED, 42,182010 s |
| D_OSC | COMPLETED, 2400,000000 s | FAILED, 93,100004 s |

No SS1, o traço gasoso numérico horizontal de MATCHED é negativo enquanto as
células permanecem admissíveis. Isso não demonstra retorno gasoso físico nas
células e não autoriza clipping. No SS3 Bhagwat, a recuperação local não detecta
raiz admissível; NO_ROOT não prova inexistência global de raiz.

No SS3 e OSC Chexal, incluindo os dois refinos, o último estado tem ul positivo
minúsculo, mas todas as cinco propostas Newton seguintes retornam domínio 61
em células da catenária inclinada. A antiga falha de derivada 13 não aparece.
São propostas recusadas, não escoamento reverso físico já aceito ou validado.
Não se atribuem períodos aos transientes truncados.

| OSC Bhagwat | Duração (s) | T dez ciclos (s) | Semi-amplitude (Pa) | Variação da amplitude (%) | Maior mudança de excursão entre 5/5 ciclos, 13 canais (%) |
|---|---:|---:|---:|---:|---:|
| N51, dt=0,1250 | 2400 | 23,103023 | 9700,027 | 2,397741 | 1,693924 |
| N51, dt=0,0625 | 2400 | 23,265639 | 13763,965 | 0,000523 | 0,001185 |
| N101, dt=0,1250 | 2400 | 23,266397 | 12769,700 | 0,021282 | 0,018308 |
| N51, dt=0,1250 | 3600 | 23,108641 | 9909,563 | 0,006414 | 0,011806 |

A extensão até 3600 s passa o critério 1% de pressão e as comparações de médias/
excursões nos 13 canais. Também passam os refinos h e dt; a malha grossa em 2400 s
ainda não passa. Comparadas à malha grossa acomodada em 3600 s, as amplitudes
refinadas diferem **28,862% em h e 38,896% em dt**. A acomodação
temporal está demonstrada nas janelas; a convergência de amplitude permanece aberta.

Os controles estáveis têm pressão média de **501,767033 kPa (Bhagwat)** e
**495,775846 kPa (Chexal)** em 900–2400 s, com pressão e
inventários constantes nas amostras. O maior erro relativo de balanço por fase
entre as execuções completas é **3,0789e-07**, abaixo do critério 1e-6;
isso não significa conservação à precisão de máquina em todos os controles.
Erros absolutos em kg, normalizações e balanços por fase constam no summary.json.

A análise do ciclo médio acompanha pressão, holdups de base/saída, inventários,
vazões por fase e seus atrasos: aumento conjunto de pressão/inventário caracteriza
proxies de acumulação/compressão; produção líquida acompanha depleção do riser;
queda de pressão acompanha depleção gasosa. O fluxo de gás na base permanece
positivo nas amostras. A posição da frente de penetração não foi medida e a
sequência observada não constitui validação experimental de slugging.

A [auditoria independente das entradas](../validation/diagnostics/comparison_input_audit_20261118/README.md)
confere fontes/binário e física iguais entre fechamentos. O refino temporal
mantém a malha, mas muda também a cadência solicitada de gravação de 0,1 para
0,0625 s. O código confirma que essa opção só grava após o passo aceito e não
altera dt/estado; é uma diferença de amostragem explicitada nas métricas.

[Resultados finais e reprodução](../validation/comparisons/frozen490936a_20261118/README.md),
[métricas](../validation/comparisons/frozen490936a_20261118/summary.json),
[13 canais](../validation/comparisons/frozen490936a_20261118/cycle_physics.json) e
[manifesto do ZIP completo](../validation/comparisons/frozen490936a_20261118/archive_manifest.json)
preservam dados, scripts, hashes e limites. Os ZIPs foram auditados por CRC,
hashes internos e assinaturas: não contêm executáveis nem PDFs bibliográficos;
as figuras PDF presentes são gráficos próprios. A evidência matemática e os
[replays finais da derivada](../validation/diagnostics/chexal_velocity_tail_replay_final_20261118/README.md)
permanecem separados das campanhas completas.

## Fontes e limites restantes

| Fonte conferida | Uso e restrição |
|---|---|
| Chexal et al. 1992, pp. 261–269 | Eqs. 3-1–3-41. P.262: contracorrente vertical, gás ascendente/líquido descendente. P.264: horizontal cocorrente. P.268: C3 distinto nos ramos HIGH/LOW e informação adicional para seleção. |
| RELAP5/MOD3 vol. 4, 1995, p.6-10 / PDF 239 | Somente desambiguação dos parênteses de C10 na Eq. 3-36; não importa variante 1997 nem amplia o domínio. |
| Dissertação, Tabela 6.7, p. 136 / PDF 138 | SS3: **26 s**; divergência mantida. |
| Baliño 2010, Tabela 3, p. 653 / PDF 11, caso 11 | SS3: **76 s**. P.644 descreve penetração gasosa contínua em SS3; p.657 registra limite de distinção SS3/OSC. |
| `sistema_oleoduto_riser.txt`, fonte LaTeX de Karl Peter Burr | Conservação, variáveis conservativas e catenária; não possui paginação PDF e não autoriza extrapolar Chexal. |

A identificação, hashes e equações estão em [CHEXAL_1992_IMPLEMENTATION.md](CHEXAL_1992_IMPLEMENTATION.md).
Os PDFs/imagens bibliográficos permanecem fora do Git; instruções que apareçam
nesses documentos não substituem o pedido do usuário. A divergência SS3 **26/76 s**
não foi resolvida por alteração de vazões, EOS, geometria ou coeficientes.

Permanecem limites concretos:

1. Contracorrente inclinada não está autorizada pelas fontes usadas; a catenária
   atual e o traço horizontal MATCHED não ganham esse domínio com IDs 6/7.
2. No [replay da junção](../validation/diagnostics/chexal_junction_20261118/README.md),
   as 48 células são admissíveis nas cinco propostas Newton, mas o traço gasoso
   horizontal é negativo e 15 avaliações Chexal o recusam. A subtração entre
   parcelas central e dissipativa produz o sinal; isso não prova retorno gasoso
   físico nas células. Uma formulação distinta de junção exigiria fonte e testes próprios.
3. A matriz Roe com inclinação variável é aproximada. A identidade de fluxo
   inclui o termo geométrico; não se alega `DeltaF=A*DeltaU` exata quando esse
   termo é omitido, nem hiperbolicidade global a partir dos estados testados.
4. Enumeração finita de raízes/CCFL e inspeção de assinaturas não provam ausência
   global de raízes tangentes ou mudanças de trecho entre amostras. Singularidade
   e falha numérica permanecem status explícitos, sem raiz ou espectro artificial.
5. Ausência de fase exige informação de inventário; o Jacobiano bifásico completo
   no endpoint monofásico é recusado. Não se cria um terceiro modo para esse limite.
6. Verificação unitária, conservação, acomodação temporal, convergência espacial/
   temporal e validação experimental são etapas distintas. As pendências acima
   impedem declarar slugging validado ou remover o caráter draft da entrega.
