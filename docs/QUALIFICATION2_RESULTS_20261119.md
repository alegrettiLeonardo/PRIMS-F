# PRIMS-F — continuação da qualificação em 19/11/2026

Continuação de `12eaa07`, em `audit/qualification2-20261119`. O operador
da campanha OSC permanece **69813ee**, binário SHA256
`7442291fd6cdc311627e79500dc5a4804813b42c7bf4daecf88b4b8b2322288a`.
A nova globalização está congelada em **a179e51**, binário
`02ce51462ec3c954cec9364b17444b7279be72ec150f87e307f1b5fc8dd0a4c3`.
São coortes separadas: resultados OSC de 69813ee não qualificam o operador
novo. Não foram misturadas as campanhas 05ad799/490936a.

Fonte experimental prioritária: [Baliño, Burr e Nemoto (2010)](https://doi.org/10.1016/j.ijmultiphaseflow.2010.04.003).

GNU Fortran 13.3.0, Fortran 2018, release `-O2`, LAPACK/BLAS, um thread
por execução. Builds, comandos completos, bibliotecas e hashes estão no
índice de evidência. EOS, geometria, propriedades, atrito, alimentação,
contornos e inicialização foram preservados nos refinamentos. `main`
local continua em `4572f95`; alterações anteriores e dados foram preservados.

## Problema → causa → correção → teste → resultado

| Requisito | Causa/evidência e referência | Implementação/verificação | Resultado |
|---|---|---|---|
| SS1/MATCHED | Dissipação de Roe cancela o transporte convectivo na face 53; propostas recusadas têm gás negativo, enquanto as células têm velocidades positivas. HLL (Harten et al. 1983; diagnóstico por derivação conservativa própria) mantém massas intermediárias positivas mas ainda dá fluxo gasoso negativo no mesmo par mapeado. | Replay do operador original bit a bit; comparação Roe/Rusanov/HLL; direção restrita sem trocar fluxo ou fechamento. | Causa localizada; trajetória **BLOCKED**. |
| Globalização SS1 | O gatilho tardio aceitava amortecimentos ~1e-5 até gastar 25 Newton, embora outra direção resolvesse o mesmo passo. | Gatilho na primeira recusa 3012/1017 com margem válida; fallback ao raio original; regressão a partir do mesmo Uold/dt falha antes e conclui depois. | Regressão **PASS**; opção experimental. |
| SS3/raízes | Propostas Newton examinadas têm limite cocorrente positivo robusto em 80 dígitos. Isso não exclui outras direções ou solução global. | L2 e minimax para uma restrição; residual/Jacobiano originais; teste de fluxo comum, massas e balanço. | Replays específicos **PASS**; trajetória **BLOCKED**. |
| SS3/fronteira | No terminal v3, J escalado tem condição 2≈1,00084. Ótimo linear∞≈1,07498e-8 excede tol 1e-8; uma direção com margem admissível também não reduz o residual suficientemente. | Bits exatos e decimais analisados separadamente; sinais sub-ULP não viram certificados de inexistência. | Limitação local demonstrada; sem prova global. |
| Contornos | Nos estados amostrados do nominal, seis offsets e N401: duas características entrantes na entrada e uma na saída, espectros reais e restrições de posto completo. | Auditoria Fortran independente das células interiores e estados impostos. N401:12.501 frames por lado, todos status 0. | **PASS** no domínio observado; inversões/ausência de fase evolutiva **NOT_VERIFIED**. |
| Conservação estável | Resíduos não nulos acumulam erro; tolerância 1e-8 já falhava no controle anterior. | Novo ACTIVE_SET/tol 1e-10 completa 2400 s; monitor, campos e ledger byte a byte iguais à referência qualificada. | Seis balanços **PASS**. |
| OSC/convergência | Amplitude varia com h e dt, apesar da acomodação em cada discretização. | Quatro malhas e três passos independentes, até 401+401 volumes. | Alvo1% de amplitude **FAIL**. |
| Sensibilidade de contorno | Desvios dos sinais da Tabela 3/caso 13 alteram fortemente o resultado numérico. | Seis offsets constantes, uma variável por vez; duas extensões 7200 s; nenhum forçamento temporal inventado. | Quantificada; validação física **NOT_VERIFIED**. |
| Eq. 15 | Forma literal com rho_l diverge da Eq. 54 com rho_g e falha em dois controles estratificados independentes. | Detector preservado; comparação independente com Taitel–Lee–Dukler 1977. | Literal **FAIL** nesses controles; AUTO/reversão **BLOCKED**. |

A numeração específica das equações HLL não foi conferida visualmente;
a fórmula do diagnóstico é apresentada por derivação conservativa própria.
A referência é SIAM Review 25(1), 35–61, DOI10.1137/1025002.

O algoritmo novo resolve subproblemas com uma desigualdade linearizada,
relacionados a Byrd–Hribar–Nocedal1999/§2, sem implementar ou reivindicar
o método completo de barreira. L2 é tentada primeiro; L∞ usa a igualdade
de Hölder, `r_i=delta*sign(q_i)/||q||1`, `q=A^-T*c`. Há até 18 propostas por
direção e18 originais, máximo 54; o limite 25 Newton permanece. Singularidade,
NaN e escala não representável retornam status. Nenhuma margem negativa
é convertida em zero. Cada proposta exige domínio físico, g>=0 e redução
do residual original. `accepted` no ledger ativo significa descida Newton,
não conclusão do passo temporal. O padrão permanece **LINE_SEARCH**.
[Derivação, snapshots, SI, lambdas e testes](../validation/diagnostics/qualification2_20261119/globalization_flux/ACTIVE_SET_MINIMAX.md)

## SS1/SS3: trajetórias e tolerância

| Caso / tol_res | Evidência final registrada | Conservação observada | Qualificação |
|---|---|---|---|
| SS1 /1e-8 | TIMEOUT externo 900 s; CPU897,75 s; ledger até 78,580750 s. Apenas ~9 ms além da parada histórica. Arquivos interrompidos têm comprimentos diferentes. | Cinco de seis normalizações por inventário falham; máximo 4,763e-6. Há falhas também por massa processada. | **FAIL** conservação; regime **BLOCKED**. |
| SS1 /1e-10 | TIMEOUT900 s; CPU899,26 s; último passo completo 2138 em 78,577569192 s; dt≈2,403e-7 s. | Seis por inventário e seis por massa processada passam; máximo 2,781e-8/2,671e-8. | Trajetória parcial conservativa; regime **BLOCKED**. |
| SS3 /1e-8 | status 6000 em 97,131852742 s, 812 passos aceitos; não completa 110 s. | Doze critérios passam; pior inventário 2,649e-9. | Regime **BLOCKED**. |
| SS3 /1e-10 | status 6000 em 97,131852191 s, 810 aceitos/24 rejeitados; CPU47,38 s. | Doze critérios passam; pior inventário 2,920e-11. | Regime **BLOCKED**. |

TIMEOUT é término de infraestrutura, não falha física demonstrada.
O custo de ~900 s de CPU para milissegundos de avanço mostra a insuficiência
prática das alternativas examinadas; não prova inexistência matemática.
O último snapshot de recusa de SS1 precede a interrupção: não é chamado
estado terminal. Apertar somente tol_res melhora os balanços e não resolve
o avanço. Não há dez ciclos acomodados de SS1/SS3.

Na célula 69, theta≈0,901239 rad desde a horizontal. Chexal1992/p. 262 admite
contracorrente vertical; p. 264 restringe horizontal a cocorrente; p. 268
discute raízes múltiplas. A catenária não está nesse domínio vertical.
[Ousaka et al. 2006](https://doi.org/10.1016/j.nucengdes.2005.12.001), Eq. 13/p. 1291, trata limite de flooding em tubos 16/26mm,
30/45/60° desde a horizontal; não fornece o fechamento completo de holdup
para este problema. [Ghiaasiaan et al. 1995](https://doi.org/10.13182/NSE95-A24084) fornece experimentos inclinados,
mas o texto integral não foi obtido para implementar uma lei verificável.
Não foram inventados coeficientes, regra de retorno, frente livre ou
continuação de Chexal inclinado. Uma alteração dessas modifica primitivas,
fluxo/Jacobianos, estacionário e acoplamento e requer nova qualificação.

## OSC: refinamentos independentes do bloqueio SS1/SS3

A referência N51/dt 0,125 já estava executada; cinco novos refinamentos foram
concluídos. Todos completaram 3600 s, zero rejeições e dez intervalos acomodados nos 13
canais predefinidos. Cada linha tem seis balanços globais e60 por ciclo
aprovados. O marcador usa pressão da **face fixa da junção**.

| Volumes por trecho | dt efetivo nos 10 intervalos(s) | Passos aceitos | T(s) | Semiamplitude(Pa) | CPU(s) |
|---:|---:|---:|---:|---:|---:|
| 51 |0,125|28822|23,108618|9909,017|821,134|
| 101|0,125|28822|23,266446|12770,664|1793,007|
| 201|0,125|28822|23,345495|14011,107|3437,280|
| 401|0,125|28822|23,385040|14596,904|6929,323|
| 51|0,0625|57618|23,265643|13763,825|1547,012|
| 51|0,03125|115215|23,344198|15347,389|3022,722|

São volumes uniformes em arco em cada trecho, h=L/N: oleoduto 57,4 m e
riser 12,546254 m. Em N401, Δs=0,143142145/0,031287416 m. Há802 volumes
físicos, além de fantasmas. A perturbação gaussiana mantém amplitude,
centro e largura físicos; sua amostragem varia com a malha. Históricos de
dt completos estão em steps.csv; a tabela não usa só o teto solicitado.

O último par espacial muda amplitude em **4,013%**; o temporal, **10,318%**,
normalizando pela maior magnitude do par. Os períodos mudam≈0,169% e
0,337%. Ordens observadas h101–401:1,0101/1,0937 para período/amplitude;
temporais:0,9992/1,2835. Não há demonstração de resolução do outro eixo,
nem GCI ou convergência conjunta. Acomodação não resolve essa lacuna.
O balanço gás/riser de N51/dt 0,0625 passa com 9,70359e-7, próximo do limite;
N401 passa com máximo 1,30834e-8. Não se adicionou tolerância aritmética
ao limite físico 1e-6.
[Índice, definições, custos e figura](../validation/diagnostics/osc_refinement_20261119/RESULTS_v2.md)

Campos são comparados em posições físicas comuns, separadamente nos dois
trechos, sem extrapolar pela junção. Entre N201/N401, erro RMS do campo
de pressão no riser≈305,93 Pa e da fração líquida≈0,004996. Isolinhas de
alpha 0,4 têm diferença RMS0,07354 m; alpha 0,5,0,26289 m com mudanças de
topologia em 4/256 fases. São cruzamentos geométricos testados com perfis
independentes, não identidades de frentes materiais. O holdup do primeiro
centro muda de posição com N; `1-alpha_riser_trace` é analisado à parte
no ponto fixo. Tempos reais e interpolação são usados: monitor solicitado
a0,03125 s apresentou também lacunas 0,0625 s. 256bins não são 256 observações
independentes por ciclo.

## Contornos, observação experimental e conservação longa

Tabela 3/caso 13, Baliño–Burr–Nemoto2010, p. 653: offsets Ps±5000 Pa,
jg0±0,06 m/s, jl±0,003 m/s. São desvios dos sinais, sem interpretação de
intervalo de confiança ou frequência de forçamento. Eqs. 78–79/p. 649:
`mug=Pstd*jg0/(Rg*Tstd)`, `mul=rho_l0*jl`. Ps varia sem mudar mug.

| Condição | Duração(s) | T(s) | Excursão média por ciclo(bar) | Resultado |
|---|---:|---:|---:|---|
| Nominal |3600|23,10862|0,19818|Acomodado; seis balanços PASS|
| Ps− |3600|23,13373|0,22712|Acomodado; seis balanços PASS|
| Ps+ |7200|23,07715|0,15875|Acomodado; seis balanços PASS;3600 provisório preservado|
| jg− |3600|23,99813|0,34503|Acomodado; seis balanços PASS|
| jg+ |7200|—|—|Sem ciclos úteis; tol 1e-8 falha cinco balanços por inventário|
| jl− |3600|23,18032|0,20977|Acomodado; seis balanços PASS|
| jl+ |3600|23,03697|0,18568|Acomodado; seis balanços PASS|

No jg+7200/tol 1e-8, erro líquido global 1,4299e-3 kg e máximo relativo ao
inventário≈1,0160e-4. Todos passam por massa processada; isso não basta.
Desde4818,93 s há passos aceitos sem nova solução linear, com residual
entre 5,45e-9 e9,9941e-9 na cauda: patamar de aceitação, não prova de
estabilidade física. A repetição 7200 s altera apenas tol_res para 1e-10:
mesmos 57622 passos/zero rejeições, seis balanços passam, máximo 3,3827e-8.
Continua sem período resolvido. O estudo de tolerância é separado da
comparação física com controles numéricos iguais.
[Sensibilidades e fonte](STATIC_BC_AND_REGIME_CRITERIA_20261119.md)

A fonte registra buffer de ar 2,57 m³, aquisição típica 1500 s/10 Hz e
contornos variáveis (pp.650–653). Faltam sinais sincronizados/leis do circuito
para reconstruir seu forçamento. Pb é pressão na base, sem associação
textual inequívoca P4/P5/coordenada do sensor. O traço é um proxy fixo;
a primeira célula é mostrada separadamente. A média/max de alpha_p exp
não tem ponderação documentada; alpha_p sim do artigo é estacionário
(p. 655), uniforme/fixo no modelo (p. 645). Não foi comparado diretamente
ao inventário médio espacial transiente como se fossem a mesma medição.

Há ainda uma diferença de EOS: o standard físico usa rho_std=1,20403095,
enquanto BAROTROPIC usa rho_g=p/343² e dá 0,86103579 kg/m³ em Pstd.
Para o mesmo mug, jg do modelo nessa pressão seria 1,622088 m/s. A massa
prescrita está mantida; essa diferença não foi corrigida por retuning.
ISOTHERMAL existe como opção, mas não integra esta campanha congelada.

## Confronto experimental e limites físicos

Referências prioritárias: SS1/caso 2,85,6 s; SS3/caso 11,**76 s**;
OSC/caso 13,19,3–33,8 s. Os26 s da dissertação não foram misturados com 76 s
nem usados como alvo; a origem da divergência continua não demonstrada.
Não há erro de período SS1/SS3 a calcular a partir destes transientes.

N401 tem excursão de pressão de base≈0,29195 bar na janela final 300 s;
a banda da curva digitalizada Fig.18(a) é0,53631–0,75978 bar, centro 0,64525.
O cálculo fica abaixo da banda, mesmo após refino. A comparação é
exploratória: a amplitude ainda não convergiu. ΔPexp=Max0,69 bar da tabela
não foi confundido com semiamplitude.

O espectro digitalizado tem bins≈0,02662/0,05324 Hz, razão 2; o calculado
tem dominante≈0,04332 Hz. A razão harmônica não identifica sozinha modos
físicos. São janelas 300 s, Hann, potência relativa por série, resolução
≈0,00333 Hz; espessura/calibração da curva e não estacionariedade acrescentam
incerteza. Os544 pontos gráficos não substituem os 1500 s/10 Hz originais
nem os períodos publicados. Sinais abaixo 100 Pa não recebem período
qualificado só porque a potência foi normalizada.

Na mesma janela de dez ciclos de N401, o detector anterior encontrou
10 episódios compatíveis com acumulação/compressão, 10 de produção líquida
e10 de blowdown; nenhum reinício após forte bloqueio. O fluxo gasoso mínimo
observado na base foi 1,31405e-3 kg/s, sempre positivo. A regra de reinício
fortemente bloqueado é apenas um motivo semelhante a SS1, não um requisito
de OSC/SS3. Limiares de análise e todos os sinais usados estão no JSON.

Pressão, holdup, inventários, fluxos assinados, campos e isolinhas estão
disponíveis. Gás positivo na base não invalida automaticamente OSC.
Acumulação/compressão, produção e blowdown são tratados como evidências
numéricas conjuntas; a identidade de frentes materiais e a sequência
experimental correspondente permanecem **NOT_VERIFIED**. Nenhum dos
três regimes foi validado experimentalmente nesta etapa.

## Detector e entrega reproduzível

Eq. 15/p. 3 de sistema_oleoduto_riser conserva `(1-hL/D)^2`, rho_l e Sl;
Eqs. 53–57/pp.8–9 mudam a densidade. A Eq. 1/p. 6 de Taitel–Lee–Dukler 1977
usa rho_g e Si=dAl/dh (Eq. 11/p. 9), sob hipótese horizontal/gás rápido.
No domínio declarado, os controles independentes e limites condicionais
estão registrados; não existe aprovação incondicional de nova correlação
nem regra reversa física. **AUTO_BURR15 continua exclusivamente diagnóstico**.

Builds debug/release a179e51,15 programas Fortran, 76 testes Python no
freeze, 21 verificações CLI e5 de telemetria passaram. A suíte do HEAD
com os analisadores tem 103 testes. Testes de CI executados
localmente são distintos das campanhas longas e não validam experimentos.

```sh
make FC=gfortran LIBS='-llapack -lblas'
make check_all FC=gfortran LIBS='-llapack -lblas'
make test FC=gfortran LIBS='-llapack -lblas'
PYTHONPATH=scripts python3 -m unittest discover -s tests -p 'test_*.py'
python3 tests/test_pipe_regime_cli.py --binary bin/prims_f_check_all --output tmp/NOVO_CLI
python3 tests/test_accepted_telemetry_cli.py --binary bin/prims_f --output tmp/NOVA_TELEMETRIA
python3 scripts/analyze_accepted_conservation.py --run-dir RUN --output NOVO_BALANCO.json
```

O índice da entrega contém comandos de campanha, decks efetivos, SHA,
compilador/flags, dependências, hashes, passos/rejeições, CPU/wall e razões
de término. Raiz externa: `/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119`.
O primeiro controle a179 falhou antes de simular por falta de newline
final no namelist gerado; a tentativa permanece preservada e r2 corrigiu
somente o gerador. PDFs bibliográficos e executáveis não estão no Git.

O GitHub registra PR#2 incorporado externamente em 19/11 às 15:38:32UTC,
merge `e4c6aeb`; comparação via conector mostrou árvore igual a ee 5810e.
O fetch tentou ler main e falhou por conexão; o checkout local foi preservado.
PR incorporado não pode continuar draft. Não houve push, merge ou
publicação remota por esta execução. Os novos commits permanecem locais.


Índice global: [manifestos, fontes e figuras](../validation/diagnostics/qualification2_20261119/delivery/index.json).
Figuras: [pressão N401](../validation/diagnostics/qualification2_20261119/delivery/figures/n401_pressure.png),
[campos espaço-tempo](../validation/diagnostics/qualification2_20261119/delivery/figures/n401_space_time.png),
[campos N201/N401](../validation/diagnostics/qualification2_20261119/delivery/figures/n201_n401_fields.png),
[sensibilidades de contorno](../validation/diagnostics/qualification2_20261119/delivery/figures/boundary_pressure.png),
[espectros](../validation/diagnostics/qualification2_20261119/delivery/figures/boundary_spectra.png).
