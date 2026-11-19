# P4/P5: evidência primária e comparação condicional

Conferência de 2026-11-19. Nenhum deck, fechamento ou parâmetro do solver foi
alterado. Os PDFs e as imagens originais continuam fora do Git. Esta nota
complementa [a leitura de Burr15](BURR15_SOURCES.md), sem escolher outra fórmula.

## P4 — controle experimental independente do detector literal

[Carraretto et al. (2020), DOI 10.3390/fluids5030101](https://www.mdpi.com/2311-5521/5/3/101)
mediram ar–água em tubo horizontal de 60 mm. As Eqs.1–2/p.5 relacionam altura
e fração de gás; pp.5–6 descrevem a sonda, repetição e precisão de altura de
0,5 mm. A Tabela4/p.6 dá incertezas relativas até 7% para JG e 3% para JL.
Na p.7, as ondas aparecem somente perto de JG=2 m/s; os dois pontos abaixo
ficam abaixo de 1 m/s mesmo com o envelope ampliado. São observações
estratificadas suaves, não pontos tabulados de transição.

Leitura visual deliberadamente ampla dos terceiros marcadores da Fig.9/p.10:

| Painel | JL nominal (m/s) | Intervalo lido de xV | Intervalo lido de alpha |
|---|---:|---:|---:|
| 9(a) | 0,03 | 0,960–0,967 | 0,76–0,80 |
| 9(b) | 0,04 | 0,948–0,954 | 0,71–0,76 |

PDF local: SHA256
`960446bbdecef9acd5687978d8918390dde00ac58f045c2a00968144abbb7887`.
Os retângulos são leituras gráficas conservadoras, não dados brutos ou
intervalos estatísticos de confiança. O diâmetro é o nominal publicado.

O [calculador independente](../scripts/assess_burr15_primary_bounds.py) usa
`xV=JG/(JG+JL)`, amplia os limites de JG por ambos os máximos de incerteza
de vazão e os de alpha por `4*0,0005/(pi*D)=0,01061033`, limite geométrico
de `|dalpha/dh|`. Inverte a área circular por altura e usa `Sl=D*phi`.
Para limitar a Eq.15 por cima, usa `(rhol-rhog)/rhol<=1`, sem assumir
densidade não documentada, e `g=9,81 m/s²`. À medida que alpha cresce,
`1-h/D`, alpha e `1/Sl` crescem; logo Ucrit cresce. Portanto basta calcular
Ucrit no maior alpha e ug mínimo no menor JG e maior alpha.

| Painel | alpha ampliado | JG ampliado (m/s) | ug mínimo (m/s) | Ucrit15 máximo (m/s) | Margem mínima ug−Ucrit (m/s) |
|---|---:|---:|---:|---:|---:|
| 9(a) | 0,749390–0,810610 | 0,649512–0,968846 | 0,801263 | 0,454787 | **+0,346475** |
| 9(b) | 0,699390–0,770610 | 0,657839–0,914264 | 0,853660 | 0,406785 | **+0,446875** |

Assim, a Eq.15 literal com perímetro de parede aponta perda de estratificação
em todo o envelope, enquanto a observação é estratificada. Isso é uma
discordância independente robusta nesse escopo horizontal/cocorrente. Não
qualifica a seleção automática, a transição reversa, uma nova correlação ou a
catenária. A [evidência reproduzível](../validation/diagnostics/primary_evidence_p4p5_20261119/burr15_primary_bounds.json)
registra os intervalos e as hipóteses; quatro testes verificam a geometria e
os limites, sem chamar Fortran.

As Eqs.53–57 de Burr e a Eq.1/p.6 do relatório Taitel–Lee–Dukler1977 usam
outra densidade; o relatório também usa `dAl/dh=Si`, em lugar de Sl. Essa
diferença permanece explícita. O resultado acima não autoriza substituir
qualquer uma dessas expressões no solver.

[AL-Dogail e Gajbhiye (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8481959/)
documentam rampas de gás nos dois sentidos, observação de regime e pressão
em tubo horizontal de 25,4 mm. O artigo foi lido, mas a planilha suplementar
retornou uma página de challenge, não XLSX; a descrição das colunas não
fornece altura/holdup. **BLOCKED dataset para qualificar Burr15/reversão:**
faltam estados simultâneos de holdup, velocidade local, regime e direção do
caminho. Não se afirma ausência universal de dados ou de leis reversas.

### Cobertura dos dois controles adicionais

1. **Perda observada de estratificação: BLOCKED dataset.** Carraretto2020,
   p.7, relata início de ondas 2D perto de JG≈2 m/s, ainda estratificado.
   A Fig.1/p.3 inclui regiões de um mapa de referência; esses nomes não
   convertem os pontos medidos em eventos observados de transição. Nas
   passagens conferidas não há par de estados locais pré/pós-transição
   com holdup, JG e perda de estratificação identificada. Não usar o início
   de ondulação como evento não estratificado ou como validação positiva
   da Eq.15. Os dois controles da Fig.9 verificam somente discordância
   onde o escoamento permaneceu estratificado.
2. **Não estratificado documentado, aplicação geométrica recusada.**
   Baliño2010, Tabela3/caso14, p.653, identifica SF; Fig.14/p.655 apresenta
   sua história experimental. Condições: Ps=2,08±0,07 bar absoluto,
   T=18,7±0,1°C, jg0 padrão=1,86±0,04 m/s, jl0=0,82±0,01 m/s e
   alpha_p exp=0,74. Esse registro de regime do sistema não estabelece
   uma interface plana em cada seção ou identifica o estado local de
   transição do oleoduto. **Não** converter 0,74 em hL/Sl nem substituir
   jg0 padrão pela velocidade local ug. Sem esses estados, não há teste
   numérico experimental legítimo de Burr15 para esse caso.

O segundo controle é uma checagem de não aplicabilidade e de configuração,
não uma avaliação da desigualdade. Em [main.f90](../scr/main.f90),
`AUTO_BURR15` exige `closure_pipe=1` (ou omissão); `FIXED_DRIFT` exige
explicitamente ID 2–7. O caso `auto_drift_id` em
[test_pipe_regime_cli.py](../tests/test_pipe_regime_cli.py), AUTO com ID5,
espera rejeição de configuração. A [evidência já executada](../validation/diagnostics/burr15_config_cli_20261119/cli_summary.json)
registra `returncode=1`, `passed=true`; foi conferida nesta revisão, sem
nova campanha. Além disso, o writer só chama o detector para célula do
oleoduto com estado válido e `cid==1`; caso contrário registra
`BURR15_NOT_EVALUATED`, sem fabricar geometria. O núcleo da fórmula recebe
apenas variáveis numéricas e não reconhece experimentalmente o regime:
a premissa de candidato estratificado é responsabilidade do chamador.
Um alpha interior, por si só, nunca comprova essa premissa. A guarda
confirma separação de modos, não valida fisicamente um fechamento drift.


## P5 — Baliño2010: identidade dos casos e condições

[Baliño, Burr e Nemoto (2010)](https://doi.org/10.1016/j.ijmultiphaseflow.2010.04.003),
PDF local SHA256
`ce04a2d1b79d3a10d442fbaeb3863e4639c8508bd34ba9e1a414ae03e5479e2b`.
As páginas impressas 650–660 correspondem a PDF8–18. Tabela3/p.653:

| Caso experimental | Regime | Ps (bar absoluto) | T (°C) | jg0 padrão (m/s) | jl0 (m/s) | ΔPexp (bar) | alpha_p exp | Texp (s) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 2 | SS1 | 2,08 ±0,03 | 19,7 ±0,1 | 0,11 ±0,01 | 0,528 ±0,008 | 0,51 | Max 0,55 | 85,6 |
| 11 | SS3 | 2,01 ±0,07 | 19,8 ±0,1 | 0,48 ±0,06 | 0,28 ±0,03 | 0,86 | Max 0,76 | **76** |
| 13 | OSC | 2,06 ±0,05 | 20,0 ±0,1 | 1,16 ±0,06 | 0,547 ±0,003 | **Max 0,69** | 0,75 | **19,3–33,8** |

Os sinais ± são os desvios reportados das histórias, sem reinterpretá-los
como incerteza instrumental ou intervalo de confiança. `Max` preserva o
qualificador da tabela. As Eqs.78–79/p.649 definem as vazões superficiais
padrão com P0=1,013 bar absoluto e T0=20°C; o deck usa os mesmos valores.
O período SS3 de 26 s da dissertação continua incompatível com os 76 s
desse caso primário; não os combinar nem usar 26 como alvo. A coluna de
erro de período da própria Tabela3/caso11 também é inconsistente: 72,2/76−1
é −5%, embora a tabela imprima −1%.

Nas pp.650–653: oleoduto 57,4 m descendente a 2°, riser catenário cerca de
9,9 m, ar de buffer de 2,57 m³, medição de vazões por turbinas de gás e
medidores magnéticos de líquido, nove transdutores de pressão, densitômetro
gama na base e 14 sondas de condutância. Aquisição típica: 1500 s a 10 Hz.
O separador controlado e a interação pressão/vazão produzem contornos
variáveis. Períodos vêm do espectro de pressão de base; OSC pode apresentar
dois picos e modulação. ΔPexp vem da história de base, com máximo quando
não há sinal claramente definido. Não é pressão média, nem meia excursão.

Comparação dos decks originais, mantidos intactos:

| Item | Fonte primária | Deck atual | Consequência |
|---|---|---|---|
| SS1 T/jl | 19,7°C / 0,528 m/s | 20,1°C / 0,53 m/s | +0,4°C; jl +0,38% |
| SS3 e OSC, médias Ps/T/jg/jl | Tabela3, casos11/13 | Coincidem | Contornos constantes não reproduzem variação experimental |
| D usado na simulação do artigo, p.654 | 0,0525018 m | 0,05254 m | +0,07276% em diâmetro |
| rho_l e viscosidades da simulação do artigo, p.654 | 1000; 0,001; 1,8e−5 em SI | rho_l0=997,05; 0,00089; 1,85e−5 | São escolhas de modelo distintas, não medidas reatribuídas ao experimento |
| Modelo | NPW, líquido incompressível, gás isotérmico, fechamento do artigo | Líquido compressível, BAROTROPIC, Bhagwat | Benchmark relacionado, sem reprodução integral do modelo de 2010 |

SHA256 dos decks conferidos: SS1
`6c7b0f9bdcc0b2d2acca212f50b760dc2613fc8f3b6474ceac7731e240d58f7d`;
SS3 `f5fe9ece006415d71541f1196b243fb047d57abd5739319cf1a81247261c86ee`;
OSC `48b02dffb9319c3786fc5d4c102c1893786e9d3ab8ec668025f1893e26272408`.

## Séries publicadas, digitalização e limite da comparação

Fig.10/p.654 é o experimento SS1/caso2; Fig.12/p.655 é SS3/caso11;
Fig.18/p.657 é pressão/FFT do caso13. A Fig.13 é OSC/caso17, não caso13.
Figs.19–20 são simulações SS1/caso2; Fig.21 é experimento caso4; Fig.22 é
um mapa de estabilidade. Não confundir essas séries ou seus identificadores.

O [script de digitalização limitada](../scripts/digitize_balino2010_case13.py)
extrai somente a curva azul de pressão de base da Fig.18(a), conferida
visualmente, do raster nativo xref263 (1256×434 px). Cobre 1200–1500 s:
544 colunas, nenhuma ausente, cerca de 0,5525 s/px e 558,7 Pa/px. Cada valor
central é o meio do intervalo azul na coluna; a faixa inclui flutuação não
resolvida/espessura, meia célula de pixel e tolerância declarada de um pixel
na calibração. Isso não recupera os dados originais de 10 Hz e não é uma
incerteza metrológica completa. Não sincronizar tempos absolutos nem tratar
esses centros como uma série experimental bruta.

CSV, metadados e dois gráficos próprios ficam somente em
`/home/leonardoab/Documentos/PRIMS-F-numerics-results-20261119/primary_evidence/balino2010_case13/`.
A origem, hash do PDF/raster/script e calibração são registrados no JSON.
O gráfico de período compara a faixa publicada com **23,10864 s**, da
extensão histórica OSC/51+51/dtmax0,125/freeze8a85e37. A inclusão nessa faixa
é **compatibilidade condicional**, não validação. A excursão comum dos
ciclos numéricos não deve ser sobreposta como equivalente ao `Max 0,69 bar`
da tabela; o trecho digitalizado de 300 s também não substitui esse máximo.

O relatório experimental original é Wordsworth et al. (1998), CALtec
CR6820, Vol.I–III, citado na p.660. Sua série bruta não foi obtida.
**BLOCKED validação quantitativa completa:** faltam a aquisição original,
incertezas instrumentais detalhadas, contornos experimentais sincronizados
e compatibilização explícita das diferenças de modelo e configuração.
