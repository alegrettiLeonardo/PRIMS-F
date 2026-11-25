# Chexal–Lellouche: domínio físico e condição de transmissão na junção

Data: 21/11/2026 (terceira rodada). Branch `implementation/chexal-riser-20261121`, base
`5028338` (confirmada); `main` em `f3717dd`. Sem push, merge ou PR. Complementa
[CHEXAL_LELLOUCHE_RISER_QUALIFICATION.md](CHEXAL_LELLOUCHE_RISER_QUALIFICATION.md), Partes I e II.
Evidências em [validation/diagnostics/chexal_riser_20261121/evidence](../validation/diagnostics/chexal_riser_20261121/evidence).

## Resumo

- **Domínio Chexal ar–água (P0.1).**
  - O núcleo tem agora um **contrato executável**: estados fora da base ar–água da fonte
    (ρg/ρl > 8,1e-3, isto é, ~6,8 bar a 20 °C, ou D fora de 0,010–0,305 m) são recusados com
    `CHEXAL_UNQUALIFIED` (68).
  - A divergência EPRI a 1000 psia foi decomposta até a primeira grandeza divergente: a
    **razão de densidades em K0 e r**, seguida do **grupo de propriedades de Vgj**. As
    linhas 2, 6 e 7 implicam os **mesmos** valores (ρg/ρf = 0,0468 ± 0,0001; grupo +2,5 %),
    com todas as formas de equação verificadas. É um conjunto de propriedades não divulgado,
    que não corresponde a um estado saturado coerente. **FAIL mantido** para uso generalizado
    em pressão; o subconjunto qualificado está na §1.
- **T1 × T2 (P0.2): recomendação T1** (G_l, G_g e F₃ contínuos), por quatro razões:
  1. as leis (72)–(75) de Burr e sua discretização nodal (128)–(129) não contêm força singular;
  2. o balanço de volume de controle em duto de área constante, com dobra, não tem força
     tangencial concentrada;
  3. Borda–Carnot é consequência da conservação do momento, não do salto T2;
  4. no estudo de Riemann local com Chexal, a T1 tem solução única em todos os 8 estados
     aceitos, inclusive no **bloqueio do SS1**, onde o gás continua passando
     (G_g* = 2,94e-3 kg/(m²·s)). A T2 não tem solução com α_R* > 0 ali: em 75 partidas,
     todas são atraídas a α_R* → 0 com G_g* ≤ 0.

  O operador Roe/MATCHED atual reproduz a T2 com diferença ≤ 0,6 %. **O operador global não
  foi alterado**: implementar a T1 exige qualificação própria (§10).
- **Banco experimental (P1.1).** Nove fontes rastreadas. Só uma foi acessada: um estudo de
  CCFL, fora do domínio. Nenhum ponto de α utilizável: **NOT_VERIFIED**.
- **Contracorrente inclinada (P1.2): BLOCKED** para o mecanismo B (líquido descendo a 65–76°).
  O mecanismo A (gás revertendo na junção horizontal) **não é exigido sob T1** no estado de
  bloqueio do SS1.
- **Campanhas (P2):** somente regressão do contrato de domínio (§5). Os gates físicos (T1 não
  implementada, P1.2 BLOCKED) impedem campanhas longas novas.

## 1. Domínio Chexal/ar–água e pendências EPRI (P0.1)

### 1.1 Matriz ramo × domínio × evidência × status

| Ramo (Chexal 1992) | Equações | Ativação | Domínio admitido no PRIMS-F | Evidência numérica | Evidência experimental | Status |
|---|---|---|---|---|---|---|
| L ar–água ascendente | 3-7 | Reg ≥ 0 | ρg/ρl ≤ 8,1e-3; D 0,010–0,305 m | núcleo = transcrição independente em 945 estados (≤ 6,9e-14); ângulos (468 verificações) | base ar–água da fonte (Tab. 5-1; só estatísticas dos autores) | PASS numérico; experimental NOT_VERIFIED |
| L ar–água descendente | 3-8 | Reg < 0 | idem | idem | idem (3 conjuntos verticais) | idem |
| L horizontal | 3-19, 3-18 | componente horizontal (Fr < 1) | fluxos ≥ 0 no eixo | idem | 5 conjuntos horizontais | idem |
| L vapor–água | 3-5, 3-6 | `fluid=STEAM_WATER` (só verificação) | não usado pelo PRIMS-F | NP-3989 linhas 1, 4 e 8 | Tab. A-1 | PASS (14,7 psia) / FAIL (1000 psia, §1.2) |
| K0, r, B1, Re (3-11 a 3-17) | — | sempre | — | NP-3989 1, 4, 8 + decomposição | — | PASS |
| C2, C5 < 1 (3-28, 3-29) | — | ρf/ρg > 18 e 150 ρg/ρf < 1 | **ramo de todos os estados do PRIMS-F** (ρg/ρl ≤ 6e-3) | NP-3989 14,7 psia reproduzido | — | PASS |
| C2 = 1 com C5 ≥ 1 | 3-28 | 150 ρg/ρf ≥ 1 | ar–água só entre ~5,7 e 6,8 bar (não atingido pelo PRIMS-F) | NP-3989 1000 psia: forma verificada **condicionada** a propriedades não divulgadas | — | FAIL (generalizado) |
| C2 = 0,4757 [ln(ρf/ρg)]^0,7 | 3-27 | ρf/ρg ≤ 18 | **recusado** (fora da base ar–água) | — | — | fora do contrato |
| C3 ascendente | 3-34 | j_l, j_g ≥ 0 | admitido | NP-3989 1, 4 | — | PASS |
| C3 descendente / CCFL HIGH | 3-35 a 3-40 | ambos < 0, ou contracorrente vertical | admitido (descendente cocorrente; contracorrente só a 90° com ramo) | NP-3989 linha 5 só na forma de 1986 (revisão demonstrada) | CCFL vertical: Bharathan, 25 e 51 mm | PASS numérico |
| C3 LOW | 3-41 | contracorrente vertical, ramo LOW | só a 90° | testes do núcleo (mpmath) | — | PASS numérico |
| C4 | 3-30 a 3-33 | D < D2 | admitido | NP-3989 1 a 4 | — | PASS |
| C9 | 3-25, 3-26 | sinal de Reg | admitido | NP-3989 | — | PASS |
| Fr (inclinação) | 3-2, 3-3, 3-21, 3-22 | 0 < θ_C < 90 | cocorrente em qualquer ângulo | ângulos (reversão de eixo exata, continuidade em θ = 0) | Beggs 1972 (D 25,4–38,1 mm; fora do D do PRIMS-F) | numérico PASS; experimental NOT_VERIFIED |
| Contracorrente inclinada | — | j_l·j_g < 0 com θ_C ≠ 0 | **recusada** (61) | testes de domínio | nenhuma na fonte (p. 262) | fora do contrato |
| Horizontal com fluxo negativo | — | θ = 0 e j < 0 | **recusada** (61) | idem | nenhuma (pp. 264–265) | fora do contrato |
| Derivadas | Richardson de 4ª ordem com guarda de quina | — | dentro de cada pedaço | sonda de fronteira (≤ 9e-10 em j_l → 0⁺; 2,2e-6 em α → 0); quina exata recusada (53) | — | PASS |

**Subconjunto efetivamente qualificado (numericamente):**
- AIR_WATER, cocorrente em qualquer ângulo, com ρg/ρl ≤ 6e-3 (ramo C5 < 1) e D = 52,54 mm;
- contracorrente só vertical, com ramo explícito (não ocorre na catenária).

**Qualificação experimental: nenhuma demonstrada** (§3).

**Contrato executável** (`scr/module_chexal_1992.f90`, `tests/test_chexal_domain_contract.f90`, 19 verificações):
- No modo AIR_WATER, `source_coordinates` e `chexal_recover_velocities` recusam com
  `CHEXAL_UNQUALIFIED` (68) se ρg/ρl > 8,1e-3 ou D ∉ [0,010; 0,305] m.
- O limite corresponde à base ar–água da fonte (Tabela I-I, p. 257: 1–6,8 bar; Tabela 5-1:
  D 11–305 mm). 6,8 bar de ar a 293 K dá ρg/ρl = 8,1e-3.
- Somam-se a isso as recusas de orientação (61) e de ramo (62).
- O teste aceita os extremos do PRIMS-F (ρg = 1,2 e 6 kg/m³, θ = 0° e 76,47°), as bordas
  da base (ρg/ρl = 8,0e-3; D = 0,010 e 0,305 m) e eixo descendente.
- Recusa ρg/ρl de 8,2e-3 e 4e-2, D de 9,9 mm e 310 mm (em coeficientes, raiz, resíduo e
  recuperação), contracorrente inclinada, gás revertendo na horizontal e gás descendo com
  líquido subindo.
- O modo vapor–água de verificação fica fora do contrato ar–água.

### 1.2 EPRI NP-3989 Tabela A-1: decomposição até a primeira divergência

[np3989_decompose.py](../validation/diagnostics/chexal_riser_20261121/np3989_decompose.py),
[tabela](../validation/diagnostics/chexal_riser_20261121/evidence/NP3989_DECOMPOSITION.md).
Os intermediários são avaliados no α **impresso**, com as equações de 1986 (a versão da
tabela; a forma ascendente é idêntica em 1992).

| Linha | P (psia), D (ft), sentido | Erro em C0 | Erro em Vgj | ρg/ρf implícito pelo C0 | Grupo de propriedades implícito pelo Vgj (calculado 0,12245) | Primeira divergência |
|---|---|---:|---:|---:|---:|---|
| 2 | 1000; 0,05; ↑ | −8,0e-4 | −2,5e-2 | 0,04678 | 0,12555 | ρg/ρf em K0 e r |
| 3 | 1000; 1,0; ↑ | −8,6e-4 | +1,4e-2 | 0,04667 | 0,12074 | idem; a linha viola a própria identidade |
| 6 | 1000; 1,0; ↓ | −2,5e-4 | −2,5e-2 | 0,04689 | 0,12555 | ρg/ρf em K0 e r |
| 7 | 1000; 0,05; ↓ | −1,9e-3 | −2,2e-2 | 0,04671 | 0,12526 | ρg/ρf em K0 e r |
| 1, 4, 5, 8 | 14,7 | ≤ 5e-5 | ≤ 3,5e-3 | = tabelas | ≤ 0,4 % | nenhuma (linha 5 na forma de 1986) |

Diagnóstico:
- As quatro linhas a 1000 psia implicam a **mesma** razão de densidades (0,0468 ± 0,0001,
  contra 0,0485 das tabelas de vapor).
- As linhas 2, 6 e 7 (dois diâmetros, dois sentidos, Re de 1,8e5 a 3,7e6) implicam o
  **mesmo** grupo [σ(ρf−ρg)/ρf²]^¼ = 0,1255 ± 0,0002.
- As formas de C2, C3 (inclusive a descendente de 1986), C4, C9 e L estão verificadas.

A divergência é, portanto, de **propriedades** usadas pelo programa de 1986 e não impressas.
Elas não correspondem a um estado saturado coerente: a razão indica ~970 psia e a σ ~273 °C.
Nenhum coeficiente foi ajustado. **Status: FAIL** para o uso generalizado em pressão (ramo
C5 ≥ 1). Esse ramo fica fora do subconjunto qualificado e só é atingido por ar–água entre
5,7 e 6,8 bar, o que o PRIMS-F não faz.

## 2. Derivação física T1 × T2 (P0.2)

### 2.1 Modelo de referência

Burr, *Modelagem unidimensional de sistema oleoduto-riser* (PDF local), eqs. (72)–(75):
∂U/∂t + ∂F/∂s = S, com

- U = [(1−α)ρl, αρg, (1−α)ρl u_l + αρg u_g]ᵀ;
- F₃ = (1−α)ρl u_l² + αρg u_g² + π_P P;
- S₃ = −π_G ρm sen θ(s) − ½ π_F ρm f_m j|j|, **sem termo singular**.

A discretização (128)–(129) integra cada intervalo [s_j, s_{j+1}] com **um único** F_{k,j} por
nó. Somados os intervalos vizinhos, F₃ se cancela em todos os nós, inclusive no de junção.

### 2.2 Balanço integral na junção

Volume de controle [s_J − ε, s_J + ε], área constante A:

d/dt ∫ U₃ ds + F₃(s_J + ε) − F₃(s_J − ε) = ∫ S₃ ds + J_M / A.

- **Termos distribuídos:** ∫ S₃ ds = O(ε) → 0, porque gravidade e atrito são limitados.
  O primeiro termo também é O(ε).
- **Dobra:** θ passa de −2° no oleoduto para 0° na base da catenária. Numa curva brusca sem
  atrito, a reação da parede é normal à bissetriz, então F_R t_R − F_L t_L = W com
  W·(t_L + t_R) = 0, e daí |F₃,R| = |F₃,L|. A dobra **não** produz força tangencial.
- **Troca de fechamento** (estratificado → Chexal): não é uma força física. As forças
  interfaciais são internas e se cancelam no momento da mistura.

Portanto, no limite ε → 0, **J_M = 0**: G_l, G_g e F₃ são contínuos. Isso é **T1**, a condição
de Rankine–Hugoniot com velocidade nula do sistema (72).

A **T2** (G_l, G_g e p contínuos) implica J_M = F₃,R − F₃,L ≠ 0 quando α_L ≠ α_R, porque
F₃ − p = G_l²/(ρl(1−α)) + G_g²/(ρg α) muda com α. Essa força não pertence a (72)–(75).

A perda de **Borda–Carnot** numa expansão é derivada justamente da conservação do momento
(T1). A dissipação aparece na energia mecânica, não como salto de F₃.

### 2.3 Estudo local independente com Chexal

Ferramenta: [junction_riemann_ref.f90](../validation/diagnostics/junction_constitutive/junction_riemann_ref.f90).
- Problema de Riemann com camada estacionária e curvas integrais dos autovetores.
- Newton em (p_L*, ln α_R*, p_M).
- Fechamento do riser via `PRIMS_PROBE_CLOSURE_RISER=5` no `offline_api`; o default 3
  preserva as ferramentas históricas.
- Fixtures de estados **aceitos** pelo solver Chexal
  ([make_probe_fixture.py](../validation/diagnostics/chexal_riser_20261121/make_probe_fixture.py)).
- Novo modo **MULTISTART** com 75 partidas por condição: uma falha de Newton não é tomada
  como inexistência.

Controles com o riser Chexal (modo TESTS):
- estado uniforme: exato;
- partição contra Roe: 1e-6 a 2e-4 (segunda ordem do salto);
- ondas pequenas: erro ∝ ε, de 1e-2 a 1e-4, nas três famílias.

Resultados ([P02_JUNCTION_BALANCE.md](../validation/diagnostics/chexal_riser_20261121/evidence/p02/P02_JUNCTION_BALANCE.md)):

| Estado aceito | G_g* T1 | G_g* T2 | G_g* Roe | T1: p_R* − p_L* | T2: J_M | \|dp/ds\| do riser | l_eq = \|J_M\|/\|dp/ds\| | Multistart (convergidas / partidas, G_g* distintos) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| OSC 60 s | 1,300 | 1,281 | 1,285 | +193 Pa | −195 Pa | 1378 Pa/m | 0,14 m = 2,7 D | T1 15/75 {1,300}; T2 15/75 {1,281} |
| OSC 135 s | 1,993 | 1,979 | 1,985 | +129 | −130 | 1387 | 1,8 D | — |
| OSC (último aceito) | 0,618 | 0,599 | 0,605 | +267 | −270 | 1367 | 3,8 D | — |
| SS3 30 s | 0,850 | 0,835 | 0,839 | +167 | −168 | 1407 | 2,3 D | — |
| SS3 (último aceito) | 0,613 | 0,599 | 0,601 | +175 | −175 | 1223 | 2,7 D | — |
| Stable 600 s | 0,851 | 0,826 | 0,831 | +289 | −291 | 1650 | 3,4 D | — |
| SS1 30 s | 0,176 | 0,162 | 0,166 | +326 | −328 | 1869 | 3,3 D | — |
| **SS1 (bloqueio, 52,94 s)** | **2,94e-3** | **sem solução** | **3,1e-10** | +355 | — | 2015 | — | **T1 65/75 {2,944e-3}; T2 0/75: 17 atraídas a α_R* → 0 com G_g* ≈ −1e-11 a −5e-9** |

(G_g* em kg/(m²·s).)

Leitura:
1. O operador atual é T2 (≤ 0,6 %).
2. T1 e T2 diferem 1–9 % em G_g* nos estados regulares, e **qualitativamente** no bloqueio
   do SS1: sob T2, a única saída é G_g* ≤ 0, a inversão de gás que o Chexal horizontal recusa
   (3061). Sob T1 existe uma solução única com gás passando.
3. A força pontual da T2 (130–330 Pa) equivale à força distribuída do próprio modelo em
   2–4 D. Uma região real de desenvolvimento de alguns diâmetros poderia produzi-la, mas
   essas forças já estão em S₃ das células. Adicioná-la como ponto seria contá-la sem base
   constitutiva própria.
4. Hiperbolicidade: as soluções T1/T2 convergidas satisfazem a configuração de ondas H4
   (λ₁ < 0 < λ₂ < λ₃; não houve status 60).
5. Limites das hipóteses:
   - H1: fontes desligadas no instante t → 0⁺;
   - H3: curvas integrais, válidas para ondas fracas; choques fortes ficam fora;
   - os iterados podem cruzar quinas do Chexal (status 153 = 100 + 53), e isso é reportado
     como falha de avaliação, não como inexistência.

**Gate P0.2: PASS (recomendação T1)**, com as hipóteses H1–H4 explícitas. A implementação
global (F₃ comum na face MATCHED; p_R* do traço do riser obtido de F₃ − cinética(α_R*)) é um
**passo seguinte** que exige Jacobiano, testes de partição e conservação, e replays curtos
dedicados. Nada foi alterado no operador global nesta rodada.

**Experimento discriminante (proposta):**
- Escoamento permanente na junção oleoduto–riser, com tomadas de pressão a −10 D e
  +2…10 D da base e holdup medido nas duas seções.
- A T1 prevê, além do gradiente distribuído, um **aumento** de pressão de
  ΔK = G_l²/ρl·[1/(1−α_L) − 1/(1−α_R)] (130–360 Pa nas condições acima). A T2 prevê salto nulo.
- A instalação de Baliño et al. (2010) tem sensores de pressão na base do riser; a separação
  entre esses dois sinais exige resolução de ~50 Pa.

## 3. Banco experimental ar–água (P1.1)

Registro rastreável: [validation/experimental_airwater/airwater_void_fraction_sources.csv](../validation/experimental_airwater/airwater_void_fraction_sources.csv)
(e `.json`), com fonte, experimento, D, θ, P, T, fluidos, j_l, j_g, α_exp, incerteza,
regime, tipo de dado, disponibilidade e relevância. A **disponibilidade foi verificada
nesta sessão**:

| Fonte | Tipo | Disponibilidade |
|---|---|---|
| Beggs (1972), via Chexal 1992 Tab. 5-1 | primária, só o resumo (621 pontos; erro médio 0,003; desvio 0,038) | pontos indisponíveis |
| Ryan & Kim (2025), NSE 199:1614–1625, doi:10.1080/00295639.2024.2385777 | primária (bolhas, inclinado) | não acessada (editor) |
| Ryan (2022), tese Purdue | primária | repositório recusou acesso automatizado (HTTP 403) |
| Ryan et al. (2024), Geoenergy Sci. Eng. | primária | não acessada |
| Dong & Hibiki (descendente inclinado; ano a confirmar, ~2023) | **compilação secundária** (2868 pontos, 12 fontes) | não acessada; não é medição primária |
| Bhagwat & Ghajar (OSU) | primária (D = 12,7 mm) | só resumos |
| Ousaka, Deendarlianto, Kariyasaki & Fukano (2006), NED 236:1282–1292 | primária, **CCFL** (16 e 26 mm; 30–60°) | **acessada** (PDF aberto); só figuras e correlação |
| Ohnuki & Adachi (1987), Trans. JSME B 53:1685–1690 | primária, CCFL horizontal + inclinado | só resumo |
| Barnea, Ben Yoseph & Taitel (1986) | primária, CCFL inclinado | não acessada |

Não há pontos de α tabulados acessíveis. Nada foi digitalizado nem inventado, e nenhum erro,
viés ou RMSE foi calculado.

A cobertura dos estados do PRIMS-F pela base da própria fonte (Parte II, §II.5) é: 0 %
validado (D fora dos dados inclinados de Beggs), 99–100 % interpolação por Fr e 0,8 %
extrapolação (SS1).

**Gate P1.1: NOT_VERIFIED.** A matriz de comparação está preparada (estrutura e
proveniência), mas os dados primários estão indisponíveis.

## 4. Mapa constitutivo da contracorrente inclinada (P1.2)

| Mecanismo | Onde | Candidato | Base | Domínio | Aplicável? |
|---|---|---|---|---|---|
| **A**: u_g < 0, u_l > 0 | traço horizontal da junção (SS1) | Chexal 1992 | horizontal só cocorrente (pp. 264–265) | — | não |
| A | idem | Ohnuki & Adachi 1987 | CCFL horizontal + inclinado | não verificado (texto completo indisponível) | não avaliável |
| A | idem | **T1 (§2)** | conservação do momento | — | **sob T1, o estado de bloqueio do SS1 não exige reversão de gás** (solução única com G_g* > 0) |
| **B**: u_l < 0, u_g > 0 | catenária 65–76° (OSC, SS3) | Chexal 1992 | contracorrente só vertical (p. 262); sem dados fora da CCFL (p. 270) | vertical | não (só sensibilidade) |
| B | idem | Ousaka et al. 2006 | início de inundação; We com altura de onda **medida** | D 16 e 26 mm; θ 30–60°; 1 atm; J_L 0,03–0,32; ±35 % | não: é limite de inundação, a correlação não é fechada, e D e θ estão fora |
| B | idem | Barnea et al. 1986 | CCFL inclinado | não acessado | não avaliável |
| B | idem | Bendiksen (fechamento 2 do repositório) | bolhas longas, inclui líquido estagnado | cocorrente | só a borda j_l = 0 |
| B | idem | Bhagwat (fechamento 3) | definido até j_l = 0 | cocorrente | só a borda |

Solver local ([inclined_cc_local.py](../validation/diagnostics/chexal_riser_20261121/inclined_cc_local.py)),
com velocidades assinadas, nos 11 estados reais de bloqueio:
- na borda j_l = 0, os candidatos diferem 2,4–3,5× em u_g;
- com j_l < 0, diferem 1,4–1,8×;
- o Chexal vertical tem inclinação infinita em j_l = 0⁻ (termo |Re_l|^0,03 da Eq. 3-36).

**Gate P1.2: BLOCKED.**
- **Mecanismo B:** não há formulação experimentalmente fundamentada para j_l → 0⁻ a 65–76°
  com D ≈ 52,5 mm. Os dados acessíveis tratam do início de inundação, em diâmetros e
  ângulos incompatíveis.
- **Mecanismo A:** depende da condição de transmissão; sob T1 deixa de ser exigido no
  estado testado.

Nenhuma física nova foi integrada.

## 5. Regressões e conservação

**Builds e testes**, com GNU Fortran 13.3.0, `-std=f2018`, LAPACK/BLAS do sistema e uma thread:
- release e `check_all` (`-fcheck=all -ffpe-trap`) sem erro;
- Fortran **29/29** (inclui `test_chexal_domain_contract`, `test_chexal_np3989`,
  `test_chexal_angles` e a partição com Chexal);
- Python `unittest`: OK;
- CLI **7/7**: `pipe_regime`, `accepted_telemetry`, `feed_contract`, `boundary_telemetry`,
  `initial_inventory`, `chexal_default` e `closure_contract`.

**Regressão do contrato de domínio:**
- Binário `prims_f_domain`, SHA-256 `0f0e82aafab52caf030760e0d24fcb8527776eb45bd67c46813292b2f49546d6`
  (fonte = `scr/` do commit `92d2d4f`).
- Os 16 decks primários (Stable, OSC, SS1, SS3 × Chexal/Bhagwat × BARO/ISO, `tol_res = 1e-10`)
  foram reexecutados e comparados com a campanha P2 da rodada anterior.
- São **idênticos byte a byte**: monitor, campos, conservação, contornos, propostas rejeitadas,
  estado rejeitado, checkpoint, `steps` sem tempo de relógio e `closure_contract.txt` sem o
  caminho do arquivo de entrada.
- Contrato `CLOSURE_CONTRACT_STATUS` PASS 16/16
  ([P3_VS_P2_BYTEWISE.json](../validation/diagnostics/chexal_riser_20261121/evidence/P3_VS_P2_BYTEWISE.json)).
- O contrato de domínio **não altera** nenhuma trajetória do PRIMS-F, e o caminho Bhagwat
  legado continua intacto. Na Parte I ele já era idêntico ao binário congelado `364b808`.

**Seis balanços** (duas fases × sistema, oleoduto e riser, com fluxos aceitos): idênticos aos
da P2. Máx. |balanço|/M₀ de 2,7e-12 a 3,0e-7 (critério 1e-6), aritmética do ledger PASS e
integral na junção sem descasamento.

## 6. Resultados globais após os gates

Os gates que condicionam campanhas novas não estão abertos:
- a T1 é recomendada mas não implementada;
- a P1.2 está BLOCKED;
- a P1.1 está NOT_VERIFIED.

Por isso não se executaram campanhas longas novas. Valem as campanhas anteriores, certificadas
pelo contrato de fechamento e idênticas byte a byte (Parte II, §II.8): Stable completo e
estável com Δt 0,5 e 0,125; OSC e SS3 `BLOCKED_CONSTITUTIVE_DOMAIN`; SS1
`BLOCKED_JUNCTION_TRACE_DOMAIN`, que este estudo atribui à convenção T2. As réplicas desta
rodada servem só de regressão do contrato de domínio (§5).

## 7. Comparação experimental

Baliño, Burr e Nemoto (2010): o PDF não está nesta máquina, e nenhuma trajetória Chexal é
completa, admissível e convergida. **NOT_VERIFIED.**

## 8. Matriz de status

| Item | Status | Evidência |
|---|---|---|
| P0.1 contrato executável de domínio ar–água | PASS | `test_chexal_domain_contract` (19); regressão §5 |
| P0.1 matriz ramo × domínio × evidência | PASS (documentada) | §1.1 |
| P0.1 EPRI 1000 psia (linhas 2, 6, 7) | **FAIL** para uso generalizado; divergência localizada em propriedades | §1.2 |
| P0.1 EPRI linha 3 | FAIL (referência inconsistente) | §1.2 |
| P0.2 recomendação T1 × T2 | **PASS: T1** | §2 |
| P0.2 implementação global de T1 | NOT_VERIFIED (não implementada) | §2.3 |
| P1.1 banco experimental | NOT_VERIFIED (dados primários indisponíveis) | §3 |
| P1.2 mecanismo A | depende de T1; não exigido sob T1 no estado testado | §4 |
| P1.2 mecanismo B | BLOCKED | §4 |
| P2 campanhas globais novas | BLOCKED pelos gates | §6 |
| P2 regressões | PASS (29/29 Fortran, Python OK, 7/7 CLI; 16 réplicas idênticas byte a byte, contrato 16/16) | §5 |
| Comparação com Baliño et al. 2010 | NOT_VERIFIED | §7 |
| SS1, SS3, OSC validados | não | — |

## 9. Commits, hashes e comandos

Commits locais desta rodada (branch `implementation/chexal-riser-20261121`, sem push nem merge; `main` = `f3717dd`):

| SHA | Conteúdo |
|---|---|
| `92d2d4f` | feat(core): contrato executável ar–água (`CHEXAL_UNQUALIFIED` = 68) e teste |
| `845325f` | diag: referência T1 × T2 com riser Chexal (fixtures, multistart, balanço) |
| `240fa1f` | diag: decomposição NP-3989; registro de fontes experimentais ar–água |
| (seguintes) | evidence (P0.2, regressão) e docs (este relatório e handoff) |

Binários: `prims_f_domain` `0f0e82aa…` (esta rodada) e `prims_f_contract` `9e54e79f…` (rodada 2).

```sh
make release && make check_all
PYTHONPATH=scripts python3 scripts/run_unit_tests.py
D=validation/diagnostics/chexal_riser_20261121
python3 $D/np3989_decompose.py OUT.json OUT.md
# T1 x T2 with the Chexal riser
python3 $D/make_probe_fixture.py RESULTS/p2/runs/D_SS1_C5_BAROTROPIC_N51_dt0.5 DUMP FIX/D_SS1_DUMP
gfortran ... offline_api.f90 junction_riemann_ref.f90 (README de junction_constitutive)
PRIMS_PROBE_CLOSURE_RISER=5 ./tmp/jc5/junction_riemann_ref FIX/D_SS1_DUMP 343 JUNCTION
PRIMS_PROBE_CLOSURE_RISER=5 ./tmp/jc5/junction_riemann_ref FIX/D_SS1_DUMP 343 MULTISTART
python3 $D/p02_junction_balance.py FIX $D/evidence/p02 OUT.json OUT.md
# regression of the domain contract
python3 $D/run_campaign.py BIN RESULTS/p3_regress --workers 6
```

Custos (CPU de uma thread):
- referência de Riemann: segundos por estado; multistart de 150 resoluções: ~10 min no SS1 e ~15 min no OSC;
- regressão: 87–547 s por execução;
- auditoria NP-3989: 3 min.

## 10. Próximos experimentos e passos discriminantes

1. **Operador T1 experimental (opt-in):**
   - na face MATCHED, F₃ comum; p_R* do traço do riser = F₃ − G_l²/(ρl(1−α_R*)) − G_g²/(ρg α_R*),
     com α_R* do Chexal em p_R* (ponto fixo);
   - Jacobiano por diferenças da face, testes de partição e conservação;
   - replay curto do SS1 a partir do estado aceito de 52,9 s.
   
   Critério: G_g* > 0 sem recusa de domínio, conservação exata e nenhuma mudança no Stable.
   Só depois disso, campanhas.
2. **Medição de pressão na junção** (§2.3): a T1 prevê aumento de 130–360 Pa, a T2 salto nulo.
3. **Dados de contracorrente inclinada** com líquido quase estagnado: α, u_g e u_l a 65–76°
   com D ≈ 50 mm (holdup e velocidades de fase), não só o início de inundação.
4. **Acesso a Ryan (2022, 2025)**, que provavelmente está mais próximo em D; a verificar:
   tabelas de α por ângulo para comparar o Chexal cocorrente inclinado.
5. **PDF de Baliño et al. (2010)** para a comparação experimental, quando houver trajetórias
   completas.
