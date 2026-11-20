# Continuação 2 da auditoria físico‑numérica — 20/11/2026

Branch local `implementation/p0p1-continuation-20261119` (worktree
`~/Documentos/PRIMS-F-p0p1c-20261119`), sobre `f3717dd`. `main`, histórico, PDFs e
dados preexistentes preservados. **Nenhum push ou merge.** Dados brutos em
`~/Documentos/PRIMS-F-p0p1c-results-20261119` (fora do Git). Relatório anterior:
`docs/QUALIFICATION_P0P1_CONTINUATION_20261120.md`.

**Conclusão.**
- **P0.1 corrigido.** A face do separador do modo reduzido passou a ser
  *well‑balanced* com G≠0. O repouso hidrostático discreto é preservado bit a bit, e o
  refluxo real sem composição declarada continua recusado.
- **P0.2.** Os três bloqueios têm certificado determinístico. Em todos, o
  bloqueio está numa **proposta de Newton** a partir de estados aceitos admissíveis.
  - SS1 e OSC‑ISO: limite do domínio constitutivo da face de junção, induzido pelo
    fluxo numérico, com as duas células ainda transportando gás para cima.
  - SS3: ausência de raiz no ramo cocorrente declarado.
- **P1.1** fica **BLOCKED por limitação documentada da formulação**: o traço de face
  não representa a passagem misto→puro. A extensão está proposta, não implementada.
- **P1.2 BLOCKED.** Nenhum fechamento qualificado admite o estado contracorrente
  inclinado de SS3. A comparação experimental é **NOT_VERIFIED**.
- **P1.3 ISOTHERMAL** (9 discretizações, N51–N401, dt 0,125–0,015625).
  - O período converge (≤ 0,2%).
  - A taxa de crescimento da instabilidade está na região assintótica, com ordem
    ≈ 1.
  - O mapa amplitude→vale do fluxo gasoso colapsa numa curva (±2,3%).
  - O instante do evento **não** converge e salta por ciclos inteiros, porque a
    amplitude‑semente é O(h^0,8).
- **BAROTROPIC** (as três discretizações acomodam nos 14 canais v2).
  - Período e média convergem em ordem ≈ 1.
  - A semi‑amplitude **não** atinge 1%: diferenças de 22% em h e de 28% → 10% em dt,
    com erro estimado de 6,8% no dt mais fino.
  - O refino combinado ainda está em execução: **NOT_VERIFIED**.
- Nenhum regime de *severe slugging* foi validado.

## 1. Ambiente e manifesto

| Item | Valor |
|---|---|
| Máquina | `leonardo-ubuntu`, Linux 7.0.0-42, 8 núcleos |
| Compilador | GNU Fortran 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) |
| Flags release | `-std=f2018 -Wall -Wextra -Wno-unused-parameter -Wno-unused-dummy-argument -fimplicit-none -ffree-line-length-none -fopenmp -Wa,--noexecstack -O2` |
| Flags `check_all` | as anteriores com `-O0 -g3 -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-999` |
| BLAS/LAPACK | `-llapack -lblas` → alternativa do sistema → **Intel MKL 2020.4.304** (`libmkl_rt.so`); `liblapack3`/`libblas3` 3.12.0 instalados, não ativos |
| Threads | `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=1` nas campanhas |
| Python | 3.12.3, NumPy 1.26.4, SciPy 1.11.4 |
| SHA congelado (regressões) | `1d71f00` → `frozen_1d71f00/`: `prims_f` `caee88e9…`, driver barotrópico `f276ff96…` |
| Binários de campanha | ISO: `prims_f_aef3774` `e33a1b8d…`; BARO: `restart_baro_813da75` `ef5c87e6…`; ambos **byte a byte equivalentes** a `1d71f00` nos decks usados (§8) |
| IC barotrópica | `common_initial_U.txt` `192c5f7e…` (projeção conservativa comum; mesma EOS) |
| Índice de execuções | `validation/diagnostics/p0p1c_20261119/evidence/RUN_INDEX.json` (hashes, deck, status, custo, balanços) |

## 2. Commits desta continuação

| Commit | Conteúdo |
|---|---|
| `0960246` | P0.1: face do separador *well‑balanced* para o repouso hidrostático com gravidade; `test_hydrostatic_outlet` |
| `1dde473` | certificados SS1/OSC‑ISO; driver da campanha barotrópica congelada |
| `9c31dea` | certificado SS3; estudo constitutivo local em velocidades assinadas; contrato da face com gás evanescente |
| `1d71f00` | análise ISOTHERMAL em janela comum e eventos por grandezas contínuas |
| `1518f5d` | evidência ISOTHERMAL N51–N401 e dt 0,125–0,015625; emenda registrada do ajuste de σ |
| `21ad166` | análise BAROTROPIC v2 (espaço, tempo e combinado) |
| (este commit) | evidências BARO/RUN_INDEX, relatório e handoff |

## 3. Matriz problema → causa → correção → teste → estado

| # | Problema | Causa demonstrada | Correção | Teste | Estado |
|---|---|---|---|---|---|
| P0.1 | repouso hidrostático com G≠0 recusado (`INITIAL_OUTLET 934`) | a densidade de face do separador vinha de `ρ_N·f_r`, e o fator de reconstrução arredondado dava `v_b = ±c·O(ε)` no repouso exato; quando negativo, virava refluxo aritmético | forma de **desvio** em relação ao equilíbrio discreto da própria célula (§4); sem tolerância no sinal | `test_hydrostatic_outlet`: 160 configurações, **1760 asserções**; `test_reduced_equilibrium_cli` (G=9,81, repouso líquido/gás) | **PASS** |
| P0.2 | bloqueios SS1, SS3 e OSC‑ISO sem classe determinística | — | certificador `certificate.f90` + `make_certificate.py`; classes do pedido, com regras e evidência registradas | certificados JSON/MD reprodutíveis | **PASS** (diagnóstico); casos BLOCKED |
| P1.1 | passagem misto→puro na junção | o traço do riser só existe para `j_g>0`; o fluxo Roe de face zera `G_g` com as duas células a `u_g>0` | **não implementada**: extensão proposta (§6) | `test_junction_vanishing_gas` (13 asserções) fixa o contrato atual | **BLOCKED** (limitação documentada) |
| P1.2 | SS3 contracorrente inclinado | o estado exige `u_l<0` a 51,6°; nenhum fechamento qualificado cobre esse ramo | solver local autônomo em velocidades assinadas (diagnóstico, não integrado) | `ss3_local_closures.f90` | **BLOCKED**; experimento **NOT_VERIFIED** |
| P1.3a | evento ISO não converge | crescimento exponencial (σ convergente, ordem ≈ 1) a partir de semente O(h^0,8); evento quantizado por ciclo | métricas pré‑registradas: janela comum, eventos contínuos, evento no espaço de amplitude | `iso_event_analysis.py` (9 discretizações, N51–N401, dt 0,125–0,015625) | período **PASS**; σ assintótico (1,2–1,5%); A* ±2,3%; instante **FAIL** |
| P1.3b | convergência BAROTROPIC | amplitude dominada pela dissipação numérica (ordem ≈ 1,3 em dt) | nenhuma: não se ajusta física para obter convergência | `run_baro_campaign.py` + `baro_refinement_v2.py` (14 canais) | acomoda **PASS**; amplitude **FAIL**; combinado **NOT_VERIFIED** |

## 4. P0.1 — saída hidrostática *well‑balanced* (modo reduzido)

Equilíbrio discreto da célula de saída N (fase k, gravidade G, cota z nos pontos de
Gauss da célula):

  ρ_eq,N = ρ_t · Σ_q w_q exp(−G (z_q − z_ref)/(ω_u² c_k²)) / Σ_q w_q

Esta é a mesma média normalizada da IC (`initialize_single_phase_hydrostatic`), de
modo que `u_N = ρ_eq,N` bit a bit no repouso. Desvio relativo:
`x = (u_N − ρ_eq,N)/ρ_eq,N`.

As relações de Riemann da saída subsônica, antes `v_b = v + c ln(ρ_f/ρ_b)` ou
`v − c(ρ_b−ρ_f)/√(ρ_bρ_f)`, passam a ser avaliadas na variável de desvio:

- `x ≥ 0`: v_b = v + c·log1p(x);
- `x < 0`: v_b = v + c·x/√(1+x).

`log1p` é implementado pelo truque de Goldberg (1991): `w = 1+x`; se `w = 1`,
devolve `x`, senão `ln(w)·x/(w−1)`. No repouso, `x = 0` exatamente, logo `v_b = 0`.
Com perturbação, o sinal de `v_b` é o da perturbação física, e a equivalência com a
forma legada vale a 10⁻⁶ relativo em todas as configurações.

Nenhuma tolerância foi introduzida no sinal da velocidade. Com refluxo real
(`x < 0`) e `outlet_reservoir_phase = 0`, o resultado continua sendo
`PHASE_DATA_REQUIRED (934)`.

Fontes (bibliográficas; páginas não conferidas):
- Käppeli & Mishra (2014), JCP 259, reconstrução por desvio de um equilíbrio
  hidrostático;
- Berberich, Chandrashekar & Klingenberg (2021), Comput. Fluids 219, métodos
  *well‑balanced* de alta ordem;
- Goldberg (1991), ACM Comput. Surv. 23(1), avaliação estável de `log1p`.

| Grade de teste | Antes (`7ba50db`) | Depois (`0960246`) |
|---|---|---|
| Fase líquido/gás × G ∈ {0; 4,905; 9,81; 19,62; 98,1} × N ∈ {4, 8, 16, 64} × geometrias (catenária; oleoduto reto ascendente β=−30°) × escalas ω ∈ {1; 3,7} | **FAIL** na 1ª configuração com G≠0 (líquido, G=4,905, N=4, catenária): repouso recusado | 160/160 repousos admissíveis; `v_b = 0` exato; Jacobiano completo; sinais ± corretos; refluxo sem composição → 934; **1760 asserções PASS** |

## 5. P0.2 — certificados determinísticos

`certificate.f90` reconstrói, no **mesmo operador**, cada estágio do último passo:
- estado aceito e última proposta;
- resíduo e Jacobiano: `dgetrf` + `dgecon`, condição em norma 1, e verificação por
  diferenças finitas centrais ou unilaterais;
- raio de Newton: primeira inadmissibilidade, mapa de raízes e margens;
- traços e fluxo de face da junção, decomposto em convectivo e dissipativo.

`make_certificate.py` aplica as regras de classificação. Os arquivos ficam em
`evidence/certificates/{SS1,OSC_ISOTHERMAL,SS3}.{json,md}`.

| Caso | Bloqueio em | Local | κ₁(J) / erro DF | Raio | Classes | Indicador |
|---|---|---|---|---|---|---|
| SS1 (78,5766 s) | `NEWTON_PROPOSAL` (aceito e tentativa admissíveis) | face 53, s = 56,84 m, −2° | 41 / ≤ 8,8·10⁻⁷ | λ* = 1,59·10⁻⁹; sem descida (razão 0,99992) | `CONSTITUTIVE_DOMAIN_LIMIT` + `UNSUPPORTED_PHYSICAL_BOUNDARY` | `NUMERICAL_FLUX_INDUCED_BOUNDARY` |
| OSC‑ISO (223,6445 s) | `NEWTON_PROPOSAL` | face 53 | 64 | λ* = 3,37·10⁻⁸ | idem | idem |
| SS3 (97,1318522 s) | `NEWTON_PROPOSAL` | célula 69, s = 61,21 m, 51,64° | 1,84 | λ* = 3,68·10⁻⁶ | `NO_ROOT_IN_DECLARED_BRANCH` + `CONSTITUTIVE_DOMAIN_LIMIT` | raiz na extremidade `u_l = 0`; margem de *drift* 9,6·10⁻¹³ já no aceito |

Classes descartadas com evidência:
- `JACOBIAN_FAILURE`: κ pequeno e Jacobiano conferido por diferenças finitas.
- `GLOBALIZATION_FAILURE`: não há descida na parte admissível do raio.
- `ROOT_SEARCH_FAILURE`: o mapa de raízes é exaustivo, com bisseção e contagem das
  recusas.
- `INADMISSIBLE_STATE`: os estados aceitos são admissíveis.

Face de junção no estado aceito, em kg/(m²·s):

| Caso | G_g convectivo | G_g dissipativo | G_g face | traço do oleoduto | traço do riser | α* do riser na face |
|---|---:|---:|---:|---:|---:|---:|
| SS1 | +0,055762 | −0,055762 | 2,6·10⁻⁸ | 0,0711 | 0,0404 | 1,2·10⁻⁸ |
| OSC‑ISO | +0,058485 | −0,058485 | 1,1·10⁻⁸ | 0,1003 | 0,0166 | 3,7·10⁻⁹ |

A dissipação Roe cancela o fluxo centrado a 10⁻⁸ relativo, levando G_g para **fora**
do intervalo dos fluxos dos dois traços, que transportam gás para cima. A
proposta seguinte pede G_g < 0, o que leva ao fechamento cocorrente do riser fora do
domínio (`3012 = FACE 2000 + 1000 + CLOSURE_DOMAIN 12`). Por isso o certificado não
conclui bloqueio físico: o limite é induzido pelo fluxo numérico. A partição
característica já registrada atribui à diferença de pressão de equilíbrio entre os
centros 61% (SS1) e 84% (OSC) da dissipação gasosa.

## 6. P1.1 — transição local misto→puro na junção

**Contrato atual** (`test_junction_vanishing_gas`: fechamento de Bhagwat no riser,
θ = 0, líquido ascendente):

| j_g de face | status | α* |
|---|---|---|
| 10⁻⁴ … 10⁻¹⁴ m/s | 0 (admissível) | α* = 1,0935·j_g [s/m], limite *drift‑flux* contínuo, razão estável a 10⁻³ |
| 0 | **21** (recusa) | — |
| < 0 | **12** (`CLOSURE_DOMAIN`) | — |

O limite α*→0⁺ existe e é contínuo, mas o fluxo nulo não seleciona o traço
monofásico líquido. Um fluxo gasoso de face negativo sob líquido ascendente está
fora do fechamento cocorrente.

A formulação atual **não permite** a passagem. Não foram introduzidos `j_g = 0`
imposto, bloqueio de gás, *clipping* ou troca de fechamento.

Literatura (bibliográfica; texto integral não lido, salvo indicação):
- Munkejord, Evje & Flåtten (2006): o MUSTA trata a transição para escoamento
  monofásico com frações positivas.
- Flåtten & Munkejord (2006, ESAIM M2AN): o Roe para *drift‑flux* exige correção na
  fase evanescente.
- Evje & Fjelde (2002, JCP): AUSMV, com fluxos de massa por fase contra a corrente
  (*upwind*) e positividade.

**Extensão proposta** (não implementada; exige aprovação):

E1. Fluxo de massa por fase na face de junção, *upwind* pela velocidade de fase
(tipo AUSMV), com o termo de pressão separado. Propriedades a demonstrar antes da
promoção:
- fluxo único por fase, compartilhado com sinais opostos;
- `G_g` no intervalo dos fluxos dos traços quando `u_g > 0` dos dois lados;
- equilíbrio estacionário preservado, reutilizando `junction_trace = 'STEADY'`;
- Jacobiano verificado por diferenças finitas;
- default `CENTER`/Roe byte a byte inalterado.

E2. Limite de face α*→0⁺ **somente** como continuação da raiz existente (α* = 0 quando
o fluxo calculado é exatamente zero e o traço adjacente é monofásico). Nunca como
regra que imponha `j_g = 0`.

**Experimento discriminante.**

1. Avaliação *offline* nos estados aceitos terminais. O E1 dá `G_g = 0,0711` (SS1) e
   `0,1003 kg/(m²·s)` (OSC), positivos e dentro dos limites. O Roe dá 2,6·10⁻⁸ e
   1,1·10⁻⁸.
2. Implementação opcional (*opt‑in*), seguida de OSC‑ISO N51/N101/N201. Se o mapa
   A*(f) do §7 mudar e o vale do fluxo gasoso deixar de ir a zero enquanto a célula
   do fim do oleoduto mantém `u_g > 0`, o bloqueio atual é artefato numérico. Se o
   vale persistir com `u_g → 0` na própria célula, é selo físico, e a passagem
   misto→puro (E2) passa a ser necessária.

## 7. P1.2 — SS3: fechamento contracorrente inclinado

Estado local (célula 69, aceito):

| Grandeza | Valor |
|---|---|
| α | 0,638 |
| p | 228 440 Pa |
| θ | 51,64° |
| G_mix | 0,774 kg/(m²·s) |

O solver local autônomo (`ss3_local_closures.f90`) fixa U₁, U₂, U₃ e a EOS e
varre `u_g ∈ [0, 2·U₃/U₂]`, com `u_l = (U₃ − U₂u_g)/U₁` assinado. Cobre assim o ramo
cocorrente e a sua extensão contracorrente.

| Fechamento | Domínio declarado | Raízes | Resultado |
|---|---|---|---|
| Estratificado de dois fluidos (id 1) | estratificado, inclinação baixa | **0** | não há equilíbrio estratificado neste estado |
| Bendiksen 1984 (id 2) | golfadas, cocorrente | 1, contracorrente (`u_l = −0,0037`, `u_g = 1,71 m/s`) | **extrapolação** fora do domínio; não qualificável |
| Bhagwat (id 3) | cocorrente | 0 | 4905/6001 avaliações recusadas pelo domínio |
| Chexal et al. 1992 | contracorrente **apenas vertical** | não aplicado | PDF local, p. 8 (p. 262 do periódico): "Countercurrent flow is only considered for vertical flows"; p. 16 (p. 270): a região contracorrente fora da linha CCFL "is not yet verifiable" |

Fontes indicadas no pedido para fechamento contracorrente inclinado (bibliográficas;
textos integrais **não** disponíveis localmente nem lidos):
- Ullmann et al. (2004); referência completa não conferida.
- Ullmann & Brauner (2006), IJMF 32:82–105: fechamentos de dois fluidos para
  escoamento estratificado.
- Barnea, Ben Yoseph & Taitel (1986), Can. J. Chem. Eng. 64:177–184: *flooding* em
  tubos inclinados de 1 a 90°.
- Ohnuki, Adachi & Murao (1988), Nucl. Eng. Des. 107:283–294: contracorrente em tubo
  horizontal com *riser* inclinado.

Também não há dados experimentais inclinados contracorrente no repositório.

Estado: SS3 **BLOCKED**; nenhum fechamento integrado. A comparação com experimento
fica **NOT_VERIFIED**. Próximo passo: obter os textos e os dados (Barnea et al.
1986; Ohnuki et al. 1988), implementar o fechamento estratificado contracorrente de
Ullmann & Brauner no solver local e só integrá‑lo após reproduzir esses pontos.

## 8. P1.3 — OSC ISOTHERMAL: janela comum e eventos

O mesmo deck e o mesmo binário foram usados nas 9 execuções:
- N ∈ {51, 101, 201, 401} com dt_max 0,125 s;
- dt_max ∈ {0,125; 0,0625; 0,03125; 0,015625} s com N51;
- o combinado N101/dt 0,0625;
- MASS × LOCAL_VOLUME.

Nenhum estado de outra EOS foi reutilizado.

Definições pré‑registradas na *docstring* de `iso_event_analysis.py` (`1d71f00`),
fixadas antes de N201, N401, dt 0,015625 e N101/dt 0,0625:
- janela comum `[0, T_c]`, com `T_c` = 197,83 s, a menor duração alcançada;
- E_g(f): primeiro instante após 20 s com `ṁ_g,junção ≤ f·ṁ_g,in`;
- E_h(h): primeiro instante com holdup líquido da base ≥ h;
- t_fail: instante da falha de Newton, separado dos eventos;
- σ: inclinação de ln A_k contra t;
- A*(f): amplitude de p_in em que o vale de ṁ_g,junção por ciclo atinge f.

**Emenda registrada:** a regra de σ ("ciclos após 60 s") foi trocada por "ciclos
iniciados em (60, 150) s" depois que as duas últimas execuções encurtaram a janela.
A regra original ajustava 7 ciclos numa execução e 6 nas outras; a emenda usa os
mesmos quatro índices de ciclo em todas.

Evidência em `validation/diagnostics/p0p1c_20261119/evidence/iso_events/`:
`ISO_EVENTS.{json,md}`, `iso_temporal.png`, `iso_growth.png` e `iso_spatial.png`.

| Execução | t_fail (s) | E_g 0,25 | E_g 0,01 | E_h 0,90 | σ (1/s) | período (s) | A₁ (Pa) | A*(0,5) | A*(0,25) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| N51 dt 0,125 | 223,645 | 197,17 | 223,41 | 197,78 | 0,02513 | 26,510 | 229 | 6712 | 11594 |
| N51 dt 0,125 LOCAL_VOLUME | 223,642 | 197,18 | 223,40 | 197,79 | 0,02513 | 26,509 | 228 | 6710 | 11591 |
| N51 dt 0,0625 | 222,518 | 194,96 | 222,38 | 195,45 | 0,02667 | 26,533 | 240 | 6420 | 11215 |
| N51 dt 0,03125 | 198,667 | 194,40 | 198,32 | 194,89 | 0,02742 | 26,558 | 247 | 6347 | 11174 |
| N51 dt 0,015625 | 197,831 | 194,18 | 197,59 | 194,68 | 0,02778 | 26,569 | 250 | 6335 | — |
| N101 dt 0,125 | 249,400 | 221,39 | 226,36 | 221,78 | 0,02685 | 26,561 | 131 | 6434 | 11331 |
| N201 dt 0,125 | 250,379 | 224,03 | 250,21 | 224,48 | 0,02776 | 26,607 | 74 | 6513 | 11187 |
| N401 dt 0,125 | 276,084 | 248,57 | 275,95 | 248,90 | 0,02821 | 26,639 | 43 | 6363 | 11068 |
| N101 dt 0,0625 | 223,588 | 197,23 | 223,43 | 197,67 | 0,02848 | 26,600 | 138 | 6442 | 11056 |

Todas as falhas são `MAX_NEWTON` (dt 0,0625: `NEWTON_GLOBALIZATION_FAILURE`), na
face de junção (53/103/203/403), com `3012`.

Leitura:

1. **Período.** Variação ≤ 0,19% entre níveis sucessivos em h e ≤ 0,09% em dt:
   **PASS** a 1%.
2. **Taxa de crescimento σ.**
   - dt: 0,02513 / 0,02667 / 0,02742 / 0,02778. Ordem observada 1,06 e 1,03;
     extrapolado 0,0281 s⁻¹; erro do nível mais fino 1,2%.
   - h: 0,02513 / 0,02685 / 0,02776 / 0,02821. Ordem 0,94 e 1,01; extrapolado
     0,0287 s⁻¹; erro 1,5%.
   - **Região assintótica demonstrada** (ordem ≈ 1, coerente com FV1 e Euler
     implícito), a 1,2–1,5% do alvo.
   - O estado inicial é linearmente instável, com amplificação ≈ 2,1–2,2 por ciclo
     no limite.
3. **Amplitude‑semente.** A₁ = 229 / 131 / 74 / 43 Pa em h, razão 1,75 / 1,77 / 1,72
   por refino, ou seja O(h^0,8). Em dt, quase constante (229–250 Pa).
   - A IC (`t = 0`) mostra um pico de u_g ≈ 1,8 m/s nos últimos metros do oleoduto,
     junto à junção, que o operador discreto relaxa (`iso_spatial.png`).
   - A semente vir desse desequilíbrio local é **inferência**, não demonstração.
4. **Evento no espaço de amplitude.** O vale de ṁ_g,junção por ciclo, como função
   da amplitude de p_in, colapsa numa única curva nas 9 execuções
   (`iso_growth.png`).
   - A*(0,25) = 11,06–11,59 kPa (±2,3%); A*(0,5) = 6,34–6,71 kPa (±2,9%).
   - Sem os dois níveis mais grossos, a dispersão fica em ±1,3%.
5. **Instante do evento.**
   - h: t_fail = 223,6 / 249,4 / 250,4 / 276,1 s.
   - dt: 223,6 / 222,5 / 198,7 / 197,8 s.
   - O instante muda por **ciclos inteiros** (≈ 26,5 s). O ciclo em que o vale chega
     a zero é um inteiro, determinado pela semente e por σ. Duas concordâncias
     (N101/N201 e dt 0,03125/0,015625) são coincidências de ciclo, não convergência.
   - **Evento ≠ falha.** E_g(0,25) precede t_fail de 3,5 a 28 s. Na N101, o vale de
     0,043·ṁ_in em 226 s recuperou‑se, e a falha veio no ciclo seguinte. A falha de
     Newton ocorre apenas quando o vale chega a zero na face.
   - Convergência do instante: **FAIL, não declarada.**
6. **Formas de onda em `[0, T_c]`.** RMS |Δp_in|/p̄ = 0,43 → 0,26 → 0,14% em dt e
   0,36 → 0,35 → 0,25% em h. MASS × LOCAL_VOLUME: < 10⁻⁴.

**Próximo experimento.** Inicializar cada discretização pelo seu **próprio** estado
estacionário discreto, mais uma perturbação declarada de amplitude fixa (≫ 250 Pa),
com a mesma EOS e o mesmo operador. Com σ e A* convergentes, o instante do evento
deve então convergir. Em paralelo, repetir com o fluxo de junção E1 (§6): se
A*(f) mudar, o selo depende do fluxo numérico.

## 9. P1.3 — BAROTROPIC: refino separado e combinado

Configuração:
- binário congelado `restart_baro_813da75` (`ef5c87e6…`), byte a byte equivalente ao
  HEAD (§10);
- IC comum `common_initial_U.txt`, com projeção conservativa;
- 900 s por execução, critérios v2 pré‑registrados (`813da75`), dez últimos ciclos e
  14 canais;
- script `baro_refinement_v2.py`; evidência em `evidence/baro_v2/BARO_V2.{json,md}`.

| Execução | passos aceitos / rejeitados | acomoda v2 (14 canais) | período (s) | semi‑amplitude (Pa) | média (Pa) |
|---|---:|---|---:|---:|---:|
| N51 dt 0,125 | 4822 / 0 (600 s) | sim | 23,108397 | 9905,70 | 266174,6 |
| N51 dt 0,0625 | 14418 / 0 | sim | 23,265321 | 13764,28 | 266573,6 |
| N51 dt 0,03125 | 28815 / 0 | sim | 23,344137 | 15354,38 | 266775,6 |
| N101 dt 0,125 | 7222 / 0 | sim | 23,266063 | 12762,72 | 266541,4 |
| N101 dt 0,0625 (combinado) | em execução (588 s de 900 às 14h24) | — | — | — | — |

| Par | período | semi‑amplitude | níveis (≤ 1%) | faixas (≤ 1%) | formas em fase fixa (≤ 2%) | Estado |
|---|---:|---:|---|---|---|---|
| h: N51 → N101 (dt 0,125) | 0,68% | **22,4%** | 13/14 (holdup da base: 1,04%) | 21–32% | 8,0–11,4% | **FAIL** |
| dt: 0,125 → 0,0625 (N51) | 0,67% | **28,0%** | 11/14 | 27–41% | 10–15% | **FAIL** |
| dt: 0,0625 → 0,03125 (N51) | 0,34% | **10,4%** | 14/14 | 9,7–17,7% | 3,8–6,4% | **FAIL** |
| combinado N51/0,125 → N101/0,0625 | — | — | — | — | — | **NOT_VERIFIED** (em execução) |

Ordem observada no eixo dt (três níveis; região assintótica não provada com três
pontos):

| Grandeza | Ordem | Extrapolado | Erro estimado do mais fino |
|---|---:|---:|---:|
| Período | 0,99 | 23,424 s | 0,34% |
| Semi‑amplitude | 1,28 | 16 469 Pa | **6,8%** |
| Média | 0,98 | 266 983 Pa | 0,08% |

Leitura: o ciclo limite barotrópico existe e acomoda em todas as discretizações, mas
a sua amplitude depende fortemente da dissipação numérica. Em dt 0,125, o erro de
amplitude estimado é de ~40%. Pela ordem observada, o alvo de 1% exigiria dt ≈ 0,007 s
com N51, além do refino em h, que tem efeito da mesma ordem (22%). O eixo h tem só
dois níveis: N201 e N101/dt 0,0625 tinham sido encerrados por custo na primeira
tentativa, e N101/dt 0,0625 foi relançada. **Convergência de amplitude BAROTROPIC:
FAIL.** Estes resultados continuam sendo regressão e comparação; não validam
*severe slugging* ISOTHERMAL.

## 10. Regressões antes/depois

| Item | Resultado |
|---|---|
| Suíte Fortran debug (`make test`) no `1d71f00` | **24/24 PASS** (22 anteriores + `test_hydrostatic_outlet` + `test_junction_vanishing_gas`) |
| `make check_all` | compila (avisos preexistentes) |
| Python `unittest discover` | **123 OK** |
| CLI | accepted telemetry, boundary telemetry, feed contract (11 execuções), initial inventory, pipe regime, reduced equilibrium: **6/6 PASS** |
| Byte a byte ISO (30 s, `prims_f_aef3774` × `prims_f_1d71f00`) | monitor, conservation, fields, boundary, checkpoint, rejected, results: **idênticos**; diferem só os tempos de CPU e o eco da chave `JUNCTION_TRACE` |
| Byte a byte BARO (N51/dt 0,125, 30 s, driver `813da75` × `1d71f00`) | todos os arquivos físicos **idênticos**; diferem só os tempos |
| Default `CENTER` vs `f3717dd` | idêntico (relatório anterior, D_OSC/D_Stable 40 s) |
| `tol_res` | `1e-10` em todos os protocolos |

## 11. Conservação

Ledger independente (`analyze_accepted_conservation.audit_ledger`): fluxos aceitos com
pesos de Euler implícito, por fase e por região (sistema, oleoduto e riser). Índice
completo, com hashes de entradas e saídas, em `evidence/RUN_INDEX.json` (80
execuções).

| Conjunto | máx \|balanço\|/M₀ | Aritmética do ledger | Descasamento da integral na junção |
|---|---:|---|---|
| ISO (9 execuções) | ≤ 4,4·10⁻¹⁰ | PASS em todas | 0 exato |
| BARO (4 concluídas) | ≤ 9,4·10⁻⁹ (N51/dt 0,03125, 28 815 passos) | PASS | 0 exato |
| SS1/SS3 (CENTER/STEADY) | ≤ 6·10⁻¹¹ | PASS | 0 exato |

O fluxo de junção é único por fase e entra com sinais opostos nas duas células; não
há fonte de massa na junção.

## 12. Estado e bloqueios remanescentes

| Item | Estado |
|---|---|
| P0.1 hidrostática *well‑balanced* com G≠0 | **PASS** |
| P0.2 certificados SS1/SS3/OSC | **PASS** (classificação); casos **BLOCKED** |
| P1.1 transição misto→puro | **BLOCKED**: limitação documentada; extensão E1/E2 proposta |
| P1.2 SS3 contracorrente inclinado | **BLOCKED**; experimento **NOT_VERIFIED** |
| P1.3 ISO: período | **PASS** (≤ 0,4%) |
| P1.3 ISO: σ | região assintótica (ordem 1,03–1,06 em dt e 0,94–1,01 em h); erro do mais fino 1,2% / 1,5%: **FAIL** estrito no alvo de 1% |
| P1.3 ISO: A*(f) | consistente a ±2,3% (não havia alvo pré‑registrado de 1% para A*) |
| P1.3 ISO: instante do evento | **FAIL** (não convergente; não declarado) |
| P1.3 BAROTROPIC: acomodação v2 | **PASS** (4/4 execuções concluídas) |
| P1.3 BAROTROPIC: convergência de amplitude a 1% | **FAIL** (h 22%; dt 28% → 10%, erro estimado 6,8%) |
| P1.3 BAROTROPIC: refino combinado | **NOT_VERIFIED** (em execução) |
| Validação de *severe slugging* | **NOT_VERIFIED** |

## 13. Próximos experimentos

1. **Fluxo de junção E1** (§6), opcional e com default inalterado byte a byte: avaliar
   SS1/OSC‑ISO e o mapa A*(f). Esse resultado decide se o selo é físico ou numérico.
2. **ISO a partir do estado estacionário discreto** de cada malha, com perturbação
   declarada: testa a convergência do instante do evento.
3. **SS3:** obter os textos e os dados de Barnea et al. (1986) e Ohnuki et al. (1988);
   qualificar no solver local um fechamento estratificado contracorrente assinado
   (Ullmann & Brauner 2006) antes de qualquer integração.
4. **BAROTROPIC:** terminar N101/dt 0,0625 (combinado). Para o eixo h, rodar N201 e
   N401 com dt pequeno; o custo de N201 com dt 0,125 por 900 s é estimado em ~4–5 h.
