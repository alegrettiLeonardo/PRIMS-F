# T1: transmissão conservativa na junção oleoduto–riser

Data: 21/11/2026. Branch `implementation/chexal-riser-20261121`, base `9c92e2f` (confirmada).
`main` em `f3717dd`. Sem push, merge ou PR. Ambiente:
- GNU Fortran 13.3.0;
- BLAS/LAPACK: `-llapack -lblas` resolvem para MKL (`/usr/lib/x86_64-linux-gnu/libmkl_rt.so`);
- uma thread.

Antecedente: [CHEXAL_PHYSICAL_DOMAIN_AND_JUNCTION_QUALIFICATION.md](CHEXAL_PHYSICAL_DOMAIN_AND_JUNCTION_QUALIFICATION.md) §2.

## Resumo

- **Implementação (P0-A):** `junction_transmission='T1'` está implementada como opção experimental, com a T2 como default.
  - Solver T1 de face com camada estacionária **na face**, na mesma estrutura da referência independente.
  - Estacionário inicial próprio da T1.
  - Contrato, checkpoint e saída `junction_transmission.csv` registram a opção.
  - A **T2 continua idêntica byte a byte** nas 16 execuções de regressão; só mudam os metadados acrescentados.
- **Referência (P0-B):** nos mesmos estados conservativos, a face T1 global coincide com a referência T1 local a ≤ 1,7e-4 em G_g* (1,1e-3 no bloqueio do SS1) e ≤ 3e-7 em G_l*, com p_R* − p_L* idêntico ao Pa.
- **Conservação (P0-C):** salto de fluxo zero nas três componentes e identidade de momento do resíduo em ~1e-16. Na T2, a mesma identidade mede exatamente a força pontual.
- **Jacobiano (P0-D):** global ≤ 1,6e-8 nas colunas da junção, inclusive no estado crítico do SS1.
- **Partição homogênea (P0-E):** transparente; diferença de segunda ordem para o Roe comum.
- **Replay e SS1 desde t = 0 (P0-F, P1-B):**
  - A T1 **atravessa o instante em que a T2 exige inversão de gás** (52,94 s) e transmite gás para a frente.
  - O fluxo de gás na junção cai continuamente até zero, na **formação do selo líquido** (~53,5 s).
  - Nesse limite a execução para (NEWTON_GLOBALIZATION_FAILURE, falhas da face 1613/1614/1615): o estado selado (G_g* = 0, α_R* = 0) e a quina do Chexal em j_g = 0 não estão representados.
  - **BLOCKED.** Nenhum ciclo foi completado, então convergência e comparação experimental ficam NOT_VERIFIED.
- **Outros casos com T1:** Stable completa 300 s e decai ao estacionário. OSC bloqueia por contracorrente inclinada no topo (163,75 s, contra 145,70 s com T2). SS3 bloqueia pelo mesmo mecanismo em 41,65 s (T2: 41,55 s), com a junção regular.


## 1. Derivação da condição T1

Burr, *Modelagem unidimensional de sistema oleoduto-riser*, eqs. (72)–(75):

- ∂U/∂t + ∂F/∂s = S;
- U = [(1−α)ρl, αρg, (1−α)ρl u_l + αρg u_g];
- F = [G_l, G_g, F₃], com F₃ = (1−α)ρl u_l² + αρg u_g² + π_P p;
- S = [0, 0, −π_G ρm sen θ − ½ π_F ρm f_m j|j|].

Integrando num volume de controle [s_J − ε, s_J + ε] de área constante, com a junção parada,
os termos de acumulação e as fontes distribuídas são O(ε). A dobra (−2° → 0°) não produz força
tangencial: a reação da parede é normal à bissetriz. A troca de fechamento não é força. No
limite ε → 0:

**G_l,L\* = G_l,R\*, G_g,L\* = G_g,R\*, F₃,L\* = F₃,R\***.

É a condição de Rankine–Hugoniot com velocidade nula. A discretização nodal de Burr,
(128)–(129), já usa um único F_{k,j} por nó. Como F₃ − π_P p = G_l²/(ρl(1−α)) + G_g²/(ρg α)
depende de α, e α_L* ≠ α_R* (fechamentos diferentes), a T1 implica
**p_R\* − p_L\* = [K_L − K_R]/π_P ≠ 0**, com K = G_l u_l + G_g u_g.

## 2. Diferenças físicas T1 × T2

| | T1 (experimental) | T2 (padrão legado) |
|---|---|---|
| Contínuos na face | G_l, G_g, F₃ | G_l, G_g, p |
| Salto de pressão | p_R* − p_L* = ΔK/π_P (+129 a +355 Pa nos estados testados) | 0 |
| Fonte singular de momento | nenhuma | J_M = F₃,R* − F₃,L* (−130 a −355 Pa) |
| Base física | eqs. (72)–(75), balanço integral | não pertence a (72)–(75) |
| Estacionário inicial | EDO com F₃ contínuo na face | EDO com p contínuo |

## 3. Implementação e arquivos modificados

| Arquivo | Mudança |
|---|---|
| `scr/module_config.f90` | `junction_transmission` (default `'T2'`) |
| `scr/main.f90` | namelist, validação (`T1`/`T2`); T1 recusada com E1, E2, BLOCKING ou traços STEADY; `junction_transmission` no `closure_contract.txt`; estacionário T1 |
| `scr/module_simulate_pipe_riser.f90` | `face_pair/basic`: solver T1 (`t1_face`, `t1_res`, `t1_eig`, `t1_field`, `t1_curve2`); cabeçalho do checkpoint; `junction_transmission.csv` (só T1); monitor |
| `scr/module_steady_state.f90` | `t1_junction`: na face da junção a EDO carrega F₃ (`t1_jump`, `t1_gpipe`) |
| `scripts/verify_closure_contract.py` | transmissão do deck = contrato = checkpoint |
| `tests/test_junction_t1.f90` | 14 asserções (§5, §7) |
| `tests/test_closure_contract_cli.py` | caso T1 verificado; transmissão desconhecida recusada |

**Solver T1 da face** (`t1_face`):
- Camada estacionária de comprimento zero **na face**, entre as ondas. Mesma estrutura e
  incógnitas da referência independente (`junction_riemann_ref.f90`):
  - U_L* na 1-onda a partir de U_L (fechamento do oleoduto), parâmetro p_L*;
  - U_M na 3-onda a partir de U_R (riser), parâmetro p_M;
  - U_R* na 2-onda a partir de U_M, parâmetro ln α_R*;
  - condição: F_L(U_L*) = F_R(U_R*).
- **Ondas acústicas** (fracas, Δp/p ~ 1e-3): linearizadas ao longo do autovetor.
- **Onda cinemática** (forte em α): curva integral RK4 com 8 passos em ln α, que mantém α > 0.
- Newton com Jacobiano por diferenças e amortecimento (|Δ ln α| ≤ 3 por iteração).
  Convergência: max|res|/escala < 1e-14, ou < 1e-13 quando a redução por iteração estagna
  (> 0,5), no arredondamento. Com 1e-11 (`7d44f5d`, binário `t1b`) ficava um piso de resíduo
  de ~3e-10 no Newton global, e o SS3 parava em MAX_NEWTON a 10,56 s com status 0. Isso foi
  corrigido em `20e00ba` (binário `t1c`).
- A configuração H4 (λ₁(U_L*) < 0 < λ₂(U_R*) < λ₃(U_M), espectros reais) é verificada na solução.
- Status explícitos 1610–1616. Não há retorno automático para a T2, e α_R* → 0 nunca é projetado.
- Os dois lados da face recebem **o mesmo** vetor de fluxo.
- O Jacobiano da face vem pela rota existente (diferenças do próprio `basic`, centradas ou
  unilaterais admissíveis).

**Implementação descartada, registrada:** a primeira versão mapeava a célula do riser pela
camada T1 **antes** do problema de Roe. Nos mesmos estados, divergia da referência em
+9 % a ×70 em G_g*: o salto de pressão do estado mapeado vira onda acústica espúria. Foi
substituída antes de qualquer campanha. A regressão iniciada com ela foi interrompida e
posta em quarentena (`t1_regress_t2_SUPERSEDED_mapping_impl/`).

## 4. Problema de Riemann local e comparação com o operador global

Referência: `junction_riemann_ref.f90` com o riser Chexal (rodada 3).
- Controles com Chexal: estado uniforme exato; partição contra Roe com erro de segunda ordem;
  ondas pequenas ∝ ε.
- Multistart de 75 partidas: solução T1 única nos estados testados. Isso é busca numérica,
  não prova global.

Hipóteses e limites:
- H1: fontes desligadas no instante t → 0⁺.
- H3: curvas integrais; válidas para rarefações e ondas fracas, enquanto os choques são
  aproximados em segunda ordem.
- Não há referência exata para choques fortes. Nos estados testados, as ondas acústicas têm
  Δp/p ≤ 1,2e-3; a onda cinemática é tratada por curva integral em todos os casos.

Mesmos estados conservativos aceitos, com fixtures do solver Chexal
([T1_FACE_VS_REFERENCE.md](../validation/diagnostics/chexal_riser_20261121/evidence/t1/T1_FACE_VS_REFERENCE.md)):

| Estado | G_g* T1 global | G_g* T1 referência | Dif. rel. G_g* | Dif. rel. G_l* | p_R* − p_L* global / referência (Pa) | G_g* Roe/T2 |
|---|---:|---:|---:|---:|---|---:|
| OSC 60 s | 1,2996 | 1,2996 | 7,5e-6 | 7,3e-8 | 193 / 193 | 1,2848 |
| OSC 135 s | 1,9933 | 1,9933 | 1,5e-5 | 5,5e-8 | 129 / 129 | 1,9853 |
| OSC (último aceito) | 0,61825 | 0,61829 | 7,4e-5 | 2,8e-7 | 267 / 267 | 0,60466 |
| SS3 30 s | 0,85011 | 0,85011 | 3,0e-6 | 1,4e-7 | 167 / 167 | 0,83907 |
| SS3 (último aceito) | 0,61269 | 0,61280 | 1,7e-4 | 0 | 175 / 175 | 0,60081 |
| Stable 600 s | 0,85079 | 0,85080 | 4,1e-6 | 8,4e-8 | 289 / 289 | 0,83078 |
| SS1 30 s | 0,17597 | 0,17596 | 4,2e-5 | 0 | 326 / 326 | 0,16564 |
| **SS1 bloqueio T2 (52,94 s)** | **2,9473e-3** | **2,9440e-3** | 1,1e-3 | 7,6e-8 | 355 / 355 | **3,1e-10** (a T2 só admite G_g* ≤ 0) |

(G_g* em kg/(m²·s).)

A diferença residual vem da linearização das ondas acústicas e dos 8 passos de RK4. Ela é
menor que a diferença Roe/T2 contra a referência T2 (≤ 0,6 %).

**Gate P0-B: PASS** nos estados testados.

## 5. Conservação de massa e momento

Para cada estado, `t1_face_check` monta o resíduo completo e verifica:

Σ w R − [Σ w (U − U_old) − Δt Σ w S + Δt (F_out − F_in)] = Δt Σ_faces (F^m − F^p).

- **T1:** salto de fluxo na face **zero nas três componentes** (bit a bit). Identidade de
  momento ≤ 2,3e-15, arredondamento, em escalas de 1,8–5,8.
- **T2:** a identidade é **exatamente** Δt(F₃,L* − F₃,R*), a fonte singular do esquema;
  as massas fecham a 1e-18.

As massas telescopam nos dois modos. As fontes de gravidade e atrito são distribuídas e
entram em S; não há força concentrada na T1.

Teste unitário (`tests/test_junction_t1.f90`, 14 asserções, PASS):
- junção heterogênea estratificado–Chexal com fluxo idêntico nas três componentes;
- G_g* > 0, p_R* − p_L* = +284 Pa;
- **no mesmo estado, a T2 propõe G_g* < 0 e é recusada (1061)**, como no SS1;
- noutro estado admissível, a T2 mantém o seu salto de momento.

**Gate P0-C: PASS** nos estados testados. Nas execuções, os seis balanços de massa e o contrato estão no §10.

## 6. Jacobianos

- **Face:** blocos J_L, J_R da face T1 contra diferenças centradas independentes do fluxo da
  face, com incremento 30× o interno: erro 1,0e-6 (`test_junction_t1`).
- **Resíduo global montado** (`global_jacobian_check` com `PRIMS_JUNCTION_TRANSMISSION=T1`):
  306 colunas por estado, diferenças centradas em h = 1e-4, 1e-5 e 1e-6 da escala, só em
  estados admissíveis.

| Estado | Erro máx. (melhor h): oleoduto / junção / riser | Colunas > 3e-4 |
|---|---|---:|
| OSC 60 s | 4,7e-8 / 9,0e-9 / 1,5e-9 | 0 |
| SS1 30 s | 3,4e-9 / 8,5e-9 / 3,2e-9 | 0 |
| **SS1 bloqueio T2 (52,94 s)** | 3,4e-9 / 1,6e-8 / 5,2e-9 | 0 |
| Stable 600 s | 2,4e-9 / 5,7e-9 / 1,5e-9 | 0 |

Pontos não diferenciáveis: os status 1613/1615 (curva integral num nó recusado) e 53 (quina do
Chexal) são recusas explícitas, não derivadas fabricadas.

**Gate P0-D: PASS** nos estados testados; ver a limitação no §8.

## 7. Partição homogênea

`tests/test_junction_t1.f90`, com o mesmo fechamento dos dois lados (estratificado) e o
caminho MATCHED ativado (θ_L = −1e-12):

- estado uniforme: F̂(U,U) = F(U) a ≤ 1e-10, p_L* = p_R*, fluxo único (sem salto artificial);
- salto pequeno ε: |T1 − Roe comum|/|F| = 3,4e-6 (ε = 1e-3) e 3,4e-8 (ε = 1e-4), **segunda ordem**.

O teste histórico `test_partition_invariance` (Bendiksen e Chexal) continua PASS.

**Gate P0-E: PASS.**

## 8. Replay diagnóstico desde 52,94 s

Fonte: o estado **aceito** do SS1 Chexal com T2 imediatamente antes da falha (dump do solver,
`last_rejected_state.dat`, colunas U_old; 102 células; SHA-256 do arquivo de U `be2c8932d12b6837…`).
Ele entrou como `initial_u_file`, sem perturbação. Mesmos EOS, fechamentos, malha N51+51,
Δt_max 0,5 s, `tol_res = 1e-10` e 200 s pedidos. A troca T2 → T1 em t = 0 do replay é uma
**mudança instantânea do problema de transmissão**, e o replay é diagnóstico.

| Execução | Binário | Resultado | t alcançado | Passos aceitos/rejeitados | Último status de proposta | CPU (s) |
|---|---|---|---:|---|---|---:|
| T2 | `prims_f_t1` (`26b18c8e…`) | MAX_NEWTON | **0** | 0/16 | 3061 (face: traço horizontal com gás reverso), 2081 vezes | 89 |
| T1 | `prims_f_t1` | COST_STAGNATION (genuína: 31 min sem avanço ≥ 1 s, sem suspensão) | 3,07 s | — | 3614 (face T1 sem convergência) | — |
| T1 | `prims_f_t1b` (`1e1c9fe3…`) | NEWTON_GLOBALIZATION_FAILURE | 3,15 s | 52/28 | 3614; recusas: 940×1613, 887×1615, 318×1614, 140× sem queda de norma | 1404 |

Série da junção (T1, `t1b`): G_g* = 2,95e-3 kg/(m²·s) em 0 s; **3,0e-2 em 0,8 s** (reajuste de
momento após a troca, com o gás passando); 1,7e-2 em 2,5 s; depois cai até o limite em que
a face T1 deixa de convergir. O salto p_R* − p_L* fica em 349–355 Pa.

**Gate P0-F:**
- **PASS** para atravessar o instante da falha T2 com conservação (fluxo único) e admissibilidade
  (G_g* > 0, sem recusa constitutiva).
- **BLOCKED** para a continuação sustentada: ela para no limite de selo (§10).

Ultrapassar o tempo histórico não é validação.

## 9. Inicialização própria da T1

Com `junction_transmission='T1'`, a EDO estacionária carrega F₃, e não p, através da face da
junção (`t1_jump`): kin_L(p_L) + π_P p_L = kin_R(p_R) + π_P p_R, com o fechamento do riser em
θ = 0 e o do oleoduto em θ = −β, resolvido por secante sem clipping. Fechamentos, EOS,
geometria, alimentação, pressão de saída e campos conservativos são os mesmos da T2. O
estado inicial é recalculado por execução, e nenhum IC da T2 é reutilizado.

p_in inicial T1 − T2 (mesmo deck): Stable −290 Pa, OSC −178 Pa, SS1 −332 Pa, SS3 −177 Pa.
É o salto ΔK da junção propagado ao oleoduto, coerente com o §4.

O primeiro passo T1 parte desse estado; a acomodação inicial do esquema FV (resíduos e
balanços) aparece nas execuções do §10, todas com contrato PASS e passos iniciais aceitos.

**Gate P1-A: PASS.** A condição inicial é admissível e consistente com o operador T1, para o
equilíbrio do perfil da EDO. O equilíbrio discreto FV acomoda-se nos primeiros segundos, como
na T2.

## 10. SS1 desde t = 0

Campanha curta discriminante, com IC própria da T1: 300 s, N51+51, Δt_max 0,5 s, BAROTROPIC,
Chexal, `tol_res = 1e-10`. Binário `prims_f_t1b`; o SS3 também com `prims_f_t1c` (ver abaixo).
Todas com contrato PASS.

| Caso (T1) | Resultado | t alcançado | Comparação com a T2 (mesmo deck) |
|---|---|---:|---|
| Stable | COMPLETED; decai ao estacionário | 300 s | estável também com T2 |
| **SS1** | NEWTON_GLOBALIZATION_FAILURE; recusas da face T1 (1080×1615, 874×1613, 389×1614) | **53,48 s** | T2: BLOCKED_JUNCTION_TRACE_DOMAIN (3061) em 52,94 s |
| OSC | BLOCKED_CONSTITUTIVE_DOMAIN (1071, célula 103, 76,08°, j_l = 3,8e-10 m/s) | 163,75 s | T2: o mesmo mecanismo em 145,70 s |
| SS3 `t1b` | FAIL:MAX_NEWTON (piso de resíduo 3,1e-10 do critério 1e-11 da face, sem recusa constitutiva) | 10,56 s | defeito numérico corrigido em `20e00ba` |
| **SS3 `t1c`** | BLOCKED_CONSTITUTIVE_DOMAIN (2617 recusas 1071, célula 86, 69,32°, j_l = 3,3e-9 m/s); junção regular (G_g* = 0,58, α_R* = 0,26) | **41,65 s** | T2: BLOCKED_CONSTITUTIVE_DOMAIN em 41,55 s (mesmo mecanismo) |

**Sequência física do SS1 com T1** (`junction_transmission.csv`):

| t (s) | G_g* na junção (kg/(m²·s)) | α_R* | p_R* − p_L* (Pa) |
|---:|---:|---:|---:|
| 40,45 | 0,179 | 0,098 | 326 |
| 48,45 | 0,0876 | 0,057 | 339 |
| 49,95 | 0,0602 | 0,042 | 344 |
| 51,45 | 0,0329 | 0,026 | 348 |
| 52,95 | 0,0075 | 0,0071 | 354 |
| 53,45 | 0,00035 | 0,00039 | 356 |

O gás deixa de entrar no riser por acúmulo de líquido na base: o selo se forma. Sob T1 isso
ocorre **sem inversão de gás**: G_g* cai continuamente para 0⁺. Sob T2, em 52,94 s a
transmissão exigia G_g* < 0, que o Chexal horizontal recusa.

O SS1 com T1 para exatamente onde o modelo precisaria do estado selado:
- G_g* = 0 com α_R* = 0, isto é, traço do riser só com líquido e gás estagnado no oleoduto;
- o Chexal **não é diferenciável em j_g = 0**: a seleção de Re (Eq. 3-15), C3 ascendente/descendente
  (3-34/3-35) e L (3-7/3-8) mudam de ramo;
- a curva integral cinemática não tem autovetor definido ali (status 1615/1613).

As fases seguintes do SS1 (acumulação, compressão do gás, penetração e blowdown) **não foram
alcançadas**.

**Gate P1-B: BLOCKED** no limite de selo. O próximo passo físico é um tratamento fechado e
fundamentado do estado selado sob T1 (§13), não um piso de α nem uma projeção.

**Balanços.** Os seis balanços de massa (duas fases × sistema, oleoduto e riser) foram calculados
pelos fluxos aceitos (`audit_ledger`, `summarize_campaign.py`), em
[T1_BALANCES.md](../validation/diagnostics/chexal_riser_20261121/evidence/t1/runs/T1_BALANCES.md) e
[T1_RUNS.json](../validation/diagnostics/chexal_riser_20261121/evidence/t1/runs/T1_RUNS.json).

| Execução válida | max \|balanço\|/M₀ (6 balanços) | ledger | contrato |
|---|---:|---|---|
| Stable T1 (300 s) | 4,0e-12 | PASS | PASS |
| SS1 T1 (até 53,48 s) | 4,3e-13 | PASS | PASS |
| OSC T1 (até 163,75 s) | 1,2e-11 | PASS | PASS |
| SS3 T1 `t1b` (até 10,56 s) | 1,9e-12 | PASS | PASS |
| SS3 T1 `t1c` (até 41,65 s) | 9,6e-13 | PASS | PASS |
| Replay T1 `t1b` (até 3,15 s) | 6,8e-13 | PASS | PASS |
| Replay T2 (0 s) | 0 | PASS | PASS |

Os balanços oleoduto e riser usam os fluxos de base à esquerda e à direita. Na T1 eles são o
mesmo vetor, então a soma oleoduto + riser fecha o sistema sem termo de junção.

**Momento.** Nas execuções não há ledger de momento integrado no tempo: o momento não é
conservado globalmente por causa das fontes de gravidade e atrito. A conservação discreta do
momento foi verificada pela identidade do resíduo do §5, que isola dt·Σ(f⁻ − f⁺). Ela vale
~1e-16 para a T1 e é exatamente a força pontual para a T2. Como a T1 entrega um único vetor
de fluxo em todas as faces, a identidade vale em cada passo aceito.

**Veredictos invalidados:** a máquina suspendeu entre 16h02 e 17h46. Três veredictos
COST_STAGNATION medidos com relógio de parede (SS1 `t1`, SS3 `t1b`, replay `t1b`) são artefatos,
marcados com `INVALID_SUSPENSION.txt`. O guarda passou a usar `time.monotonic()` (commit
`0b1e01d`), e as execuções afetadas foram refeitas (`t1b_rerun`).

## 11. Convergência espaço–tempo

Não há trajetória sustentada com ciclos completos sob T1. O SS1 para no selo, OSC e SS3 nos
bloqueios do §10, e o Stable não oscila. Não houve refino de espaço ou de tempo de períodos e
amplitudes, e nenhuma convergência obtida com Bhagwat/T2 foi transportada.

**Gate P1-C: NOT_VERIFIED.**

## 12. Comparação experimental

Baliño, Burr e Nemoto (2010): o PDF não está nesta máquina; há só o registro anterior (SS1
85,6 s, SS3 76 s, OSC 19,3–33,8 s). Além disso, não há ciclo completo, admissível e convergido.

**P2: NOT_VERIFIED.**

## 13. Matriz de status

| Item | Status | Evidência |
|---|---|---|
| P0-A implementação T1 opt-in, T2 default preservada | **PASS** | §3; contrato e checkpoint; 9 casos de contrato no CLI |
| P0-B Riemann local × face T1 global × Roe/T2 | **PASS** | §4; G_g* ≤ 1,7e-4 (1,1e-3 no bloqueio do SS1) |
| P0-C conservação discreta e integral (massa e momento) | **PASS** | §5 (identidade ~1e-16); seis balanços ≤ 1,2e-11·M₀ (§10) |
| P0-D Jacobiano × FD independente | **PASS** | §6; face 1,0e-6, global ≤ 1,6e-8 |
| P0-E partição homogênea / junção heterogênea | **PASS** | §7; `test_junction_t1` 14/14 |
| P0-F replay de 52,94 s | **PASS** atravessar com G_g* > 0 / **BLOCKED** continuação (selo, 3,15 s) | §8 |
| P1-A IC própria da T1 | **PASS** | §9 |
| P1-B SS1 desde t = 0 | **BLOCKED** (selo em 53,48 s; G_g* = 0 / α_R* = 0 e a quina do Chexal em j_g = 0 não representados) | §10 |
| P1-C convergência espaço–tempo | **NOT_VERIFIED** (nenhuma trajetória com ciclos) | §11 |
| P2 Baliño, Burr e Nemoto (2010) | **NOT_VERIFIED** | §12 |
| Stable com T1 | **PASS** (300 s, decai ao estacionário) | §10 |
| OSC com T1 | **BLOCKED_CONSTITUTIVE_DOMAIN** 163,75 s (T2: 145,70 s) | §10 |
| SS3 com T1 | **BLOCKED_CONSTITUTIVE_DOMAIN** 41,65 s (T2: 41,55 s) | §10 |
| Regressão T2 / Bhagwat | **PASS** (16/16 byte a byte em `dc92c80`; réplicas de 30 s SS1-C5 e Stable-C3 idênticas a `20e00ba`, exceto as colunas de tempo de CPU) | `T2_REGRESSION_VS_P3.json` |
| Suítes | **PASS**: Fortran 30/30, Python 123 OK, CLI 8/8 (`closure_contract` 9 casos) | `make test`, `unittest`, `tests/*_cli.py` |
| Contrato de fechamento das execuções | **PASS** em todas | `CLOSURE_CONTRACT_STATUS.json` |

**Conclusão.** A T1 é conservativa e consistente nos estados testados (fluxo único, Jacobiano
correto e partição transparente). Ela elimina o bloqueio da T2 no SS1: não há inversão de gás.
O SS1 com T1 alcança a formação do selo, mas o modelo não representa o estado selado.

**Próximo passo físico:** representar G_g* = 0 / α_R* = 0 na onda cinemática e no Chexal no
limite j_g → 0⁺ com derivação fundamentada, sem piso de α. Também continua aberta a recusa
de contracorrente inclinada (OSC e SS3), independente da transmissão.

Cada falha é limitação do modelo, não validação.

## 14. Commits e comandos

Branch `implementation/chexal-riser-20261121`, a partir de `9c92e2f`. **Sem push nem merge;
`main` intacta em `f3717dd`.**

| Commit | Conteúdo |
|---|---|
| `dc92c80` | feat: `junction_transmission='T1'` experimental (face, estacionário, contrato, checkpoint, teste) |
| `d9b18e1` | diag: `t1_face_check`, transmissão no Jacobiano global e no lançador de campanhas |
| `4f28e12` | fix: transmissão nas verificações de contrato; a chave acrescentada ao namelist preserva o newline final |
| `7d44f5d` | fix: Newton da face até o arredondamento; `junction_transmission.csv` com flush |
| `0b1e01d` | fix(guard): relógio monotônico (suspensão não é estagnação) |
| `2b1c5be` | evidence: face × referência local, identidade de momento, Jacobianos globais, regressão T2 |
| `20e00ba` | fix: critério de estagnação da face 1e-13 (SS3 passa de 10,56 s) |
| `a5cf21e` | evidence: ledger das execuções T1, seis balanços, séries da junção |
| (este commit) | docs: este relatório e o handoff |

**Ambiente congelado:**
- GNU Fortran 13.3.0;
- `-llapack -lblas` → MKL `libmkl_rt.so`;
- OMP 1 thread;
- `tol_res = 1e-10`, N51 + 51, Δt_max 0,5 s, BAROTROPIC, Chexal ar–água (C5) no riser e estratificado (C1) no oleoduto.

| Binário (`…-results-20261121/bin/`) | Fonte | SHA-256 | Usado em |
|---|---|---|---|
| `prims_f_t1` | `dc92c80` | `26b18c8e58da34c4…` | regressão T2 16/16, replay T2 e T1 (`t1_replay`) |
| `prims_f_t1b` | `7d44f5d` | `1e1c9fe310d0f8e4…` | `t1b_short`, `t1b_rerun` |
| `prims_f_t1c` | `20e00ba` (idêntico ao build do commit) | `91a41988594dac41…` | `t1c_ss3` |

Campanhas T1 e T2 separadas por diretório. Nenhuma campanha mistura Chexal e Bhagwat.

**Comandos principais:**
```
make -j4 FC=gfortran LIBS='-llapack -lblas'; make test FC=gfortran LIBS='-llapack -lblas'
PYTHONPATH=scripts:validation/diagnostics/chexal_riser_20261121 python3 -m unittest discover -s tests -p 'test_*.py'
PYTHONPATH=scripts python3 tests/<nome>_cli.py --binary bin/prims_f --output tmp/cli_<nome>
python3 validation/diagnostics/chexal_riser_20261121/run_campaign.py ... --transmission T1 --time 300
validation/diagnostics/chexal_riser_20261121/t1_face_check <fixture> <c_gas> [TRIAL]
python3 validation/diagnostics/chexal_riser_20261121/summarize_campaign.py OUT.json OUT.md RUN_DIR...
```
O caminho de `--output` dos CLI precisa ser curto, porque o Fortran trunca o nome do
arquivo de saída. Os resultados brutos ficam em `~/Documentos/PRIMS-F-chexal-riser-results-20261121/`,
fora do repositório.
