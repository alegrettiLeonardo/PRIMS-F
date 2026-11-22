# PRIMS-F — Qualificação numérica T1 + LIQUID_BRIDGE: operador, Jacobiano implícito e desempenho (P0/P1/P2)

Rodada de 22/11/2026. Branch `implementation/t1-seal-jacobian-20261122`, worktree
`~/Documentos/PRIMS-F-t1seal-jac-20261122`, criada a partir de `56998bb`
(`implementation/chexal-riser-20261121`). Sem push, merge nem MR; `main` intacta (`f3717dd`).
Evidência bruta: `~/Documentos/PRIMS-F-t1seal-jac-results-20261122/` (abaixo `RES/`); cópias
leves em [validation/diagnostics/t1_seal_jac_20261122/evidence/](../validation/diagnostics/t1_seal_jac_20261122/evidence/).

**O que não mudou:**
- física: modelo conservativo de Burr, EOS, geometria, fechamentos estratificado (1) e
  CHEXAL_1992 (5), critérios de recusa constitutiva;
- `tol_res = 1e-10`, tolerância da face 1e-11, Armijo, número máximo de backtracks e
  troca para o Jacobiano completo por diferenças;
- binário de referência `prims_f_seal5_46e8830` (SHA-256 `e9af384c…`).

**O que mudou:** três opções novas de namelist, desligadas por default. Com os defaults, o
binário novo reproduz byte a byte T2, Bhagwat, T1 e a seal legada (§8).

## Resumo executivo

1. **Operadores (P0.1, PASS).**
   - B (H4_ONLY) e C (LIQUID_BRIDGE) coincidem bit a bit no ramo H4.
   - A curva cinemática nova, em (p, U₃) com passos fixos, converge com ordem ≈ 4.
   - A(128) e C(128) concordam a 3e-12 kg m⁻² s⁻¹ em G_g e a 4e-9 Pa em F₃.
   - **Os 273 Pa entre T1 e seal4 vêm inteiramente do histórico de Δt.** Com o mesmo Δt,
     T1 (A8) e C16 diferem ≤ 0,05 Pa até 43 s. Com C16 e Δt = 0,0625 s a diferença é
     35,8 Pa em 20 s, contra 34,7 Pa da seal4 (Δt ≈ 0,065 s).
2. **Jacobiano implícito (P0.2, PASS).**
   - dF*/dW é obtido do sistema convergido H4/H7 pelo teorema da função implícita, com
     solução linear com pivoteamento e sem inversas explícitas.
   - Erro ≤ 7,6e-6 contra Richardson da face completa em 9 estados, 3 deles H7. O FD legado
     erra até 6,6e-5 e custa 20–50× mais.
   - Troca H4/H7: fluxo contínuo e Lipschitz numérico; J₊ ≠ J₋, uma quina genuína.
     **Semissuavidade não demonstrada.**
3. **Newton (P0.3).**
   - **PASS em H4.** Benchmark isolado SS1 0→3 s: 0 rejeições contra 8; 2,5 contra 4,9
     Newton/passo; 6,4 contra 77 faces/passo; **17,3 contra 862 CPU-s por s simulado (50×)**.
   - Com o mesmo operador, FD e implícito dão a mesma solução (2,9e-10 Pa); o implícito é 5× mais barato.
   - **Diagnosticado, não resolvido, em H7.** O piso de resíduo ≈ 1e-10 é de representação:
     a pressão de célula quase líquida obtida de U₁ tem quantum nCl²·ε, e o resíduo o vê
     com fator π_P Δt/Δs. Força Δt ≈ 0,085 s em H7.
4. **Frente e quantidade de movimento (P0.5, PASS).**
   - Rankine–Hugoniot da frente móvel: |s[U]−[F]| ≤ 8,6e-12 relativo em 37 estados H7.
   - A frente não é a face fixa: G_g absoluto no lado do duto é ≠ 0, e só a face tem G_g* = 0.
   - Balanço integral de momento ≤ 5e-15 do bruto (duto, riser, sistema), inclusive durante o selo.
5. **SS1 desde t = 0 (P1.1, PASS parcial).**
   - A T1 original falha em 53,4753 s. A T1 revisada sem selo (B) recusa em 53,472 s (1628).
   - **A T1 revisada com selo e Jacobiano implícito (C) completa:**
     - SEAL em 53,45 s;
     - PEN (penetração do gás) em 67,30 s;
     - BLOW (blowdown) em 84,93 s;
     - novo SEAL em 91,93 s.
     Leva 2073 CPU-s; a seal5 legada não passava de 27,6 s em 6 h.
   - A réplica independente desde 53,4753 s reproduz os eventos a ≤ 0,23 s.
   - Seis balanços de massa ≤ 9e-13·M₀.
6. **Ciclos (P1.2, BLOCKED).**
   - Há um único intervalo SEAL→SEAL: 38,48 s, Δp = 37,7 kPa.
   - Em 92,2 s a queda de líquido no trecho **inclinado** da catenária é contracorrente.
     O CHEXAL_1992 implementado recusa esse estado por domínio da fonte (61 → 1071).
   - Não há segundo ciclo; nenhum critério de recusa foi removido.
7. **Convergência (P1.3, NOT_VERIFIED para período e amplitude).**
   - O erro temporal domina: centenas de Pa a Δt = 0,125 s, com ordem observada 0,6–0,9
     (pré-assintótica).
   - O erro da curva é desprezível. O erro de malha em p_in tem ordem ≈ 1 (0,1–1,2 kPa).
   - **O instante do primeiro selo não converge em malha:** 53,2 s (51 nós), 94,4 s
     (101 nós), e nenhum selo até 92 s (201 nós, provisório).
8. **Baliño 2010 (P2, NOT_VERIFIED).** Não existe ciclo-limite; o intervalo de partida não é
   comparável a T_exp = 85,6 s.

---

## 1. Operadores comparados (P0.1)

### 1.1 Modelo e auditoria textual de Burr

Variáveis conservativas, fluxo e fonte, como no código:

```
U = [(1−α)ρ_l, αρ_g, (1−α)ρ_l u_l + αρ_g u_g]
F = [U₁u_l, U₂u_g, U₁u_l² + U₂u_g² + π_P p],   π_P = ω_P/(ω_ρ ω_u²)
S₃ = −ρ_m (G L_r/ω_u²) sinθ − (2L_r/D) ρ_m f_F j|j|
```

- **Pressão no fluxo.** As eqs. (74) e (184) de Burr incluem +p em F₃; a (135) omite o
  termo. A omissão é tipográfica: sem p o sistema não teria a onda acústica que as próprias
  eqs. (136)–(137) usam. O código segue (74)/(184) (`pip*p`, `module_simulate_pipe_riser.f90`).
- **Atrito.** A eq. (75) traz −½ π_F ρ_m f_m j|j| com f_m de Darcy; a (137) imprime −π_F,
  sem o ½. O código usa `drag = 2 L_r/D` com fator de Fanning f_F = f_D/4, ou seja
  2(L/D)(f_D/4) = ½ π_F f_D. **Coincide com (75).** A (137) é considerada inconsistente
  textualmente; o código não foi alterado.
- **Chexal em α → 0.** u_g = C₀ j + V_gj com C₀ → 0 e V_gj → 0,34419 m/s, e F_g de classe C¹.
  Foi verificado na rodada anterior ([PRIMSF_T1_CHEXAL_SEAL_LIMIT_QUALIFICATION.md](PRIMSF_T1_CHEXAL_SEAL_LIMIT_QUALIFICATION.md) §3)
  e não foi extrapolado a outras condições. Nesta rodada só se usa o limite já verificado.

### 1.2 Operadores

| Operador | Opções | Formulação da face (células 53 \| 54) |
|---|---|---|
| **A** T1 original | `junction_seal='OFF'`, `junction_t1_curve_steps=8` | Newton 3×3 em (p_L*, ln α_R*, p_M). Curva cinemática RK4 em ln α com N passos fixos (8 = original); ondas acústicas por curva integral; fluxo no traço esquerdo |
| **B** revisado sem selo | `junction_seal='H4_ONLY'` | Solver de selo restrito a H4; sem raiz que transmita gás, **recusa com 1628**, sem fallback |
| **C** revisado com selo | `junction_seal='LIQUID_BRIDGE'` | H4 (gás passa) ou H7 (ponte líquida, frente móvel no duto), selecionado pelo sinal de g₀ |
| C, modo de curva | `junction_seal_curve_steps = 0` | Legado seal5: regra adaptativa, curva em (U₁, U₃) com reconstrução da pressão p = P₀ + Cl²(U₁/(1−α) − ρ₀) |
| C, modo de curva | `junction_seal_curve_steps = N` (8–256) | **Novo:** curva em z = √α nas variáveis (p, U₃), com U₁ = (1−α)ρ_l(p) e U₂ = αp/Cg², bem condicionadas. O traço puramente líquido usa a pressão da própria curva: F = (U₃, 0, U₃²/U₁ + π_P p) |
| Jacobiano | `junction_seal_jacobian = 'FD' \| 'IMPLICIT'` | FD legado ou implícito (§2) |

Por que o modo novo existe: no modo legado, o quantum de arredondamento de U₁ vira um quantum
de pressão nCl²·ε ≈ 4,9e-12 (adimensional), porque nCl² = 22 081. Isso produz saltos de F₃ na
face de ±1,67e-12 relativo (§5.1).

### 1.3 Mesmos estados, três operadores

Os estados são congelados de soluções aceitas nas células 53 | 54 (fixtures em `RES/fixtures/`).
A ferramenta é `seal_jac_check OPS`. Referência: C com 128 passos e tolerância da face 1e-12.
**A solução de 8 passos não é usada como referência.**

Diferenças contra a referência, em ΔG_g* (kg m⁻² s⁻¹) / ΔF₃* (Pa) (`RES/p01_ops/P01_TABLE.md`):

| Estado | Config. | ref G_g* | ref F₃* (Pa) | A(8) | A(16) | A(32) | A(128) | C legado | C(8) | C(16) | C(32) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rep_t15 | H4 | 4,793e-2 | 304 519,778 | −1,2e-10 / 1,1e-7 | −7,9e-12 / 7,2e-9 | −5,5e-13 / 7,0e-10 | −2,4e-14 / 1,7e-10 | −8,6e-8 / 7,6e-5 | 6e-15 / 8e-10 | −9e-15 / 2,9e-10 | −3,5e-14 / −1,1e-9 |
| rep_t30 | H4 | 0,8410 | 292 845,004 | ~1e-14 / 6e-11 | idem | idem | idem | −3,8e-11 / −3,8e-7 | ≤1e-14 / ≤3e-10 | idem | idem |
| ss1_dump43.4 | H4 | 0,1843 | 295 021,588 | falha 1614 | 1,7e-14 / 2,9e-10 | idem | falha 1614 | 8,9e-12 / −2,7e-7 | 2e-14 / 2,9e-10 | 3,7e-14 / 5,2e-10 | 4e-14 / 5,2e-10 |
| ss1_t30 | H4 | 0,1612 | 297 478,906 | falha 1614 | 4,4e-14 / −5,8e-11 | idem | −1,1e-13 / 0 | −3,9e-12 / −1,9e-7 | ≤1e-14 / ≤2e-10 | idem | idem |
| ss1_t51.95 | H4 | 2,413e-2 | 296 587,270 | 1,2e-11 / −9,8e-9 | 8,3e-13 / −2,9e-10 | falha 1614 | 1,9e-13 / 2,9e-10 | 5,2e-9 / −4,3e-6 | −3e-15 / 3,5e-10 | 8e-16 / 2,3e-10 | 2,4e-14 / −1,7e-10 |
| **ss1_t53.45** | H4 perto da troca | 3,490e-4 | 297 293,326 | falha 1614 | 1,5e-8 / −1,3e-5 | 8,6e-10 / −7,6e-7 | 2,9e-12 / −3,8e-9 | **3,4e-7 / −3,0e-4** | −1,9e-10 / 1,7e-7 | −1,2e-11 / 1,1e-8 | −1,1e-12 / 2,2e-9 |
| **ss1_T1_dump53.48** | **H7** | 0 | 297 306,264 | **raiz H4 espúria** 3,8e-7 / −0,10 | falha 1614 | falha 1614 | falha 1614 | 0 / −3,0e-4 | 0 / 4,2e-6 | 0 / 4,5e-7 | 0 / 4,7e-8 |

Leitura da tabela:

- **Ordem do RK4 no modo novo.** Perto da troca (53,45 s), o erro de G_g cai de 2,9e-10 (8)
  para 1,2e-11 (16) e 1,1e-12 (32): razões 24 e 11, ordem ≈ 3,5–4,6. Em H7 (F₃), de
  4,2e-6 para 4,5e-7 e 4,7e-8: razão ≈ 9,5. Longe da troca o erro fica no piso da
  tolerância da face (~1e-10 Pa).
- **Normalização.** ε_F por componente: G em kg m⁻² s⁻¹ com escala max(|G_g*|, ρ_g(P_s) j_g0)
  para o gás, e F₃ em Pa. Os erros relativos ficam ≤ 1e-12 fora da vizinhança da troca.
- **A(8)** falha com 1614 em três estados H4 e produz a raiz H4 espúria em 53,48 s, onde
  não há raiz H4 (§3 da rodada anterior). A(16 e acima) falha em 53,48 s com 1614, que é o
  comportamento correto do operador sem extensão: não há raiz.
- **B ≡ C bit a bit no ramo H4** (`test_junction_seal` item 9).
- **Viés do C legado:** até 3,4e-7 kg m⁻² s⁻¹ em G_g, ≈ 1e-3 relativo perto da troca, e
  3e-4 Pa em F₃. O modo novo elimina esse viés.
- **Traços, autovalores e configurações** estão por estado nos logs `OPS2_*`, colunas
  `pL*, pR*`, `aL*, aR*`. As condições de onda (λ₁(U_L1) < 0, λ₂ > s, λ₂,₃(U*) > 0) são
  verificadas pelo próprio solver (status 1621 se violadas): nenhuma violação.

### 1.4 Origem dos 273 Pa (SS1 desde 0, T1 × extensão da rodada anterior)

Diferença p_in(op) − p_in(C16, Δt_max = 0,5) em Pa. As execuções antigas estão em
`~/Documentos/PRIMS-F-chexal-riser-results-20261121/`.

| t (s) | A8 T1 orig. (t1b) | A8 bin. seal0 | C legado seal4 | C16 Δt 0,25 | C16 Δt 0,125 | C16 Δt 0,0625 | C32 Δt 0,5 |
|---|---|---|---|---|---|---|---|
| Δt médio até 27,6 s | 0,327 | 0,316 | 0,065 | 0,203 | 0,114 | 0,059 | 0,327 |
| 1 | 0,0 | −0,2 | −1,8 | 0,0 | −0,4 | −2,1 | 0,0 |
| 10 | 0,0 | −0,1 | −8,6 | −4,9 | −7,5 | −8,9 | 0,0 |
| 20 | 0,0 | 0,2 | 34,7 | 18,6 | 29,6 | 35,8 | 0,0 |
| 30 | 0,0 | 0,5 | 121,2 | 63,9 | 104,2 | 127,8 | 0,0 |
| 43 | 0,0 | −0,7 | −248,0 | −151,8 | −223,0 | −257,0 | 0,0 |
| 50 | 0,0 | −2,0 | — | −433,3 | −667,9 | −787,6 | 0,0 |

Decomposição:

- **(i) Integração da curva:** C16 contra C32 dá 0,0 Pa até o selo e ≤ 2 Pa depois.
- **(ii) Mudança do operador H4:** A8 contra C16, com Δt idêntico porque não há rejeições
  em nenhuma das duas, dá ≤ 0,05 Pa até 43 s e ≤ 2 Pa até 50 s. **As duas H4 são
  fisicamente equivalentes.**
- **(iii) Troca H4/H7:** nenhuma antes de 53,35 s.
- **(iv) Histórico de Δt mais amplificação dinâmica:** explica integralmente os 273 Pa.
  - A seal4 (Δt ≈ 0,065 s) acompanha a C16 com Δt = 0,0625 s: 34,7 contra 35,8 Pa em 20 s.
  - Acompanha também a C16 com Δt = 0,125 s em 30 s: 121 contra 104 Pa.
  - Na rodada anterior a seal4 tinha Δt pequeno por causa das rejeições do Jacobiano por
    diferenças (≈ 1e-4).
  - A diferença de discretização temporal cresce com a dinâmica: 2 Pa em 1 s, 433 Pa em
    50 s (Δt 0,5 contra 0,25).

**Gate P0.1: PASS.** A H4 é fisicamente equivalente nos três operadores, e a fonte da
diferença está quantificada como erro temporal (§11).

---

## 2. Equações do Jacobiano (P0.2)

A face convergida é o sistema implícito

```
E(v, W) = 0,   W = (U_L, U_R) ∈ ℝ⁶,   F* = Φ(v, W)
```

com as incógnitas internas v, as equações E e o mapa de fluxo Φ abaixo:

| Configuração | v | E(v, W) | Φ = F* |
|---|---|---|---|
| H4 | (p_L*, p_M, z), z = √α_R* | [F₁(U_L*) − F₁(U_R*), F₃(U_L*) − F₃(U_R*), F_g(U_L*) − F_g(U_R*)] | F_R(U_R*) |
| H7 | (p_L1, p_M) | [s(U₁* − U₁,L1) − (F₁* − F₁,L1), s(U₃* − U₃,L1) − (F₃* − F₃,L1)] | (U₃*, 0, U₃*²/U₁* + π_P p*) |

Mapas elementares e definições:

- U_L* ou U_L1 = 𝒞₁(U_L; p): onda 1 linearizada ao longo de r₁(U_L).
- U_M = 𝒞₃(U_R; p_M): onda 3.
- U_R* = 𝒦(U_M; z): curva cinemática RK4 com N passos.
- s = F_g(U_L1)/U₂(U_L1) é a velocidade da frente.
- No H7, U* = 𝒦(U_M; 0), o extremo α = 0.

**Derivada.** Sem inversas explícitas:

```
E_v X = −E_W            (eliminação de Gauss com pivoteamento parcial; 1629 se singular)
dF*/dW = Φ_W + Φ_v X
```

**O que é diferenciado numericamente.** Só os mapas elementares explícitos (`t1s_side`):
lado 1 (U_L, p) ↦ (U, F_pipe(U)); lado 2 (U_R, p_M, z) ↦ (U, F_riser(U)). Usa-se diferença
central com passo

```
h_k = ε · min(|U_k|, p/|∂p/∂U_k|),   ε = 1e-5   (junction_seal_ijac_eps)
```

- O limite p/|∂p/∂U₁| é necessário porque dp/dU₁ = c_l²/(1−α) em células líquidas.
- Se uma amostra é inadmissível, tentam-se h/4, 4h e por fim o estêncil unilateral de 2ª
  ordem no lado admissível. A causa típica é o teste de consistência da derivada do Chexal
  (`CLOSURE_DERIVATIVE = 13`) em α ≲ 1e-10.
- Enquanto se diferencia, o número de passos da curva fica congelado (`nz_frozen`). Assim
  a derivada é a do **mesmo mapa RK4 discreto** do integrador, e não a da curva exata.

**O que é exato.** E, Φ, a regra da cadeia e a solução linear. A sensibilidade da curva é a
derivada do RK4 discreto: o integrador é o mesmo e o número de estágios está congelado.

**Status de falha (nunca substituído por outra derivada):** 73 falha de avaliação de mapa
elementar, 74 E_v singular, 75 não calculado.

**Correção `ececf68` (estênceis adicionais).**
- **Sintoma.** As execuções com Δt = 0,125 e 0,0625 s pararam em 69,18 e 68,51 s com
  `ASSEMBLY`, status 2073 = FACE + 73.
- **Causa.** A célula 54 do riser tem α ≈ 2e-14, e o fechamento do riser recusa amostras
  **isoladas** dos mapas elementares (teste de consistência da derivada). A sequência h,
  h/4, 4h e unilaterais pode ser toda atingida. O Jacobiano falhava com ε = 1e-5 e 1e-7 e
  funcionava com 1e-4, 1e-6 e 1e-8 (`evidence/p02_jac/JAC_fail_T95_dt0.0625_fixed.log`).
- **Correção.** Os passos h·{1/16, 16, 1/2, 2, 1/64, 64} são tentados **só depois** de toda
  a sequência original falhar.
  - Toda avaliação que antes tinha sucesso fica inalterada: réplica de 0,3 s byte-idêntica;
    T95_dt0.125 idêntica à jac1 até o passo 674, e o passo 675 é exatamente o primeiro em
    que a jac1 rejeitava.
  - No estado da falha, o implícito com ε = 1e-5 concorda a 5,8e-7 com Richardson.
  - Nenhum fechamento, tolerância ou critério de recusa foi alterado.
- **Campanhas.** As execuções da jac1 que não encontraram o status 73 são válidas para o
  binário corrigido (SS1 C16 até 92,2 s, C32, réplica, Δt 0,25 e 0,5, malha 101).

### 2.1 Verificação contra diferenciação independente

Referência independente: diferenças centrais da **face completa**, solução iterativa inteira
com tolerância da face 1e-12, em dois passos h e h/2 com extrapolação de Richardson. Toda
amostra fica na configuração base: o passo de cada coluna é reduzido (2e-4/4^k) até não
cruzar a troca. Erro máximo da matriz escalada J_ik U_k/F_s,i
(`RES/p02_jac_final/JAC_*.log`, ferramenta final, 16 passos):

| Estado | Config. | implícito ε = 1e-5 | implícito ε = 1e-6 | FD, passos fixos, tol 1e-12 | FD legado seal5 | CPU: implícito / FD fixo / seal5 (s) |
|---|---|---|---|---|---|---|
| rep_t15 | H4 | 1,26e-6 | 3,3e-8 | 4,5e-8 | 2,6e-5 | 0,40 / 4,3 / 10,0 |
| rep_t30 | H4 | 4,0e-7 | 3,4e-8 | 2,0e-8 | 2,2e-6 | 0,45 / 5,1 / 13,4 |
| ss1_dump43.4 | H4 | 6,0e-7 | 2,9e-8 | 4,0e-8 | 1,5e-6 | 0,42 / 4,7 / 12,8 |
| ss1_t30 | H4 | 6,4e-7 | 3,6e-8 | 3,3e-8 | 1,6e-6 | 0,36 / 3,9 / 11,0 |
| ss1_t51.95 | H4 | 1,56e-6 | 5,8e-7 | 4,6e-7 | 4,8e-6 | 0,36 / 5,1 / 14,2 |
| **ss1_t53.45** | H4 perto da troca | **7,6e-6** | 2,8e-4 | 6,8e-6 | 6,6e-5 | 0,92 / 9,9 / 21,6 |
| rep_tol12_dump14.02 | **H7** | 8,1e-7 | 8,3e-8 | 2,7e-7 | 9,7e-8 | 0,20 / 7,2 / 9,9 |
| smoke_rep_ij_dump | **H7** | 8,5e-7 | 1,6e-8 | 1,9e-7 | 1,6e-6 | 0,20 / 3,2 / 6,3 |
| ss1_T1_dump53.48 | **H7** | 8,2e-7 | 2,3e-8 | 1,6e-7 | 2,5e-6 | 0,18 / 2,9 / 5,7 |

- ε = 1e-5 é o único valor admissível em todos os estados. Com ε = 1e-6, o estado de
  53,45 s degrada para 2,8e-4: o passo encontra o ruído da curva perto de z → 0.
- O teste unitário (`test_junction_seal` item 13) repete a verificação no estado de selo e
  num estado H4 independente: 5,8e-6 (H7) e 9,6e-5 (H4), com limiar 2e-4.

**Gate P0.2: PASS**, incluindo o estado crítico junto à troca e três estados H7.

## 3. Derivadas implícitas H4 e H7

- **H7 é derivado separadamente.** Tem duas incógnitas e duas equações de Rankine–Hugoniot,
  com a sensibilidade da velocidade da frente:

  ```
  ds/dW = (∂F_g,L1/∂W)/U₂,L1 − F_g,L1 (∂U₂,L1/∂W)/U₂,L1²
  ```

  Também entra o extremo U* = 𝒦(U_M; 0). **O Jacobiano H4 não é reaproveitado.**
- **Seleção na saída do solver.** Aplica-se o Jacobiano da configuração que gerou o fluxo:
  - H4 na saída da raiz (label 650);
  - H7, ou seja J₋, nas saídas H7 e na cláusula de arredondamento em que g₀ está no nível
    da tolerância. Ali H4(z = 0) = H7(s = 0) é o limite comum, e usa-se a derivada
    unilateral H7.
- **Consistência com o fluxo** (§2.1): H4 em 6 estados ≤ 7,6e-6; H7 em 3 estados ≤ 8,5e-7.

## 4. Troca H4/H7

O parâmetro é a quantidade de movimento do duto U₃,L. A troca U₃c (g₀ = 0) é obtida por
bisseção até 4ε. Os desvios ±δ U₃c vão de 1e-4 a 1e-11 (`RES/p02_switch/`).

| δ | ‖F₊ − F₋‖ (F₁ kg m⁻²s⁻¹, G_g, F₃ Pa) | G_g*(+) | max\|J₊ − J₋\|/max\|J₊\| |
|---|---|---|---|
| 1e-4 | 6,3e-2, 1,6e-4, 26,2 | 1,6e-4 | 0,20 |
| 1e-6 | 6,4e-4, 1,6e-6, 0,26 | 1,6e-6 | 0,20 (53,48) / 0,44 (53,45) |
| 1e-8 | 6,4e-6, 1,6e-8, 2,6e-3 | 1,6e-8 | 0,83 / 0,79 |
| 1e-10 | 6,4e-8, 1,6e-10, 2,6e-5 | 1,6e-10 | 0,85 / 0,82 |

O que está demonstrado (numericamente, nos dois estados):

- **Continuidade do fluxo:** ‖F₊ − F₋‖ → 0 linearmente em δ ao longo de 7 décadas, com
  constante estável (≈ 2,6e5 Pa por unidade de δ em F₃).
- **Lipschitz numérico.**
- **G_g*(+) → 0⁺** do lado H4.

O que não está demonstrado:

- **J₊ ≠ J₋.** Para δ ≥ 1e-6 a diferença relativa é estável em 0,20. A quina é genuína e
  F não é C¹ na troca.
- **Abaixo de δ ≈ 1e-6, J₊ não converge** (0,4–1,0 errático). O Jacobiano H4 em z → 0
  envolve termos α^0,45 do Chexal (F_g é C¹ mas não C² em α = 0), e as diferenças da
  curva perto de z = 0 ficam no ruído.
- Portanto **semissuavidade não demonstrada** e **subdiferencial de Clarke não usado**.
- O Newton global usa o Jacobiano da configuração ativa, uma escolha de elemento do
  B-subdiferencial **não provada**.

**Status: CONTÍNUO E LIPSCHITZ (numérico); C¹ e semissuavidade NOT_VERIFIED.**

## 5. Diagnóstico do Newton (P0.3)

Globalização inalterada:

- função de mérito ‖R‖∞ escalada por [max(U₁,old, 1e-3), max(U₂,old, 1e-7), max(|U₃,old|, 1)];
- Armijo ‖R(x + λd)‖ < ‖R(x)‖(1 − 1e-4 λ), até 18 backtracks;
- Jacobiano completo por diferenças a partir da 7ª iteração;
- nenhum amortecimento ou limite artificial de passo foi introduzido.

### 5.1 Piso de resíduo na face (H4): causa e correção

- **Causa na seal5.** O piso ≈ 1e-10 do resíduo de momento nas células 53/54 vinha da
  reconstrução da pressão na curva cinemática legada, p = P₀ + Cl²(U₁/(1−α) − ρ₀).
  - O quantum de U₁ (ε) vira nCl²·ε ≈ 4,9e-12 em p.
  - F₃ na face salta ±1,67e-12 relativo (`RES/p03_noise/NOISE_*`).
  - A busca linear rejeitava em normas 1,0–2,0e-10: 3850 rejeições de status 0 no SS1 seal5.
- **Correção:** curva em (p, U₃) e fluxo puramente líquido pela pressão da curva (§1.2).
  - Ruído de F₃ na face ≈ 1e-16, com resposta linear a perturbações de 1e-15 a 1e-13
    (`test_junction_seal` item 12: 3,8e-16 em H7, 1,9e-16 em H4).
- **Precisão do Jacobiano:** a fonte do erro ~1e-4 era o FD da face completa, que
  diferencia através do ruído do solver interno (tol 1e-11). O implícito só diferencia
  mapas explícitos suaves (§2), e o erro cai para ≤ 7,6e-6.

### 5.2 Efeito no Newton (benchmark isolado, §6)

| SS1 0 → 3 s | seal5 (C legado + FD) | A (T1, 8 passos) | C16 + FD | **C16 + implícito** |
|---|---|---|---|---|
| passos aceitos / rejeitados | 60 / 8 | 38 / 0 | 36 / 0 | **36 / 0** |
| Newton por passo | 4,90 | 3,29 | 2,50 | **2,50** |
| tentativas de backtrack por Newton | 2,40 | 1,37 | 1,00 | **1,00** (passo pleno aceito sempre) |
| montagens de resíduo | 999 | 296 | 180 | **180** |
| soluções de face por passo | 77,0 | 70,3 | 36,4 | **6,4** |
| CPU por s simulado | 862 | 384 | 87,6 | **17,3** |

- **Mesma solução:**
  - C16 + FD e C16 + implícito diferem 2,9e-10 Pa em p_in. A sequência de iterados e de
    Δt é a mesma, porque os dois Jacobianos têm precisão suficiente para o Newton.
  - seal5 e C16 implícito diferem 5,5 Pa em 3 s. A razão é o histórico de Δt (60 contra
    36 passos, §1.4) mais o viés legado de 1e-7.
- **Rejeições em H4:** 0. As 2 recusas de backtrack pré-selo da trajetória longa estão a
  menos de 0,1 s da troca.

### 5.3 Piso de representação em células líquidas (H7): diagnosticado, não resolvido

Na ponte líquida a trajetória longa tem 1777 rejeições de backtrack, todas de status 0, e
Δt médio de 0,085 s (§9). Diagnóstico offline com Newton denso a partir do estado aceito da
réplica em 9,5708 s (t_abs = 63,05 s), com o Δt = 0,0830 s que falhou. Ferramentas:
`newton_floor.f90` e `ckpt_fixture.py`; log em `RES/floor/newton_floor_9.57.log`.

| Iteração | ‖R‖ | entradas dominantes | max\|dx\|/escala |
|---|---|---|---|
| 0 | 7,2e-2 | R₂, células 58–62 | 6,7e-2 |
| 1 | 2,6e-4 | R₃, células 57–61 | 6,9e-6 |
| 2 | 5,5e-10 | R₃, células 57–62 | 4,3e-12 |
| 3–10 | **1,08e-10 – 3,3e-10, sem decrescer** | **R₃ das células 54–59** (α ≈ 6e-11 – 2e-7) | 3–15e-13 |

- O Newton é quadrático até o piso, e depois os passos ficam no nível de arredondamento.
- **O piso é da célula, não da face.** Em célula quase líquida a pressão sai de U₁ com
  quantum δp = nCl²·ε·U₁ ≈ 4,9e-12. O resíduo de momento o enxerga como

  ```
  δR₃ ≈ (Δt/Δs) · π_P · δp ≈ 0,18 · 361,8 · 4,9e-12 ≈ 3e-10   (um ulp de U₁)
  ```

  Na escala de R₃ (≈ 1), isso fica na ordem de `tol_res`. O piso é proporcional a Δt e,
  nesta malha, cai abaixo de 1e-10 só para Δt ≲ 0,05–0,08 s. Esse é o Δt observado em H7.
- **Não é** erro do Jacobiano, troca H4/H7 na busca linear nem precisão da face.
- **Não foi corrigido nesta rodada.** As correções possíveis mudam incógnitas ou critério:
  - formular as células com a perturbação de densidade δρ_l = ρ_l − ρ₀ ou com p como
    incógnita (mudança das variáveis conservativas);
  - um critério de aceitação que reconheça o piso de representação calculado (mudança do
    critério de convergência).

  Ambas estão fora do mandato ("não alterar tol_res", "não modificar a EOS"). Fica
  registrado como **recomendação**.

### 5.4 Falha em 92,2 s

Depois do novo SEAL (91,93 s), a última fase tem 637 recusas 1071 na célula 94 e 80 recusas
3132/3133 de contorno. O código 1071 é CELL com status 10 + 61, `CHEXAL_DOMAIN`. O líquido
cai no trecho superior **inclinado** da catenária enquanto o gás sobe: o estado é
contracorrente, e a fonte do Chexal só o admite em tubo vertical (`module_chexal_1992.f90`,
`source_coordinates`).

- É uma **recusa constitutiva por desenho**, não removida.
- A réplica falha no mesmo ponto físico: 38,75 s + 53,4753 s = 92,22 s.

**Gate P0.3:**

- **PASS em H4.** Há menos faces, backtracks, rejeições e CPU, com a mesma solução.
- **PARCIAL em H7.** O piso foi diagnosticado e quantificado, mas não removido.

## 6. Profiling e benchmark isolado (P0.4)

**Ambiente** (`RES/bench/environment.txt`):

- CPU i7-1185G7, com governor `powersave` e turbo até 4,8 GHz. Não foi fixado: exige sudo.
  Fica como incerteza de frequência.
- gfortran 13.3.0, `-O2`; MKL `libmkl_rt`.
- `taskset -c 3`, `OMP_NUM_THREADS=MKL_NUM_THREADS=1`.
- **Um processo por vez** (`run_bench.sh`, `run_bench2.sh`).
- Malha 51 + 51, Δt_max 0,5 s, `tol_res` 1e-10, tolerância da face 1e-11, mesmo estado inicial.

**Binários:**

- `prims_f_seal5_46e8830` (`e9af384c…`);
- `prims_f_jac1_191d722` (`bab94fb0…`, ver §14); `prims_f_jac2_ececf68` (`a6f52b16…`) para as reexecuções de §11.

| Execução | Operador | Término | Aceitos/rej. | Newton/passo | backtr./Newton | Jac. | Resíduos | faces/passo | CPU/face (s) | CPU total (s) | **CPU/s sim.** | RSS (kB) | máx. balanço de massa/M₀ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ss1_s5 | seal5 legado | COMPLETED 3 s | 60/8 | 4,90 | 2,40 | 294 | 999 | — | — | 2587 | 862 | 26 764 | 1,1e-12 |
| ss1_L | jac1, defaults legados | COMPLETED 3 s | 60/8 | 4,90 | 2,40 | 294 | 999 | 77,0 | 0,552 | 2581 | 860 | 26 880 | 1,1e-12 |
| ss1_A | A (T1, 8) | COMPLETED 3 s | 38/0 | 3,29 | 1,37 | 125 | 296 | 70,3 | 0,428 | 1151 | 384 | 26 860 | 3,5e-13 |
| ss1_C16fd | C16 + FD | COMPLETED 3 s | 36/0 | 2,50 | 1,00 | 90 | 180 | 36,4 | 0,197 | 263 | 87,6 | 26 752 | 1,5e-13 |
| **ss1_C16ij** | **C16 + implícito** | COMPLETED 3 s | 36/0 | 2,50 | 1,00 | 90 | 180 | **6,4** | 0,207 | **52,0** | **17,3** | 26 864 | 1,8e-13 |
| rep_s5 | seal5 legado | COMPLETED 2 s | 48/4 | 4,33 | 1,71 | 208 | 563 | — | — | 650 | 325 | 26 976 | 8,0e-14 |
| rep_L | jac1, defaults legados | COMPLETED 2 s | 48/4 | 4,33 | 1,71 | 208 | 563 | 65,7 | 0,199 | 645 | 322 | 26 988 | 8,0e-14 |
| rep_C16fd | C16 + FD | COMPLETED 2 s | 33/0 | 2,48 | 1,00 | 82 | 164 | 36,8 | 0,113 | 142 | 70,8 | 26 872 | 3,5e-15 |
| **rep_C16ij** | **C16 + implícito** | COMPLETED 2 s | 33/0 | 2,48 | 1,00 | 82 | 164 | **6,3** | 0,118 | **28,7** | **14,4** | 26 956 | 2,0e-15 |

Diferenças de solução (máx. |Δp_in| na janela):

| Par | Δp_in |
|---|---|
| ss1_s5 : ss1_L | **0** |
| rep_s5 : rep_L | **0** |
| ss1_C16fd : ss1_C16ij | 2,9e-10 Pa |
| rep_C16fd : rep_C16ij | 2,9e-10 Pa |
| ss1_A : ss1_C16ij | 1,2 Pa |
| ss1_s5 : ss1_C16ij | 5,5 Pa |

- O erro de fluxo de cada operador está em §1.3; o de conservação, na última coluna da tabela.
- Os pares s5 : L zerados mostram que o binário novo com defaults reproduz o legado (§8).
- Fonte: `RES/bench/P04_TABLE.md` e `P04_SUMMARY.json`, gerados por `bench_summary.py`.
- O tempo de face inclui as avaliações internas do Jacobiano implícito. Por face resolvida,
  o custo é igual nos modos C (≈ 0,2 s): o ganho vem do número de faces por passo, 6 contra 36.

**COST_STAGNATION (auditoria).**

- O guarda (`run_guarded.py`) agora lê o tempo simulado do último passo **aceito**
  (`steps.csv`); antes lia `monitor.csv`, que pode atrasar até um intervalo de saída.
- Usa relógio monotônico.
- Grava em `GUARD_STATUS.json` a comparação exata que disparou: janela, tempos e
  `advance_s < min_advance_s`.
- Reauditoria da seal4: avançou 1,29 s em 1833 s de parede, uma parada limítrofe,
  agravada pela granularidade do monitor.

**Gate P0.4: PASS.**

## 7. Frente de Rankine–Hugoniot e balanço de quantidade de movimento (P0.5)

### 7.1 Frente móvel × face fixa

Na configuração H7 a face fixa (53 | 54) transporta F* = (U₃*, 0, F₃*), com G_g* = 0 exato.
A frente U_L1 → U* é uma descontinuidade **móvel** dentro do duto, com velocidade s ≤ 0.

Verificação independente (`seal_jac_check RH`, `RES/p05_rh/`):

- o solver registra (s, U_L1, U*, p_L1, p*) em `junction_seal_front`, só escrita e nunca
  lida pela numérica;
- F(U_L1) é recalculado por `point_state` com o fechamento estratificado;
- F(U*) é o fluxo líquido puro.

| Grandeza (37 estados H7: 8 fixtures; U₃,L × 1; 0,99; 0,97; 0,95; 0,9; 0,8) | Resultado |
|---|---|
| s | −8,1e-5 a −0,22 m/s |
| \|s[U₁] − [F₁]\| / escala | ≤ 8,6e-12 |
| \|s[U₂] − [F₂]\| / escala | 0 (por construção de s) |
| \|s[U₃] − [F₃]\| / escala | ≤ 7,8e-13 |
| fluxo de gás através da frente, F_g,L1 − s U₂,L1 | 0 |
| G_g **absoluto** do lado do duto, F_g,L1 | −1,1e-4 a −0,30 kg m⁻² s⁻¹ (≠ 0, para trás com a frente) |
| G_g* na face fixa | 0 |
| F*(pressão da curva) − F*(`point_state`) | ≤ 1e-12 relativo (o quantum c_l² de §5) |

- **G_g = 0 não é imposto na frente.** Só a face fixa tem G_g* = 0. Na frente vale
  [F] = s[U], e o gás é transportado com ela.
- `test_junction_seal` item 14 repete a verificação no estado de selo do SS1: 8,6e-15 e
  3,6e-15 relativo; F_g,L1 = −1,3e-4 kg m⁻² s⁻¹. Foram 3 asserções novas; o total é 38.

### 7.2 Balanço integral de momento por segmento

Ledger por passo aceito (`momentum_ledger=.true.`, `momentum_ledger.csv`), com todos os
termos avaliados no estado aceito t^{n+1} (Euler implícito):

```
duto : I_p^{n+1} − I_p^n = Δt(F₃,in − F₃,jL + S_grav,p + S_fric,p) + Δt R_p
riser: I_r^{n+1} − I_r^n = Δt(F₃,jR − F₃,out + S_grav,r + S_fric,r) + Δt R_r
```

- I = Σ Δs U₃ é o armazenamento.
- R = Σ Δs R₃/Δt é o resíduo do iterado aceito.
- A junção T1 tem um único fluxo: F₃,jL = F₃,jR.

Normalização: massa bruta transportada e produzida, Σ Δt(|F₃,in| + |F₃,out| + |S_grav| + |S_fric|).
Ferramenta: `momentum_balance.py`.

| Execução | Janela | E_P/bruto duto | riser | sistema | máx. por passo | (E_P − ∫R)/bruto |
|---|---|---|---|---|---|---|
| ss1_C16ij | 0 → 92,2 s (H4, H7, PEN, BLOW, SEAL) | 4,8e-15 | −3,5e-15 | 1,4e-15 | 2,4e-10 | ≤ 1,5e-17 |
| rep_C16ij | 53,48 → 92,2 s | 3,6e-15 | −2,8e-15 | 7,8e-16 | 1,7e-11 | ≤ 2,7e-17 |
| ss1_C32ij | 0 → 60 s | 2,6e-15 | −6,4e-16 | 2,0e-15 | 1,2e-12 | ≤ 7,5e-18 |
| ss1_C16ij Δt 0,125 | 0 → 60 s | 2,3e-15 | −3,2e-16 | 2,0e-15 | 1,2e-12 | ≤ 6,5e-18 |

- Salto na junção: Σ Δt(F₃,jL − F₃,jR) = 0 exatamente.
- Termos integrados ∫(·)dt no SS1 0 → 92,2 s, em Pa·s = kg m⁻¹ s⁻¹: F₃,in = 2,67e7;
  F₃,out = 1,92e7; gravidade no duto +7,98e5 (declive de 2° para baixo) e no riser −8,02e6;
  atrito no duto −2,23e5 e no riser −9,45e4.
- Variação do armazenamento: duto +70,2 e riser −5562 kg m⁻¹ s⁻¹.
- O fechamento é a identidade discreta, até o resíduo do Newton (≤ tol_res).
- O máximo por passo (2,4e-10) ocorre na fase final de Δt ~ 1e-6, onde o bruto do passo é minúsculo.

**Gate P0.5: PASS.** Balanços fechados; nenhuma lacuna a classificar.

## 8. Regressões

| Verificação | Resultado |
|---|---|
| jac1 (defaults) × referência seal5: T2 SS1 30 s, T2 Stable (Bhagwat) 30 s, T1 SS1 5 s, réplica seal legada 1 s | **byte-idênticos**: `monitor`, `conservation`, `fields`, `boundary_conditions`, `rejected_proposals`, `junction_transmission`, `checkpoint`, colunas 1–12 de `steps`, contrato (exceto `input_file`) |
| benchmark ss1_L × ss1_s5, rep_L × rep_s5 (jac1 com defaults × seal5) | Δp_in = 0; contadores idênticos |
| release final (com `junction_seal_front`) × jac1 congelado, réplica implícita 0,3 s | **byte-idênticos** em 8 arquivos e em `steps` 1–12 (`RES/regr2/`) |
| `ececf68` (estênceis adicionais) × jac1, réplica implícita 0,3 s | **byte-idênticos** em 8 arquivos; `make test` 31/31 |
| Fortran (`make test`) | **31/31 PASS**; `test_junction_seal` com 38 asserções |
| Python (`PYTHONPATH=scripts python3 -m unittest discover -s tests`) | **123 OK** |
| CLI (`tests/*_cli.py`) | **8/8 PASS** (após limpar diretórios de saída preexistentes; a primeira tentativa falhou só por `FileExistsError`) |
| `make check_all` (`-O0 -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan`) | build OK; réplica com selo, Jacobiano implícito e ledger, 0,05 s: COMPLETED, 14 passos, sem armadilha (`RES/debug/`) |

## 9. SS1 desde t = 0 (P1.1)

Todas as execuções partem do mesmo estado inicial D_SS1 (Chexal ID 5, barotrópico, 51 + 51,
Δt_max 0,5 s). O horizonte curto (3 s, §6) foi feito antes do longo.

| Operador | Binário | Execução | Término | Alcance | Passos aceitos/rej. | CPU (s) |
|---|---|---|---|---|---|---|
| A, T1 original | `prims_f_t1b` (7d44f5d) | `seal_p01_t1b` (rodada anterior) | NEWTON_GLOBALIZATION_FAILURE | **53,4753 s** | 152/28 | — |
| A, T1 com 32 passos | jac1 | `ss1_A32` | interrompida pelo operador | 7e-3 s | 5/0 | 114 (≈ 65 CPU-s por passo de 2e-3 s; > 20 h projetadas) |
| **B**, revisada sem selo | jac1 | `ss1_B16ij` | NEWTON_GLOBALIZATION_FAILURE, recusa 1628 na face 53 | **53,472 s** | 146/21 | 1443 |
| C legado | seal5 | `seal5_ss1` (rodada anterior) | TIMEOUT 6 h | 27,61 s | 417/94 | — |
| **C16 + implícito** | jac1 | `ss1_C16ij` | NEWTON_GLOBALIZATION_FAILURE, Chexal 1071, §5.4 | **92,198 s** | 389/72 | **2073** |
| C32 + implícito | jac1 | `ss1_C32ij` | COMPLETED 60 s | 60 s | 170/8 | 1153 |
| C16 + implícito, réplica desde 53,4753 s | jac1 | `rep_C16ij` | NEWTON_GLOBALIZATION_FAILURE, Chexal 1071 | 38,75 s (92,22 s abs.) | 296/79 | 1708 |

- A réplica desde 53,4753 s é só verificação cruzada; não substitui a campanha desde t = 0.
- **B ≡ C até o selo:** B16ij e C16ij são idênticos em `steps.csv` (colunas 1–12) e em
  `monitor.csv` nos primeiros 135 passos, até 52,948 s.
  - No passo 136, C aceita Δt = 0,5 s terminando em H7 (53,448 s).
  - B recusa (1628), como especificado (sem raiz H4 e sem H7), reduz Δt e falha em 53,472 s.
- **C atravessa o selo** com só 2 falhas de avaliação da face (3625 = FACE, 2000 + 1625) em
  169 passos de H7.

Eventos e custo por fase (ss1_C16ij):

| Fase | Intervalo (s) | Passos | Newton/passo | Δt médio (s) | CPU (s) | CPU/s sim. | Recusas de backtrack |
|---|---|---|---|---|---|---|---|
| H4 pré-selo | 0 – 53,35 | 135 | 2,85 | 0,392 | 393 | 7,4 | 2 (status 0) |
| ponte H7 | 53,35 – 67,30 | 169 | 3,25 | 0,085 | 1014 | 70,6 | 1777 (status 0, piso §5.3) + 2 × 3625 |
| penetração e produção (H4) | 67,30 – 84,93 | 62 | 3,52 | 0,284 | 255 | 14,4 | 216 (status 0) |
| blowdown | 84,93 – 91,93 | 14 | 4,14 | 0,500 | 61 | 8,7 | 15 + 117 × 1071 |
| novo selo e queda de líquido | 91,93 – 92,20 | 9 | 4,78 | 0,030 | 351 | 1291 | 1743 (status 0) + 637 × 1071 + 86 × 3132/3133 |

Concordância entre SS1 desde 0 e a réplica (tempos absolutos):

| Evento | SS1 desde 0 | Réplica | Diferença |
|---|---|---|---|
| PEN | 67,30 s | 67,53 s | 0,23 s |
| BLOW | 84,93 s | 84,77 s | 0,16 s |
| novo SEAL | 91,93 s | 91,77 s | 0,16 s |
| falha | 92,20 s | 92,22 s | 0,02 s |

Os históricos de Δt diferem: a réplica parte de um estado T1 de 8 passos.

**Gate P1.1: PASS parcial.**

- A trajetória desde t = 0 passa pela formação do selo, bloqueio, penetração, blowdown e
  novo selo.
- Massa e momento se conservam a ≤ 1e-12 e ≤ 5e-15.
- Não chega a um segundo ciclo, pela restrição de domínio do fechamento (§5.4).

## 10. Blowdown e ciclos (P1.2)

`cycle_metrics.py` usa eventos pré-declarados:

- **SEAL:** G_g na base = 0 depois de amostra > 0;
- **PEN:** G_g > 0 depois de SEAL;
- **BLOW:** α_out cruza 0,05 depois de PEN.

O monitor é amostrado a 0,1 s.

| ss1_C16ij | Valor |
|---|---|
| eventos | SEAL 53,45 → PEN 67,30 → BLOW 84,93 → SEAL 91,93 s |
| intervalo SEAL→SEAL T₁ | 38,48 s |
| duração da ponte | 13,85 s |
| SEAL → BLOW | 31,48 s |
| p̄_in no intervalo | 289,35 kPa |
| p_max (em 71,93 s) / p_min | 297,84 / 260,14 kPa |
| Δp = p_max − p_min (A_p = Δp/2) | 37,70 kPa (18,85 kPa) |
| inventário no SEAL₁ | M_l 79,64 kg, M_g 0,178 kg |

- **Há um único intervalo, a partir da partida.** Não se declara periodicidade, e a
  dispersão e a deriva de inventário entre ciclos não são calculáveis.
- **Status: BLOCKED** pelo domínio do fechamento CHEXAL_1992 em contracorrente inclinada,
  em 92,2 s. Este é o próximo bloqueio físico. Removê-lo exige estender o fechamento, o que
  está fora do mandato desta rodada.

## 11. Convergência (P1.3)

**Estado em 22/11/2026, 20:03 (-03):** três execuções com o binário corrigido `ececf68`
continuam rodando sob o guarda (máx. 6 h). Os números marcados *provisório* vêm de prefixos
dessas execuções.

### 11.1 Tempo

Δt_max = 0,5; 0,25; 0,125; 0,0625 s, mesma malha (51 + 51). Em H4 o Δt fica em Δt_max
(sem rejeições), então o refinamento é limpo até o selo. Depois do selo o Δt é adaptativo
(piso de §5.3), e o refinamento deixa de ser controlado.

| t (s) | p_in Δt 0,5 (Pa) | d(0,5 − 0,25) | d(0,25 − 0,125) | d(0,125 − 0,0625) | p_t (0,5; 0,25; 0,125) | p_t (0,25; 0,125; 0,0625) |
|---|---|---|---|---|---|---|
| 10 | 290 387 | 4,9 | 2,7 | 1,3 | 0,88 | 1,01 |
| 20 | 290 730 | −18,6 | −11,0 | −6,2 | 0,76 | 0,81 |
| 30 | 290 379 | −63,9 | −40,3 | −23,6 | 0,67 | 0,77 |
| 45 | 288 173 | 242,4 | 126,1 | 65,9 | 0,94 | 0,94 |
| 50 | 288 758 | 433,3 | 234,7 | 119,7 | 0,88 | 0,97 |
| 53 | 289 889 | 468,5 | 277,1 | 152,0 | 0,76 | 0,87 |
| 60 | 293 821 | 222,6 | 147,0 | 80,1 | 0,60 | 0,88 |

- **A ordem temporal observada tende a 1** (Euler implícito): 0,8–1,0 na terna mais fina.
- **Erro estimado (Richardson, p = 1) da solução com Δt = 0,0625 s:** ≈ 120 Pa em 50 s e
  ≈ 150 Pa em 53 s.
- **Eventos:**

  | Evento | Δt 0,5 | Δt 0,25 | Δt 0,125 | Δt 0,0625 |
  |---|---|---|---|---|
  | SEAL | 53,45 s | 53,42 s | 53,18 s | 53,33 s |
  | PEN | 67,30 s | 69,3 s | não alcançado (69,9 s, provisório) | não alcançado (61,9 s, provisório) |
  | BLOW | 84,93 s | 86,7 s | — | — |
  | novo SEAL | 91,93 s | 93,2 s | — | — |

  - SEAL varia de modo não monotônico, com amostragem de 0,1 s.
  - Intervalo SEAL→SEAL: 38,48 s (Δt 0,5) contra 39,78 s (Δt 0,25), **diferença de 3,4 %**.
    Δp: 37,7 contra 40,8 kPa (8 %).
- A execução com Δt 0,25 s também termina por CHEXAL_DOMAIN (1071, célula 92) em 93,3 s.

### 11.2 Malha

Malhas de 51, 101 e 201 nós por segmento, Δt_max = 0,125 s. Até o selo o histórico de Δt é
idêntico nas três (Δt = Δt_max).

| t (s) | p_in 51 | d(51 − 101) | d(101 − 201) | p_s observado |
|---|---|---|---|---|
| 5 | 290 297 | 93,8 | 43,0 | 1,13 |
| 10 | 290 380 | 140,4 | 63,6 | 1,14 |
| 20 | 290 759 | 321,3 | 155,6 | 1,05 |
| 30 | 290 483 | 131,0 | 85,2 | 0,62 |
| 40 | 288 708 | −794,8 | −348,8 | 1,19 |
| 45 | 287 805 | −1166,6 | −573,1 | 1,03 |
| 50 | 288 090 | −843,0 | −519,2 | 0,70 |

- **Ordem espacial ≈ 1 em p_in** até 50 s, como esperado para Roe de 1ª ordem
  (`use_muscl=.false.`). Erro estimado da malha de 201 nós ≈ 0,3–0,6 kPa em 40–50 s.
- **O primeiro selo não converge em malha:**
  - 51 nós: SEAL em 53,2 s;
  - 101 nós: SEAL só em **94,4 s**, e falha por CHEXAL_DOMAIN (1071, célula 187) em 94,9 s;
  - 201 nós: **nenhum selo até 92,06 s** (`N201_dt0.125_T130`, provisório, em execução até 130 s).
- O mínimo da vazão de gás na base até 60 s é 0 (51 nós), 1,8e-5 (101) e 1,4e-4 kg/s (201).
- O inventário de líquido diverge 1,1 kg entre as malhas em 60 s. Na malha grossa, a difusão
  numérica antecipa o acúmulo na base até o limiar de selagem.

### 11.3 Conclusão de P1.3

- Δt: ordem observada → 1 em p_in; erro ainda O(100 Pa) com Δt = 0,0625 s.
- Malha: ordem ≈ 1 em p_in; **o instante de selagem desde t = 0 não está convergido**
  (diferença de 41 s entre 51 e 101 nós).
- O período e a amplitude a 1 % **não podem ser verificados**: há um só intervalo de
  partida, que difere 3,4 % entre Δt 0,5 e 0,25, e a sequência de eventos muda com a malha.
- **Status:**
  - **NOT_VERIFIED** para período e amplitude;
  - **FAIL (não convergido)** para o instante do primeiro selo em malha;
  - **PASS com ordem ≈ 1** para p_in no regime H4 (0–50 s).

## 12. Comparação experimental (P2, Baliño, Burr & Nemoto 2010)

Valores do registro `B_SS1` em [docs/benchmark_sources.json](benchmark_sources.json). Eles vêm da
Tabela 3, caso 2, p. 653, e foram verificados visualmente no PDF fornecido pelo usuário em
18/11/2026. **O PDF não estava disponível nesta sessão;** a nova busca web também não encontrou
uma cópia livre. Por isso a verificação nesta rodada é por esse registro.

| Grandeza | Experimental (Wordsworth et al. 1998, via BBN2010) | Simulação publicada (BBN2010) | Esta rodada |
|---|---|---|---|
| P_s (bar abs.) | 2,08 ± 0,03 | — | 2,08 (deck D_SS1) |
| j_g0 padrão (m/s) | 0,11 | — | 0,11 |
| j_l0 (m/s) | 0,528 | — | 0,53 (D, Tab. 6.2; diferença mantida) |
| T (s) | **85,6** | 80,4 (informado no pedido; **não consta do registro**, não verificado) | sem ciclo-limite; intervalo de partida 38,48 s |
| ΔP pico a pico (bar) | **0,51** | — | 0,377 (só no intervalo de partida) |

- 80,4 s **não** foi usado como valor experimental. A meta experimental é 85,6 s / 0,51 bar.
- O intervalo de partida é um transiente que começa em t = 0, não um ciclo-limite; compará-lo
  com T_exp seria indevido. Nenhum coeficiente foi ajustado.
- **Status: NOT_VERIFIED.** A comparação exige ciclos estabelecidos, bloqueados por §10.

## 13. Matriz de status

| Item | Status | Evidência |
|---|---|---|
| P0.1 Três operadores, mesmos estados; B ≡ C em H4; convergência do RK4 | **PASS** | §1.3; `test_junction_seal` 9–10 |
| P0.1 Origem dos 273 Pa | **PASS**: histórico de Δt + amplificação; operador ≤ 0,05 Pa | §1.4 |
| P0.2 Jacobiano implícito (H4 e H7, derivação separada) | **PASS** ≤ 7,6e-6 contra Richardson; estado α ≈ 2e-14 corrigido em `ececf68` (5,8e-7) | §2–3; item 13 |
| P0.2 Troca: continuidade e Lipschitz | **PASS (numérico)** | §4 |
| P0.2 Troca: C¹, semissuavidade e Clarke | **NOT_VERIFIED**: J₊ ≠ J₋, J₊ não converge para δ < 1e-6 | §4 |
| P0.3 Newton em H4 (menos faces, backtracks, rejeições e CPU; mesma solução) | **PASS**: 50× em CPU/s sim. contra seal5 | §5.2, §6 |
| P0.3 Piso de resíduo em H7 | **FAIL (diagnosticado)**: piso de representação ∝ Δt; correção fora do mandato | §5.3 |
| P0.4 Benchmark isolado e auditoria do COST_STAGNATION | **PASS** (governor não fixado: registrado) | §6 |
| P0.5 Rankine–Hugoniot da frente móvel | **PASS** ≤ 8,6e-12 | §7.1; item 14 |
| P0.5 Balanço de momento por duto, riser e sistema | **PASS** ≤ 5e-15 | §7.2 |
| Regressões T2 / Bhagwat / T1 / seal legada; suítes; check_all | **PASS** | §8 |
| P1.1 SS1 desde t = 0, três operadores | **PASS parcial**: A falha em 53,4753 s; B recusa em 53,472 s; C chega a 92,2 s com selo, PEN, BLOW e novo SEAL | §9 |
| P1.2 Ciclos completos, T̄, Δp, dispersão, deriva | **BLOCKED**: contracorrente inclinada fora do domínio do CHEXAL_1992 (1071) em 92,2 s | §5.4, §10 |
| P1.3 Convergência Δt / malha (p_in, H4 0–50 s) | **PASS**: p_t → 0,8–1,0, p_s ≈ 1; erro O(100 Pa) | §11.1–11.2 |
| P1.3 Instante do primeiro selo em malha | **FAIL (não convergido)**: 53,2 / 94,4 / > 92 s (provisório) | §11.2 |
| P1.3 Período e amplitude a 1 % | **NOT_VERIFIED** (um só intervalo) | §11 |
| P2 Baliño 2010 | **NOT_VERIFIED** | §12 |

## 14. Commits, hashes e comandos

Branch `implementation/t1-seal-jacobian-20261122` (base `56998bb`):

| Commit | Conteúdo |
|---|---|
| `191d722` | feat(seal): Jacobiano implícito da face, curva cinemática de passos fixos bem condicionada, operador B (H4_ONLY), ledger de momento, contadores de face, guarda auditável |
| `0c392c6` | test(seal): registro da frente H7 e teste de Rankine–Hugoniot, ferramentas offline (Newton denso, fixture de checkpoint), evidência P0/P1 |
| `ececf68` | fix(seal): estênceis adicionais do Jacobiano implícito só depois de todos falharem |
| HEAD da branch publicada | docs: este relatório, handoff e evidência de §11 |

Binários, em `~/Documentos/PRIMS-F-t1seal-jac-results-20261122/bin/`, não versionados:

| Binário | SHA-256 |
|---|---|
| `prims_f_jac1_191d722` | `bab94fb08f82e4897fd58f8cc988f1c91ba9bb62ecbea313b0affac41b0b91c0` |
| `prims_f_jac2_ececf68` | `a6f52b16a98cd8495777d2fb71964fc8b149f2ec15144570366b6accbd4738ef` |
| referência `prims_f_seal5_46e8830` | `e9af384c98aebc337e1a3f1599b4356f9dae5f0c4e7659c2f45a4333542a8601` |

O patch contra `56998bb` está em `RES/bin/source_diff_jac1_191d722_vs_56998bb.patch`.

**Comandos** (no worktree; `RES = ~/Documentos/PRIMS-F-t1seal-jac-results-20261122`):

```bash
make release && make check_all && make test                       # Fortran 31/31
PYTHONPATH=scripts python3 -m unittest discover -s tests -p 'test_*.py'   # 123
for t in tests/*_cli.py; do b=$(basename $t .py); rm -rf tmp/fv_$b; python3 $t --binary bin/prims_f --output tmp/fv_$b; done   # 8/8
# ferramentas offline (face, Jacobiano, troca, ruído, RH, piso)
D=validation/diagnostics/t1_seal_jac_20261122
$D/build.sh $D/seal_jac_check.f90 $RES/tools/seal_jac_check
$D/build.sh $D/newton_floor.f90 $RES/tools/newton_floor
export PRIMS_PROBE_CLOSURE_RISER=5 PRIMS_JUNCTION_TRANSMISSION=T1 PRIMS_JUNCTION_SEAL=LIQUID_BRIDGE PRIMS_SEAL_STEPS=16
$RES/tools/seal_jac_check {OPS|JAC|SWITCH|NOISE|RH} $RES/fixtures/<estado>
python3 $D/ckpt_fixture.py <checkpoint.dat> $RES/fixtures/smoke_rep_ij_dump <saida>; $RES/tools/newton_floor <saida>
# benchmark isolado (um processo por vez, núcleo 3)
bash $RES/bench/run_bench.sh; bash $RES/bench/run_bench2.sh
python3 $D/bench_summary.py $RES/bench/runs P04.json P04.md ss1_C16fd:ss1_C16ij ss1_s5:ss1_L rep_C16fd:rep_C16ij
# campanhas (decks em $RES/camp/decks; cópias em evidence/p1_campaign)
python3 validation/diagnostics/junction_constitutive/run_guarded.py <binario> <deck.nml> <run_dir> --max-wall-s 21600 --window-s 1800 --min-advance-s 1
python3 $D/cycle_metrics.py <run_dir> ; python3 $D/momentum_balance.py <run_dir>
```

Opções de namelist novas (defaults entre parênteses):

- `junction_seal`: OFF | LIQUID_BRIDGE | H4_ONLY (OFF);
- `junction_seal_curve_steps`: 0 | 8–256 (0);
- `junction_seal_jacobian`: FD | IMPLICIT (FD);
- `junction_t1_curve_steps`: 4–256 (8);
- `momentum_ledger` (.false.).

O deck das campanhas é D_SS1 com `junction_transmission='T1'`, `junction_seal='LIQUID_BRIDGE'`,
`junction_seal_curve_steps=16`, `junction_seal_jacobian='IMPLICIT'` e `momentum_ledger=.true.`.
