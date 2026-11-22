# PRIMS-F — T1 + Chexal: limite de selo líquido na junção

Data: 21–22/11/2026 · Branch `implementation/chexal-riser-20261121` (a partir de `538df1e`) · `main` = `f3717dd`, intacta · sem push nem merge.

Este documento qualifica a falha do SS1 com T1 em 53,48 s e a extensão experimental
`junction_seal='LIQUID_BRIDGE'`. Antecedentes: [PRIMSF_T1_JUNCTION_QUALIFICATION.md](PRIMSF_T1_JUNCTION_QUALIFICATION.md).
Evidências: [validation/diagnostics/t1_seal_20261121/](../validation/diagnostics/t1_seal_20261121/).

## Resumo

- **Causa da falha em 53,48 s (P0.1, PASS): ausência local da fase gasosa na face.**
  - Não é degeneração constitutiva do Chexal nem falha primária do método não linear.
  - O mapa local independente mostra que, a partir de 53,45–53,475 s, **não existe raiz da configuração
    de transmissão de gás (H4)** para nenhum α_R* > 0: g₀ troca de sinal de +0,020 para −6,6e-3 kg/(m²·s).
  - A "raiz" aceita pela T1 em 53,4753 s era artefato de truncamento da curva cinemática de 8 passos.
  - A falha do Newton (1613/1614/1615) é **consequência**: Δt até 3e-8 s não muda nada.
- **Limite Chexal α → 0 (P0.2, PASS):** C₀ → 0, V_gj → 0,34419 m/s, fluxo de gás C¹, espectro
  estritamente hiperbólico. A fase ausente exige mudança de configuração, não piso de α.
- **Extensão `junction_seal='LIQUID_BRIDGE'` (P1.1–P1.3, PASS, experimental):**
  - **H7, ponte líquida:** G_g* = 0 exato, com frente de acumulação no oleoduto (Rankine–Hugoniot,
    s = u_g ≤ 0).
  - **Seleção** por sinais de g; fluxo contínuo na troca H4 ↔ H7.
  - **Conservação:** sem piso, clipping, reservatório, fluxo imposto nem reutilização da E2.
  - **Preservação:** OFF é byte-idêntico a T2, Bhagwat e T1.
  - **Testes:** 26 asserções de fase ausente e reaparecimento.
  - **Tolerância:** tol da face 1e-11, com piso de ruído em ~1e-12.
- **Execuções (P2.1):**
  - Na réplica do estado aceito, a T1 sem extensão não aceita nenhum passo.
  - Com a extensão, o SS1 **atravessa o selo**:
    1. 14 s em H7, com compressão do gás no oleoduto (p_in 290 → 298 kPa);
    2. reaparecimento do gás em G_g* = 0⁺;
    3. produção de líquido;
    4. **início do blowdown** em 31,5 s de réplica, ≈ 85 s absolutos (α_saída 0,005 → 0,23).
  - Dois binários com algoritmos de raiz distintos concordam a ≤ 3 Pa.
  - Seis balanços ≤ 5,8e-11·M₀.
- **Limites:**
  - O ciclo não se fecha antes do guard (COST_STAGNATION/TIMEOUT).
  - O SS1 desde t = 0 com extensão não chega ao selo em 6 h. Precisa de ≈ 5× mais passos que a T1
    A causa provável é o Jacobiano da face por diferenças, que só tem precisão ~1e-4.
  - Período, amplitude, convergência e Baliño 2010 (P2.2): **NOT_VERIFIED**.

## 1. Certificado da falha em 53,48 s (P0.1)

### 1.1 Reprodução congelada

| Item | Valor |
|---|---|
| Binário | `prims_f_t1b` (fonte `7d44f5d`), SHA-256 `1e1c9fe310d0f8e4…`, GNU Fortran 13.3.0, `-llapack -lblas` → MKL `libmkl_rt.so`, OMP 1 |
| Deck | `D_SS1_C5_BAROTROPIC_N51_dt0.5_T1_t300` (BAROTROPIC, N51+51, Δt_max 0,5 s, `tol_res=1e-10`, estratificado no oleoduto, Chexal no riser, T1) |
| Original (`t1b_short`) | NEWTON_GLOBALIZATION_FAILURE em **53,4753 s**, 152/28 passos |
| Reprodução 1 (`seal_p01_t1b`, `field_interval_s=15`) | idêntica: 53,4753 s |
| Reprodução 2 (`seal_p01_dense`, `field_interval_s=0,25`, só saída) | idêntica: **53,475299860598831 s** |
| Binário atual sem extensão (`prims_f_seal0_538df1e`, critério da face 1e-13) | **COST_STAGNATION**: último passo aceito em **53,473 s** (Δt = 2,8e-3 s), 93 tentativas globais de passo recusadas sem status de face (0), mais recusas da face 366×3614, 79×3613, 21×3615; 7021 s de parede. Mesmo instante físico, mesma assinatura de face |

### 1.2 Último estado aceito e proposta recusada (SI)

Estados das células adjacentes à junção no último estado aceito (`t=53,4753 s`, dump do solver;
sonda `seal_face_probe`, linha `CELL`):

| Célula | θ (°) | U₁ | U₂ | U₃ (kg/(m²·s)) | p (Pa) | α | ρ_g | ρ_l | u_g (m/s) | u_l (m/s) | j_g | j_l | C₀ | V_gj | F₁ | F₂ | F₃ (Pa) | λ₁, λ₂, λ₃ (m/s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 53 (oleoduto) | −2,000 | 438,934 | 1,4113 | 528,375 | 296 603 | 0,55981 | 2,5211 | 997,137 | 0,05461 | 1,2036 | 0,03057 | 0,52981 | — | — | 528,298 | 0,07707 | 297 239 | −244,0; 1,613; 245,1 |
| 54 (riser) | 2,334 | 990,267 | 0,017391 | 524,916 | 296 970 | 0,006890 | 2,5242 | 997,137 | 0,41795 | 0,53007 | 0,002879 | 0,52642 | 0,14220 | 0,34268 | 524,909 | 0,007268 | 297 248 | −205,4; 0,4635; 206,2 |

Face T1 do solver no estado aceito: status 0, G_l* = 528,099, **G_g* = 3,78e-7 kg/(m²·s)**,
α_oleoduto* = 0,5598, **α_riser* = 4,36e-7**, p_L* = 296 671 Pa, p_R* = 297 026 Pa (salto 356 Pa).

Proposta recusada (`last_rejected_iteration.txt`): iteração de Newton 6, 18 retrocessos,
λ = 7,63e-6, status de proposta 3615 (face: falha de avaliação), norma NaN. Recusas da execução:
1080×1615, 874×1613, 389×1614 (face sem passo admissível / sem convergência / avaliação).

Réplica a partir do estado aceito (`seal_replay/…_REPLAY53.4753_T1`, `prims_f_t1b`): **nenhum passo
aceito**; 16 passos recusados com Δt de 1e-3 até **3,05e-8 s**; 1489 recusas da face (685×3615,
623×3613, 180×3614). A falha **não depende do passo de tempo**.

### 1.3 Mapa local independente

Sondas independentes do solver (só compartilham os fechamentos): `seal_face_probe.f90` (curvas
integrais RK4 em ln α e em p) e `seal_config_probe.f90` (curva cinemática em z = √α com U₂
reconstruído da EOS). Para cada α_R* fixo, resolve-se a continuidade de G_l e F₃ em (p_L*, p_M)
e mede-se o descasamento g(α_R*) = F_g(U_L*) − F_g(U_R*).

No estado aceito de 53,4753 s (`probe_t1b_53.48_w0.log`, `conv_t1b_53.48_w0.log`):

| α_R* | g, curvas convergidas (64–256 passos) | g, 8 passos (como o solver) |
|---:|---:|---:|
| 2,8e-2 | −1,59 | — |
| 1e-3 | −5,9e-2 | −5,9e-2 |
| 1e-5 | −7,2e-3 | −6,6e-3 |
| 4,36e-7 (raiz do solver) | **−6,69e-3** | −2,7e-3 |
| 1e-9 | −6,67e-3 | — |

- g é negativo e monótono em todo α_R* ∈ [1e-13; 2,8e-2]: **não existe raiz H4**.
- O limite fechado (`seal_config_probe`, α_R* = 0 exato) converge em nz:
  g₀ = −1,90e-3 (8), −6,79e-3 (16), −7,082e-3 (32), −7,1009e-3 (64) kg/(m²·s). O mesmo valor
  sai da parametrização em ln α com 256 passos: −7,10e-3.
- Ondas acústicas linearizadas × integrais: −6,67e-3 × −7,10e-3. O sinal é o mesmo; a linearização
  é de segunda ordem (Δp/p ≈ 2,5e-4).
- **A "raiz" aceita pelo solver em 53,4753 s é artefato de truncamento** da curva cinemática com
  8 passos RK4 sobre ln(α_M/α_R*) ≈ 9,7. RK4 em ln α com h ≈ 1,2 erra 0,8 % por passo, e
  dG_g/dp_L* ≈ −1,1e-3 kg/(m²·s) por Pa amplifica esse erro.
- O espectro nos traços é real e estritamente hiperbólico; o lado do riser tende a
  λ = (−1497; 0,3431; 1499) m/s, ver §3.

### 1.4 Trajetória antes da falha

Estados aceitos salvos pela execução densa e reconstruídos (`make_probe_fixture.py`), avaliados
com a face com extensão (configuração e g₀) e com a face T1:

| t (s) | G_g* T1 | G_g* com extensão | Configuração | g₀ (kg/(m²·s)) |
|---:|---:|---:|---|---:|
| 50,45 | 5,104e-2 | 5,104e-2 | H4 | +2,16 |
| 50,95 | 4,190e-2 | 4,190e-2 | H4 | +1,86 |
| 51,45 | 3,290e-2 | 3,290e-2 | H4 | +1,46 |
| 51,95 | 2,413e-2 | 2,413e-2 | H4 | +1,09 |
| 52,45 | 1,564e-2 | 1,564e-2 | H4 | +0,738 |
| 52,95 | 7,533e-3 | 7,533e-3 | H4 | +0,390 |
| 53,20 | 1614 no estado reconstruído | 3,888e-3 | H4 | +0,210 |
| 53,45 | 1614 no estado reconstruído | 3,49e-4 | H4 | +0,0201 |
| 53,4753 | 3,8e-7 (espúria) | **0** | **H7** | **−6,6e-3** |

Até 53,45 s a raiz T1 era genuína e as duas faces coincidem. **g₀ troca de sinal entre 53,450 e
53,475 s.** Em 53,20 e 53,45 s a T1 sem extensão nem inicializa nos estados reconstruídos de
`fields.csv`: diferenças de arredondamento bastam para derrubá-la, o que confirma o mau
condicionamento perto do selo.

### 1.5 Classificação

| Hipótese | Veredito | Evidência |
|---|---|---|
| Falha de busca (Newton da face) | **não é a causa primária** | o mapa independente não tem raiz H4 para nenhum α_R* > 0 com curvas convergidas |
| Falha de globalização do Newton global | consequência | não há fluxo de face consistente na vizinhança; Δt até 3e-8 s não ajuda |
| Perda de diferenciabilidade do Chexal | **não** | limite α → 0⁺ regular e C¹, com espectro hiperbólico (§3) |
| Perda de hiperbolicidade | **não** | espectros reais nos traços e na célula |
| Limite constitutivo do Chexal | indireto | a única forma de manter H4 seria G_g* < 0 com líquido subindo na base do riser, fora do domínio a 2,3° (status 61) |
| **Ausência local da fase gasosa na face** | **causa** | a solução local física é a ponte líquida (H7): G_g* = 0, traços da face só com líquido, frente de acumulação no oleoduto com s ≤ 0; a T1 não contém essa configuração |
| Solução fisicamente inadmissível | não | H7 satisfaz Rankine–Hugoniot e as condições características (§4) |

Achado adicional, que afeta o T1 e qualquer solver de face: o fechamento estratificado do oleoduto é
só C⁰ na troca laminar/turbulento da tensão de parede do gás (Re_g = 2000, `fanning_wall_shear`).
Ali o Jacobiano A e os autovetores são descontínuos, e uma curva integral RK4 com número fixo de
estágios fica **descontínua no parâmetro p** quando um estágio cruza Re_g = 2000. Observado: salto
de U₃ de 6e-5 e troca de sinal de F_g numa variação de 1,8e-8 em p_L*, com estagnação do Newton
interno. Além disso, u_g é mal condicionado em U₃: du_g/dU₃ ≈ 1/U₂ ≈ 700.

**Gate P0.1: PASS.** Causa demonstrada: perda da configuração de transmissão de gás (H4) no início
do selo. Ela foi mascarada por um passo pela raiz espúria da curva grosseira da T1.

## 2. Últimos estados aceitos e propostas recusadas

Arquivos congelados (resultados fora do repositório, `~/Documentos/PRIMS-F-chexal-riser-results-20261121/`):

- **Último estado aceito e tentativa recusada:** `t1b_short/runs/D_SS1_…_T1_t300/last_rejected_state.dat`.
  Colunas U_old e U_trial das 102 células; t_nd = 2,2589937 ↔ 53,4753 s.
- **Fixture offline:** `seal_fx/SS1_T1_53.48_DUMP`.
- **U inicial da réplica:** `seal_replay/U_SS1_T1_accepted_53.4753s.txt`, SHA-256 `b5f575f14d0503ce…`.
- **Propostas recusadas:** `rejected_proposals.csv` de cada execução, com status composto da proposta
  (2000 + status da face), estágio (FACE/CELL/BOUNDARY/NONE), λ e normas.

Na tentativa recusada a configuração da junção é a mesma do estado aceito: g₀ = −7,1e-3 e H7 com
s = −8,6e-5 m/s. O passo de Newton era quase nulo (λ = 7,6e-6), e a recusa veio da face T1
(3614/3615), não do resto do domínio.

Resíduo local, no estado aceito:
- **T1:** o mapa independente dá |g| ≥ 6,7e-3 kg/(m²·s) para todo α_R* > 0 (§1.3).
- **Extensão:** o problema H7 converge com resíduo interno < 1e-11 relativo.

Resíduo global do último passo aceito: 3,75e-8 (adimensional, `last_rejected_state.dat`,
cabeçalho), com `tol_res` = 1e-10 do critério de aceitação.

## 3. Limite Chexal–Lellouche quando α_g → 0 (P0.2)

Equações implementadas (`module_chexal_1992.f90`, ar–água), com
u_g = C₀(j_l + j_g) + V_gj e j_g = α·u_g:

- **L(α):** L = 1,15 α^0,45 no ramo vertical (Eq. 3-7); L = 1,125 α^0,6 no ramo horizontal (Eq. 3-19).
- **C₀:** C₀ = L/(K₀ + (1 − K₀) α^r) (Eq. 3-4a); no horizontal, multiplicado por (1 + α^0,05 (1 − α)²) (Eq. 3-18).
- **V_gj:** V_gj = 1,41 [(ρ_l − ρ_g) σ g/ρ_l²]^0,25 C₂ C₃ C₄ C₉, com C₉ = (1 − α)^B₁.
- **Interpolação no ângulo:** fr = (90 − θ_Chexal)/90 (Eqs. 3-2, 3-21). Na base do riser θ = 2,334°, ou seja θ_Chexal = 87,67° e fr = 0,026.
- **Re (Eq. 3-15):** Re_g só vence quando Re_g > Re_l ou j_g < 0. Com j_l = 0,526 m/s, Re_l ≈ 3,1e4, enquanto Re_g ≈ 7,1e3·j_g. Não há troca de ramo para j_g → 0.

Limites (derivação conferida numericamente, `chexal_alpha0_limit.log`):

1. **C₀ → 0** como α^0,6 → α^0,45. O expoente local medido cai de 0,61 para 0,52 em α ∈ [1e-3; 1e-12]. dC₀/dα é singular (~α^−0,4 a α^−0,55).
2. **V_gj → V_gj(0) = 0,34419 m/s**, contínuo: em α = 0 exato, C₀ = 0 e V_gj = 0,3441864.
3. **u_g → V_gj(0)**, então j_g = α u_g ∝ α (expoente local → 1,0000). No ramo co-corrente ascendente, **j_g = 0 ⇔ α = 0**.
4. **Gás presente e parado** (α > 0, u_g = 0) exigiria j = −V_gj/C₀ < 0, isto é, líquido descendo. É contracorrente, fora do domínio a 2,3° (`source_coordinates`, status 61).
5. **j_g → 0⁻ com j_l > 0:** CHEXAL_DOMAIN (61). O gás reverso na base inclinada do riser não é suportado pela correlação.
6. **Continuidade do fluxo:** F_g = ρ_g α u_g é C¹ em α = 0, com dF_g/dα(0⁺) = ρ_g V_gj(0). Não é C², porque a segunda derivada vai como α^−0,4.
7. **Jacobiano:** α·du_g/dα ~ α^0,45 → 0, então A = dF/dU tem limite finito. No traço do riser, λ → (−1497; 0,3431; 1499) m/s, com λ₂ → V_gj(0)·(1 + O(α^0,45)). O sistema é **estritamente hiperbólico no limite**.
8. **Velocidade da fase ausente:** em α = 0, `chexal_recover_velocities` e `compute_primitive_variables` gravam u_g = 0 por convenção de armazenamento, diferente do limite V_gj(0). O fluxo α u_g não muda. O Jacobiano em U₂ = 0 não é avaliável: `DFduj` divide por U₂. A extensão **nunca** avalia fechamento nem Jacobiano em α = 0, só o fluxo de líquido puro.
9. **U₂ conservativo × α:** U₂ = α ρ_g(p). A curva cinemática é integrada em z = √α sobre (U₁, U₃), com U₂ = z² p(U₁, α)/C_g² exato pela EOS (ρ_l = U₁/(1 − α) dá p). Todo estágio tem U₂ > 0, e o extremo α = 0 tem campo nulo (dU/dz = 2z·dU/dα, com dU/dα finito), sem chamada ao fechamento.
10. **T1 no limite:** a configuração H4 permanece bem definida em α_R* = 0 (G₀). Quando g₀ < 0 o limite não fecha o balanço de gás, e a solução física é outra configuração (§4).

**Gate P0.2: PASS.** Limite regular, contínuo, C¹ e hiperbólico do lado do riser. Fase ausente na
face exige mudança de configuração, não um piso de α.

## 4. Formulação da extensão (P1.1)

Opção `junction_seal = 'OFF'` (default) ou `'LIQUID_BRIDGE'`.

- **Restrições:** experimental; exige `junction_transmission='T1'` e o par de fechamentos (1, 5);
  `junction_seal_face_tol` ∈ [1e-14; 1e-9], default 1e-11.
- **Registro:** contrato e cabeçalho de checkpoint gravam a opção **só quando ativa**, e
  `junction_transmission.csv` ganha as colunas `configuration, front_speed_m_s, g0_kg_m2_s`.
- **Preservação:** com OFF, T1 e T2 continuam idênticos (§6).

**Configurações** (todas as ondas 1 e 3 linearizadas ao longo de r(U_célula), em p; onda 2 do riser
pela curva integral em z, RK4 com nz = max(16, min(32, ⌈32·|Δz|/max(z, z₁)⌉)) passos, isto é,
proporcional ao comprimento relativo do trecho):

- **H4 (gás passa):** onda 1 no oleoduto | camada estacionária | ondas 2 e 3 no riser. Incógnitas
  (p_L*, z_R*, p_M). F_L(U_L*) = F_R(U_R*) nas três componentes. Admissibilidade:
  λ₁(U_L*) < 0 < λ₂(U_L*), λ₂(U_R*) > 0 e λ₃(U_M) > 0.
- **H7 (ponte líquida):** onda 1 no oleoduto até U_L1; frente de acumulação (onda 2 do oleoduto)
  de U_L1 até um estado de líquido puro U*; U* é o extremo α = 0 da curva cinemática do riser a
  partir de U_M = curva 3(U_R; p_M).
  - **Rankine–Hugoniot da frente:** gás s(0 − U₂,L1) = 0 − F₂,L1, logo **s = u_g(U_L1)**;
    líquido s(U₁* − U₁,L1) = U₃* − F₁,L1; momento s(U₃* − U₃,L1) = F₃(U*) − F₃,L1.
  - **Incógnitas:** (p_L1, p_M).
  - **Fluxo de face:** F(U*) = (U₃*, 0, F₃*).
  - **Admissível se:** s ≤ 0 (a frente fica no oleoduto); λ₂(U_L1) > s (o gás entra na frente);
    λ₁(U_L1) < 0; λ₂(U_M) > 0 e λ₃(U_M) > 0 (rarefação cinemática do riser de α = 0 a α_M).

**Seleção só por sinais:**

1. Avalia-se g(z_c) na célula do riser. Se g(z_c) ≥ 0, a raiz está acima: H4.
2. Se g(z_c) < 0, o estágio G₀ resolve a H4 em α_R* = 0, e decide o sinal de g₀:
   - g₀ > 0: raiz H4 em (0, z_c);
   - g₀ ≤ 0: H7.
3. Se |g₀| ≤ tol·max(|F₁|, |g₀|), z = 0 já resolve as equações dentro da tolerância de
   convergência. Esse é o limite comum H4 = H7 (s = 0), devolvido com G_g* = 0.

**Continuidade:** em g₀ = 0 tem-se s = u_g(U_L1) = 0, U_L1 = U_L* e as equações da H7 se reduzem às
da H4 em α_R* = 0. O fluxo é contínuo na troca: verificado com salto O(afastamento), §5.

**Raízes:**

- **Problema interno** (z fixo): Newton 2×2 em (p_L* ou p_L1, p_M) com busca linear. Convergência
  por |r| < 1e-14, ou por |r| < tol com estagnação. Falha da busca linear só é aceita abaixo de tol
  (piso de ruído). Reinícios determinísticos: chute anterior, pressões das células, solução de G₀.
  Status: 1622 (singular), 1623 (sem passo amortecido), 1624 (sem convergência).
- **Problema externo em z:**
  - se g(z_c) ≥ 0, secante salvaguardada para cima (extrapolação com 1 % de folga; em falha
    interna recua para o último z válido; 1627 se não houver bracket);
  - se g₀ > 0, Illinois no bracket (0, z_c);
  - se G₀ sai do domínio do Chexal, z_c é reduzido por fatores de 4, até 6 vezes.
- **Aceitação:** |g| ≤ tol·max(|F_g|, |F₁|); se g(z_c) já está nesse intervalo, a raiz é a própria célula.
- A versão anterior (`68608e1`, `e5cf228`) usava Newton 3×3 em (p_L*, z, p_M). Em estados com muito
  gás ele não convergia e custava ~98 s por passo global. Foi substituído em `46e8830` (≈ 26 s por
  passo, fluxos iguais dentro de 1e-5, §8).

**Diferenças finitas:** em modo selo, amostras na outra configuração contam como inadmissíveis, e o
estêncil unilateral de 2ª ordem do lado da configuração base dá um elemento consistente do
subdiferencial B (Newton semissuave).

**Escolhas descartadas, com causa:**

- **Ondas acústicas por curva integral RK4:** o fechamento estratificado é C⁰ em Re_g = 2000 e a
  curva discreta fica descontínua em p.
- **Mismatch de gás medido relativo a |F_g|:** o ruído dos fechamentos é relativo ao fluxo da
  mistura.
- **Busca que exigia G₀ sempre:** para estados com muito gás, a curva até α = 0 sai do domínio do
  Chexal.
- **Fluxo H4 avaliado no traço do oleoduto:** G_g ali é amplificado por du_g/dU₃ ≈ 1/U₂.

**Não reutiliza a E2 (Bhagwat).** Não há piso de α, clipping, reservatório nem G_g imposto. O gás
que chega do oleoduto fica nas células do oleoduto, que pertencem ao domínio distribuído.

**Gate P1.1: PASS** (formulação conservativa, admissibilidade explícita, comportamento definido nos limites).

## 5. Testes de fase ausente e reaparecimento (P1.2)

`tests/test_junction_seal.f90` (26 asserções, PASS; estados SI do SS1 em 53,4753 s):

| # | Caso | Resultado |
|---|---|---|
| 1 | OFF: operador T1 puro, sem configuração de selo | PASS (T1 dá 1613 nesse estado) |
| 4 | Estado de selo do SS1 | H7, G_g* = 0 exato, s = −9,15e-5 m/s, g₀ = −7,5e-3, um único fluxo, traços α* = 0 |
| 5 | Fluxo de gás → 0 e reaparecimento (parâmetro U₃ do oleoduto) | troca em U₃c = 528,4182 kg/(m²·s); H4 de um lado e H7 do outro; salto do fluxo 1,59e-5 (afastamento 1e-5) e 1,52e-7 (1e-7), linear; G_g* (lado H4) = 1,58e-5 → 1,52e-7 → 0⁺; sem gás reverso |
| 2 | Gás presente, co-corrente | H4, G_g* = 5,49895e-4 (T1: 5,49890e-4); máx. diferença relativa do fluxo 5,7e-6, atribuída à curva cinemática da T1 (8 passos em ln α) |
| 6 | Jacobiano × diferenças centradas independentes (incrementos 3e-5 e 3e-6) | H4: 2,0e-4 / 1,2e-4; H7 longe da troca: 7,6e-5 / 7,6e-5 (limite 5e-3, ver §7) |
| 7 | Densidade do gás +20 % | admitido, fluxo único, G_g* = 0 ⇔ H7 |
| 3 | Partição homogênea (estratificado dos dois lados) | H4, \|F̂ − F\|/\|F\| = 7,8e-16, sem salto de pressão |
| 8 | Célula do riser só com líquido (U₂ = 0) | recusa explícita (status 12): Jacobiano do fechamento em U₂ = 0 não suportado |

Valores da tabela: binário de teste compilado de `46e8830` (log `tmp/unit_tests/logs/test_junction_seal.log`,
22/11/2026 07:25).

Varredura de 51 estados em torno do selo, com três parâmetros (momento do oleoduto, gás do
oleoduto, gás do riser), executada com a versão intermediária (Newton 3×3 e nz adaptativo com mínimo
4), **não repetida** com o binário final: **51/51 convergidos**. H7 com G_g* = 0 quando g₀ < 0 e H4 com G_g* > 0
quando g₀ > 0, com G_g* de 4,2e-3 a 1,1e-1 kg/(m²·s). Antes das correções do §4, 4 dos 51
falhavam (1623).

**Não testado / fora do escopo** (P1.2 itens 3, 5, 8):

- **Gás presente e estacionário no riser:** fora do domínio do Chexal (§3.4).
- **Fase exatamente ausente numa célula:** recusada explicitamente.
- **Gás reverso na base do riser:** não suportado pela correlação a 2,3°.

## 6. Preservação de T1/T2 e regressões

Binário final `prims_f_seal5_46e8830` (SHA-256 `e9af384c98aebc33…`). As três comparações foram feitas
com `seal4_e5cf228` e **repetidas com `seal5_46e8830`**, com o mesmo resultado byte a byte
([seal5_vs_reference_compare.txt](../validation/diagnostics/t1_seal_20261121/evidence/regression/seal5_vs_reference_compare.txt): 22/22 arquivos idênticos):

| Comparação | Resultado |
|---|---|
| T2 SS1-Chexal, 30 s: `prims_f_t1` (dc92c80) × extensão OFF | **idêntico byte a byte**: monitor, conservação, campos, contorno, recusas, checkpoint e steps (col. 1–12). No contrato muda só `input_file` |
| T2 Stable-Bhagwat (C3, ISOTHERMAL), 30 s | **idêntico** (idem) |
| T1 SS1, 5 s: `prims_f_seal0_538df1e` (T1 antes da extensão) × extensão OFF | **idêntico**, inclusive `junction_transmission.csv` |
| Suíte Fortran | 31/31 PASS (`test_junction_seal`: 26 asserções) |
| Suíte Python | 123 OK () |
| Suítes CLI | 8/8 PASS (contrato: 11 casos, 4 recusas; caso com selo encurtado para 0,01 s por custo) |
| Builds | release e `check_all` (`-std=f2018 -fcheck=all -ffpe-trap=invalid,zero,overflow`) sem erros |

## 7. Jacobianos e tolerâncias do Newton da face (P1.3)

**Tolerância da face com extensão** (`seal_tol_study.f90`, 10 chamadas por estado):

| tol | H4 longe do selo (51,95 s) | H4 perto do selo (53,45 s) | H7 (53,475 s) | CPU por face |
|---|---|---|---|---:|
| 1e-10 | 0 | 0 (dif. 5,9e-12 em G_g, 2,8e-8 Pa) | 0 | 0,39–0,52 s |
| **1e-11** (default) | 0 | 0 (idem) | 0 | 0,38–0,54 s |
| 1e-12 | 0 | 0 (referência) | 0 | 0,38–0,56 s |
| 1e-13 | **1623** | **1623** | **1623** | 1,6–2,2 s |

**Leitura da tolerância:**

- O piso de ruído de cada avaliação de fluxo é ~1e-12 relativo. Vem da tolerância de recuperação
  de velocidades dos fechamentos (`tol` = 1e-12 no deck).
- 1e-13 fica abaixo desse piso, e o Newton interno não consegue reduzir o resíduo.
- Entre 1e-10 e 1e-12 os fluxos coincidem dentro do ruído; o default 1e-11 (10 × `tol`) fica no meio.
- **Mesma explicação para a T1:** com 1e-11 (`t1b`), o SS3 tinha MAX_NEWTON por um piso de 3e-10 no
  resíduo global. Com 1e-13 (`t1c`), o SS1 acumula recusas 1614 desde o início: 27 antes de 0,5 s
  e 314 até 53,4 s. Nenhum dos dois valores foi escolhido por tempo físico alcançado.

**Jacobiano:**

- A precisão por diferenças é limitada pelo ruído da face. Na coluna da massa de gás (U₂ ≈ 1e-3),
  tol·|F₁|/(1e-6·U₂) ≈ 4e-3. É precisão de Newton inexato, documentada.
- Junto à troca H4/H7 o fluxo é C⁰ (quina). Diferenças centradas misturavam as duas derivadas
  laterais: no estado de 53,4753 s, dF₁/dU₁,R valia −128 (misto) contra −203,2 (H7).
- Na réplica, o Newton global precisava de 17–25 iterações por passo. Com as diferenças
  unilaterais do lado base (`e5cf228`) passou a 3–6 e o passo voltou a crescer.

**Gate P1.3: PASS**, com a limitação explícita de precisão do Jacobiano (~1e-3 na coluna de U₂).

## 8. Replays e SS1 desde t = 0 (P2.1)

Execuções com o deck `D_SS1_C5_BAROTROPIC_N51_dt0.5` (BAROTROPIC, N51+51, Δt_max = 0,5 s,
`tol_res = 1e-10`), isoladas pelo guard `run_guarded.py`:

- relógio monotônico, que não conta suspensão da máquina;
- teto de 21 600 s de parede;
- COST_STAGNATION quando a simulação avança menos de 1 s em 1800 s de parede.

O guard encerra a execução com sinal: não há `run_status.txt`, e o `summarize_campaign.py`
classifica essas execuções como `FAIL:NO_STATUS/0`. A classificação válida é o veredito do guard.

Até 6 execuções rodaram ao mesmo tempo numa máquina de 8 núcleos, então **os tempos de parede não
são comparações controladas de custo**.

Nas réplicas, t = 0 corresponde a 53,4753 s do SS1 (estado aceito congelado, SHA-256 `b5f575f1…`).

### 8.1 Três operadores, sem mistura

| Operador | Binário | Execução | Resultado | t alcançado | Passos aceitos / recusados | Propostas recusadas (status) |
|---|---|---|---|---:|---|---|
| T2 legado (Chexal no riser) | — | [relatório T1, §10](PRIMSF_T1_JUNCTION_QUALIFICATION.md) | BLOCKED_JUNCTION_TRACE_DOMAIN (3061) | **52,94 s** | — | a T2 exige G_g* < 0 na base a 2,3° |
| T1 sem extensão | `prims_f_t1b` (7d44f5d) | `t1b_short` | NEWTON_GLOBALIZATION_FAILURE | **53,4753 s** | 152 / 28 | 1080×1615, 874×1613, 389×1614 |
| T1 sem extensão | `prims_f_seal0_538df1e` | `seal_p01` | COST_STAGNATION | **53,473 s** | 143 / 2 | 366×3614, 79×3613, 21×3615, 93×0 |
| T1 sem extensão, réplica | `prims_f_t1b` | `seal_replay/…_T1` | NEWTON_GLOBALIZATION_FAILURE | **0 s** (nenhum passo; Δt até 3,05e-8 s) | 0 / 16 | 685×3615, 623×3613, 180×3614 |
| **T1 + LIQUID_BRIDGE, réplica** | `seal4_e5cf228` | `seal_replay/…_T1SEAL` | COST_STAGNATION | **33,91 s** (87,39 s abs.) | 512 / 112 | 4534×0, 13×3625, 3×3627 |
| **T1 + LIQUID_BRIDGE, réplica** | `seal5_46e8830` | `seal5_replay` | TIMEOUT (6 h) | **31,85 s** (85,32 s abs.) | 485 / 104 | 4318×0, 23×3625, 5×3627 |
| T1 + LIQUID_BRIDGE, réplica, tol da face 1e-12 | `seal5_46e8830` | `seal5tol12_replay` | COST_STAGNATION | 14,02 s | 243 / 56 | 2361×0, 33×3625, 19×3623, 14×3627 |
| **T1 + LIQUID_BRIDGE, SS1 desde 0** | `seal4_e5cf228` | `seal_ss1` | TIMEOUT (6 h) | **43,36 s** | 649 / 154 | 6143×0 |
| T1 + LIQUID_BRIDGE, SS1 desde 0 | `seal5_46e8830` | `seal5_ss1` | TIMEOUT (6 h) | 27,61 s | 417 / 94 | 3850×0 |
| T1 + LIQUID_BRIDGE, SS1 desde 0, tol 1e-12 | `seal5_46e8830` | `seal5tol12_ss1` | COST_STAGNATION | 2,02 s | 44 / 4 | 199×0 |

Status 0 = a proposta foi recusada pela globalização do Newton global (busca linear sem
decréscimo suficiente), sem falha de avaliação da face.

**Nenhuma execução com extensão chegou a 53,48 s desde t = 0**, porque o custo é alto. A travessia do
selo só foi demonstrada pelas réplicas, que partem do estado aceito exato.

### 8.2 Física da réplica (seal4; a seal5 coincide, §8.3)

| Fase (t da réplica) | Configuração | Observado |
|---|---|---|
| 0 → 14,03 s | **H7** (105 faces amostradas) | G_g* = 0 exato. A frente de acumulação fica no oleoduto com s ≤ 0: de −8e-5 m/s até **−0,19 m/s em 8,2 s**, depois volta a ≈ −0,035 m/s. α na base do riser → 0 (1e-8 em 10 s); α no fim do oleoduto 0,56 → 0,50. p_in 290,1 → 297,1 kPa. M_l 79,64 → 81,39 kg; o gás que entra (2,87e-4 kg/s) fica retido no oleoduto: M_g,oleoduto 0,1737 → 0,1777 kg, M_g,riser 4,4e-3 → 7,7e-4 kg |
| 14,03 s | troca H7 → H4 | g₀ troca de sinal; gás volta a penetrar no riser com G_g* = 0⁺, sem salto de fluxo e sem gás reverso |
| 14 → 18,6 s | H4 | α na base do riser cresce, máximo de **p_in = 297,94 kPa em 18,55 s** |
| 18,6 → 31,5 s | H4 | o gás sobe pelo riser: α_base 0,076 → 0,30, ṁ_l,saída 1,23 → 2,0 kg/s, p_in cai para 284 kPa e M_l para 76 kg; α_saída ainda ≤ 0,005 |
| 31,5 → 33,9 s | H4 | **o gás chega ao topo (blowdown):** α_saída 0,005 → 0,23, ṁ_l,saída até 2,44 kg/s em 32,39 s, p_in 272,5 kPa em 33,91 s, M_l 72,6 kg |

A sequência é a do golfão severo:

1. bloqueio da base;
2. acumulação de líquido no riser com compressão do gás no oleoduto;
3. penetração do gás;
4. produção de líquido;
5. início do blowdown.

O modelo T1 sem extensão não passa da fase 1. **O ciclo não se fecha**: as execuções param no início
do blowdown, sem nova selagem. Período e amplitude não podem ser medidos.

### 8.3 Reprodutibilidade entre binários e comparação com a T1

- **Réplica seal4 × seal5** (algoritmos de raiz diferentes, mesma formulação):
  - p_in difere ≤ 3 Pa até 31 s e 10 Pa em 31,7 s, já no início do blowdown;
  - α na base do riser difere ≤ 1e-4;
  - troca H7 → H4 em 14,027 s × 14,042 s.
  
  O resultado não depende do algoritmo de raiz.
- **SS1 desde 0, seal4 × seal5:** p_in difere ≤ 1,2 Pa até 27,6 s. Nesse trecho só há H4 (318/318 faces
  amostradas).
- **SS1 desde 0, extensão × T1 (`seal0`), operadores distintos:** a diferença em p_in cresce de
  1,6 Pa (1 s) para 34 Pa (20 s), 121 Pa (30 s) e **273 Pa (43 s)**; α_base difere até 8,6e-3. Mesmo
  só em H4, os dois operadores não são iguais:
  - a T1 integra a curva cinemática com 8 passos RK4 em ln α e as ondas acústicas por curva integral;
  - a extensão usa 16–32 passos em z e ondas acústicas linearizadas;
  - na mesma face co-corrente o fluxo difere 5,7e-6 relativo (§5).
  
  O desvio acumulado é da ordem do erro de truncamento da T1 já documentado em §1.3. **Os resultados
  dos dois operadores não são intercambiáveis.**

### 8.4 Custo

CPU (s) até o mesmo instante simulado; carga concorrente variável:

| Execução | Intervalo | Passos | Newton total | CPU | CPU/passo | Δt médio |
|---|---|---:|---:|---:|---:|---:|
| T1 sem extensão (`seal0`) | SS1 0 → 27,6 s | 87 | 296 | 2933 | 33,7 | 0,32 s |
| extensão seal4 | SS1 0 → 27,6 s | 423 | 2044 | 14 494 | 34,3 | 0,065 s |
| extensão seal5 | SS1 0 → 27,6 s | 416 | 2248 | 21 501 | 51,7 | 0,066 s |
| extensão seal4 | réplica 0 → 31,7 s | 484 | 2809 | 13 398 | 27,7 | — |
| extensão seal5 | réplica 0 → 31,7 s | 483 | 2856 | 21 108 | 43,7 | — |

- **Mais passos que a T1:** a extensão precisa de ≈ 5× mais passos no SS1 em regime puramente H4. As
  recusas são de globalização (status 0), não da face.
  - Causa provável: o Jacobiano da face por diferenças é preciso só a ~1e-4 (§5, §7) contra ~1e-8 na
    T1. O Newton global fica inexato e a busca linear recusa mais passos.
  - O Δt é então limitado pelo controlador de recusas, não pela física.
- **seal5 não ficou mais barato que seal4 em produção**, apesar do custo menor por face nos estados
  rígidos isolados (§4). Sob carga concorrente, a CPU por passo foi 1,5× maior.
- **Conclusão de custo: inconclusiva quanto ao binário, negativa quanto à extensão em H4.** Fica
  registrado como limitação; não foi ajustado nenhum parâmetro para acelerar.

## 9. Conservação e convergência

### 9.1 Seis balanços de massa

Fonte: `scripts/analyze_accepted_conservation.audit_ledger`, que integra os fluxos **aceitos**.
São duas fases × sistema/oleoduto/riser; a face da junção entra nos balanços do oleoduto e do riser.
Arquivo: [SEAL_BALANCES_final.md](../validation/diagnostics/t1_seal_20261121/evidence/SEAL_BALANCES_final.md).

| Execução | t alcançado (s) | max \|balanço\|/M₀ (6 balanços) | Ledger |
|---|---:|---:|---|
| réplica seal4 (H7 + H4 + início do blowdown) | 33,91 | 5,8e-11 | PASS |
| réplica seal5 | 31,77 | 5,9e-12 | PASS |
| réplica seal5, tol da face 1e-12 | 14,01 | 7,0e-14 | PASS |
| SS1 seal4 desde 0 | 43,36 | 1,1e-11 | PASS |
| SS1 seal5 desde 0 | 27,61 | 6,0e-12 | PASS |
| SS1 seal5, tol 1e-12 | 1,96 | 2,8e-12 | PASS |
| SS1 T1 sem extensão (`seal0`), referência | 53,42 | 8,8e-13 | PASS |

Nas fases H7 o fluxo de gás da face é **zero exato**. O gás que chega do oleoduto fica nas células do
oleoduto, e M_g,oleoduto cresce como a integral de ṁ_g,in (§8.2). Os balanços fecham sem fonte,
sumidouro ou reservatório.

O resíduo de 5,8e-11 da réplica seal4 está quatro ordens abaixo do `tol_res` relativo e é coerente
com a acumulação de 512 passos.

### 9.2 Quantidade de movimento

- **Face:** o fluxo de face é um único vetor (F₁, F₂, F₃) usado pelas duas células. A identidade
  discreta de conservação na face (`t1_face_check`, linha BAL: soma dos fluxos dos dois lados menos o
  fluxo único) foi medida em todos os estados da trajetória 50,45–53,4753 s, com e sem extensão.
  - Resíduo máximo: 1,9e-15 na componente de momento (escala 4,48, adimensional) e ≤ 2e-18 nas componentes de massa.
  - Na H7 (53,4753 s): 1,4e-16.
- **Rankine–Hugoniot da frente H7:** as três relações, inclusive a de momento
  s(U₃* − U₃,L1) = F₃* − F₃,L1, entram no sistema interno com resíduo < tol (1e-11 relativo).
- **Ledger de momento em nível de execução:** não existe. O momento não é conservado globalmente,
  porque gravidade, atrito e termos geométricos são fontes. **NOT_VERIFIED.**

### 9.3 Convergência (P2.2)

**NOT_VERIFIED.** Nenhuma execução produziu um ciclo completo: as réplicas param no início do
blowdown, e o SS1 desde 0 não chega ao selo dentro do teto de 6 h. Sem trajetória completa não há
período nem amplitude para refinar em Δt ou malha, nem para comparar com Baliño et al. (2010).

## 10. Comparação experimental

Baliño, Burr e Nemoto (2010): o PDF não está no ambiente. Além disso, a convergência depende de
haver uma trajetória completa (§9). **NOT_VERIFIED.** Nenhuma qualificação Bhagwat/T2 foi transferida.

## 11. Matriz de status

| Item | Status | Evidência / restrição |
|---|---|---|
| **P0.1** Causa da falha em 53,48 s | **PASS** | Ausência local da fase gasosa na face (perda da configuração H4 no início do selo). Não é degeneração do Chexal nem perda de hiperbolicidade. O Newton falha porque não há raiz H4; a "raiz" de 53,4753 s é artefato da curva de 8 passos (§1) |
| **P0.2** Limite Chexal α → 0 | **PASS** | C₀ → 0, V_gj → 0,34419 m/s, F_g C¹ (não C²), espectro hiperbólico (§3) |
| **P1.1** Extensão `LIQUID_BRIDGE` | **PASS (experimental)** | Opcional, só T1 + (1, 5); seleção por sinais; sem piso, clipping, reservatório, fluxo imposto nem E2 (§4) |
| **P1.2** Fase ausente e reaparecimento | **PASS**, com itens fora do domínio | 26 asserções; troca contínua com salto O(afastamento). Não suportados e recusados explicitamente: gás estacionário, célula sem gás, gás reverso a 2,3° (§5) |
| **P1.3** Tolerância do Newton da face | **PASS**, com limitação | 1e-10…1e-12 equivalentes; 1e-13 abaixo do piso de ruído (1623). Default 1e-11; `tol_res = 1e-10` mantido. Jacobiano por diferenças preciso a ~1e-4 (§7) |
| Preservação T2 / Bhagwat / T1-OFF | **PASS** | Byte a byte com seal4 e seal5 (§6) |
| Builds debug (`check_all`) e release | **PASS** | release reproduz o binário final byte a byte (SHA-256 `e9af384c…`); `check_all` (`-std=f2018 -fcheck=all -ffpe-trap=invalid,zero,overflow`) sem erros; [final_verify/](../validation/diagnostics/t1_seal_20261121/evidence/final_verify/) |
| Suítes Fortran / Python / CLI | **PASS** | Fortran 31/31; Python 123 OK com `PYTHONPATH=scripts` (sem ele, 3 módulos de `scripts/` não importam); CLI 8/8 |
| Seis balanços de massa | **PASS** | ≤ 5,8e-11·M₀ em todas as execuções (§9.1) |
| Momento | **PASS** (face e frente) / **NOT_VERIFIED** (ledger global) | Identidade da face ≤ 1,9e-15 (§9.2) |
| **P2.1** Réplica a partir de 53,4753 s | **PASS** (travessia do selo) / **BLOCKED** (ciclo) | A T1 sem extensão não aceita nenhum passo. Com extensão: 33,9 s de réplica atravessando H7 → H4 → blowdown, encerrada por COST_STAGNATION |
| **P2.1** SS1 desde t = 0 com extensão | **BLOCKED** (custo) | TIMEOUT em 43,4 s (seal4) e 27,6 s (seal5); ≈ 5× mais passos que a T1 (§8.4) |
| **P2.2** Convergência e Baliño 2010 | **NOT_VERIFIED** | Sem ciclo completo; PDF de Baliño 2010 ausente (§10) |
| Critério final (degeneração × fase ausente × método) | **Respondido** | Fase ausente na face. O método não linear falha em consequência, e a perda de constitutividade é só indireta (gás reverso fora do domínio) |

## 12. Commits e comandos

### 12.1 Commits (branch `implementation/chexal-riser-20261121`, sem push nem merge; `main` = `f3717dd` intacta)

| Commit | Conteúdo |
|---|---|
| `538df1e` | base (relatório T1, rodada 4) |
| `68608e1` | feat: extensão experimental `junction_seal='LIQUID_BRIDGE'` (H4/H7, contrato, checkpoint, CSV, testes) |
| `e5cf228` | fix: diferenças unilaterais da face na configuração base (Newton semissuave); diagnósticos de tolerância e do limite Chexal |
| `1d56984` | evidence: mapas locais do P0.1, sonda de configuração na trajetória, face OFF × selo, estudo de tolerância, limite α → 0 |
| `46e8830` | perf: secante salvaguardada em z, nz proporcional (16–32), raiz no valor da célula; caso CLI encurtado |
| ⟨THIS⟩ | docs/evidence: este relatório, balanços finais, comparação de regressão, handoff |

### 12.2 Binários congelados

Diretório `~/Documentos/PRIMS-F-chexal-riser-results-20261121/bin/`. Compilador GNU Fortran 13.3.0
(Ubuntu 13.3.0-6ubuntu2~24.04.1); `-llapack -lblas` resolvido para MKL `libmkl_rt.so`; `OMP_NUM_THREADS=1`.

| Binário | Fonte | SHA-256 |
|---|---|---|
| `prims_f_t1b` | `7d44f5d` (T1, face 1e-11) | `1e1c9fe310d0f8e4451bd3cf91e09c3db6a1b86070bf1f2fe560296f199aec2a` |
| `prims_f_seal0_538df1e` | `538df1e` (T1, face 1e-13) | `91a41988594dac4108e920ee23a50b0e3c18399f046117dd3a39a6da2326718a` |
| `prims_f_seal3_68608e1` | `68608e1` | `cc4c96139a41b0f1efd4bfb9554dc71a27f7bd532d3176ee66d27c3df1d76403` |
| `prims_f_seal4_e5cf228` | `e5cf228` | `92ff1b7a20aca8a7c8824498dcc71b1792c1487ca351eb304575607713d02785` |
| `prims_f_seal5_46e8830` | `46e8830` (**final**) | `e9af384c98aebc337e1a3f1599b4356f9dae5f0c4e7659c2f45a4333542a8601` |
| `prims_f_seal1`, `prims_f_seal2` | intermediários, **SUPERSEDED** | `ed8ba765…`, `52198524…` |

### 12.3 Entradas

| Arquivo | SHA-256 |
|---|---|
| U inicial da réplica `seal_replay/U_SS1_T1_accepted_53.4753s.txt` | `b5f575f14d0503ce2ca02c7ba5e7088600a6692c1f9d70da497f507f64f48ed9` |
| deck réplica T1 | `59691e54…` |
| deck réplica seal4 | `626e83cc…` |
| deck réplica seal5 | `5179e1ea…` (difere do seal4 só em `out_dir`) |
| deck réplica tol 1e-12 | `0534769b…` (+ `junction_seal_face_tol=1e-12`) |
| deck SS1 seal4 | `05e73dab…` |
| deck SS1 seal5 | `8e00caea…` |
| deck SS1 tol 1e-12 | `48c956db…` |

### 12.4 Comandos

```bash
# build e suítes (worktree no HEAD final)
make -j4 release && make -j4 check_all
FC=gfortran LIBS='-llapack -lblas' python3 scripts/run_unit_tests.py --timeout 900
PYTHONPATH=scripts python3 -m unittest discover -s tests -p 'test_*.py'
for t in tests/*_cli.py; do python3 $t --binary bin/prims_f --output tmp/fv_$(basename $t .py); done

# execuções isoladas pelo guard (teto 21 600 s, janela 1800 s, avanço mínimo 1 s)
python3 validation/diagnostics/junction_constitutive/run_guarded.py BIN DECK.nml RUN_DIR \
    --max-wall-s 21600 --window-s 1800 --min-advance-s 1

# balanços
python3 validation/diagnostics/chexal_riser_20261121/summarize_campaign.py OUT.json OUT.md RUN_DIR...

# sondas offline (validation/diagnostics/t1_seal_20261121/)
#   seal_face_probe, seal_config_probe, seal_tol_study, chexal_alpha0_limit; logs em evidence/
```
