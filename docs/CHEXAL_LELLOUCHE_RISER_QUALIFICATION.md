# Chexal–Lellouche (1992) como fechamento de todo o riser: qualificação

Data: 21/11/2026. Worktree `PRIMS-F-chexal-riser-20261121`, branch
`implementation/chexal-riser-20261121`, base `d2c8101` (E1/E2 experimentais
preservados, desligados). Resultados brutos: `~/Documentos/PRIMS-F-chexal-riser-results-20261121`.
Ferramentas e evidências: [validation/diagnostics/chexal_riser_20261121](../validation/diagnostics/chexal_riser_20261121/README.md).
Sem push, merge ou alteração da `main`.

## Resumo

O Chexal–Lellouche (1992, ar–água, ID 5, ramo AUTO) é agora o fechamento padrão
de todo o riser na branch de desenvolvimento. O Bhagwat (ID 3) continua selecionável
e reproduz byte a byte o binário congelado anterior. Três verificações independentes
mostram que o caminho numérico inteiro usa o Chexal:

- **Equações:** o núcleo Fortran concorda com uma terceira transcrição independente
  das equações em 945 estados (C0 até 2,4e-16, Vgj até 6,9e-14). Essa transcrição
  reproduz a Tabela A-1 do EPRI NP-3989 nos casos de baixa pressão, que são o
  regime do PRIMS-F.
- **Jacobiano global:** o Jacobiano montado do domínio acoplado coincide com
  diferenças finitas do resíduo completo (306 colunas, erro ≤ 4,3e-8) em estados
  aceitos dos quatro casos.
- **Estados aceitos:** todos os ~18 mil estados de riser gravados nas execuções
  Chexal satisfazem a relação Chexal com resíduo ≤ 8,3e-14. Nas execuções Bhagwat,
  usadas como controle, o resíduo é de 6e-2 a 6,6e-1.

O resultado físico da campanha:

- **Stable:** completa 2400 s com as duas EOS.
- **OSC e SS3:** param ao exigir **contracorrente inclinada** no riser (líquido
  descendo com gás subindo, a 68–76° da horizontal). O Chexal 1992 não cobre esse
  estado (p. 262), e a recusa é explícita: `BLOCKED_CONSTITUTIVE_DOMAIN`.
- **SS1:** para quando o fluxo de gás na junção precisa inverter sobre o traço
  horizontal, que a fonte só admite cocorrente (pp. 264–265):
  `BLOCKED_JUNCTION_TRACE_DOMAIN`.

Esses eventos persistem sob refino de malha e passo. Nenhum caso com Chexal produziu
ciclos completos, portanto não há validação experimental. A promoção a produção
continua bloqueada pela lacuna constitutiva de contracorrente inclinada.

## 1. Mapeamento das substituições Bhagwat → Chexal

| Local | Antes | Agora | Natureza |
|---|---|---|---|
| `scr/main.f90`, default do namelist | `closure_riser = 2` (Bendiksen) | `5` (CHEXAL_1992 AUTO) | mudança intencional do default |
| `scr/module_config.f90`, default do tipo | `2` | `5` | idem (testes que constroem `sim_config_t()`) |
| `input/benchmarks/D_*.nml` (8 decks primários) | `closure_riser=3` | `closure_riser=5` | só ID e linha de proveniência; geometria, EOS, alimentação e controles idênticos |
| `input/benchmarks/bhagwat/` | — | 8 decks Bhagwat arquivados byte a byte | baseline de regressão |
| `input/benchmarks/transition_to_chexal.json` | — | hashes arquivados/primários e as 2 linhas alteradas por deck | rastreabilidade |
| `scripts/run_validation.py --closure` | default 3 | default 5 | idem |
| `scripts/run_validation.py assign()` | editava a 1ª ocorrência, inclusive em comentário | só linhas não comentadas | correção de robustez (ver §8) |
| ACTIVE_SET (`admissibility_margin`, `constrained_direction`) | ficava **silenciosamente inativo** com outro fechamento | `error stop` se `closure_riser≠3` | margem derivada só para Bhagwat |
| E1 `junction_flux='PHASE_UPWIND'` | aceito com qualquer fechamento | `error stop` com ID 5–7 | estudado só com Bhagwat |
| E2 `junction_riser_limit='ALPHA0'` | idem | `error stop` com ID 5–7 | idem |
| BLOCKING (projeção de barreira) | só IDs 2/4 | inalterado (com Chexal falha com status explícito) | — |
| `offline_api.f90` (diagnóstico de 19/11) | Bhagwat fixo | inalterado | ferramenta histórica, fora do solver |

Pontos do solver que já usavam o fechamento do riser de forma genérica (`cfg%closure_riser`
ou `cell_closure`) e passam a usar Chexal sem mudança de código:

- recuperação de primitivas (`module_find`, `chexal_recover_velocities` com dica de continuação);
- estado pontual, fontes e Jacobiano (`point_state`, `DFduj`, derivadas Chexal por Richardson);
- matriz de Roe por quadratura no caminho;
- traço do riser na junção MATCHED (`solve_alpha_from_superficials` em θ=0);
- traços de contorno (`module_boundary_trace`, `module_boundary_conditions`);
- perfil estacionário (`EstadoEstacionario_ndim_simp` → `chexal_solve_alpha`).

O estado inicial de cada execução é recalculado pelo solver com o fechamento selecionado;
nenhum estado inicial Bhagwat é reutilizado.

## 2. Equações e domínio físico

Fonte: Chexal, Lellouche, Horowitz e Healzer, *Prog. Nucl. Energy* 27(4) 255–295, 1992
(PDF local, SHA-256 `b0b6b68f…`, página PDF = impressa − 254). Implementação:

- identidade drift-flux `R = α[C0(jl+jg) + Vgj] − jg` (Eq. 3-1 corrigida, sem dividir por j);
- C0 ponderado por Fr (3-2/3-3);
- L ar–água (3-7/3-8) e componente horizontal (3-18/3-19);
- B1, K0 e r com Reynolds assinados (3-11–3-17);
- Vgj = Vgj⁰·C9 com C2, C3 e C4 (3-21–3-38);
- C10 com a parentização do RELAP5/MOD3 1995;
- ramos HIGH/LOW e CCFL (3-39–3-41).

Domínio da fonte, que o código respeita:

| Situação | Fonte | Código |
|---|---|---|
| Cocorrente ascendente/descendente, 0–90° | p. 261–265 | aceito |
| Contracorrente (gás ↑, líquido ↓) | **só vertical**, p. 262 | aceito só em θ_C = 0, com ramo explícito |
| Contracorrente inclinada | não coberta | `CHEXAL_DOMAIN` (61) |
| Horizontal com fluxo negativo | não coberta, pp. 264–265 | `CHEXAL_DOMAIN` (61) |
| Gás ↓ com líquido ↑ | não coberta | `CHEXAL_DOMAIN` (61) |

A catenária D vai de 0° (base, horizontal) a 76,47° no topo e não tem trecho vertical.
Logo, **nenhuma** célula do riser admite contracorrente com a correlação de 1992.

## 3. Auditoria da implementação existente

As equações 3-2 a 3-41 foram relidas no PDF e comparadas linha a linha com `component` e
`evaluate` em `scr/module_chexal_1992.f90`:

- B1 = min(0,8, A1);
- seleção Re = Reg se Reg > Ref ou Reg < 0;
- os dois ramos de C2 (ρf/ρg ≤ 18 e C5);
- C4 com D2 = 0,09144 m;
- C3 ascendente max(0,5, 2e^{−|Ref|/60000}) e descendente/CCFL (3-35–3-38);
- C9 com mínimo 0,7 para Reg < 0;
- 3-4b como máximo para Reg < 0;
- Fr para Reg < 0 (3-3b) e Vgj descendente com (Fr−1) (3-22);
- conversão de ângulo θ_C = 90 − θ·180/π;
- fluxos transformados para eixo descendente.

**Nenhuma divergência.** Itens conferidos: sinais (Reynolds assinados, gravidade pelo
ângulo), limites α → 0/1 (fase ausente só por inventário nulo, nunca por fluxo nulo),
contracorrente vertical com ramo, recuperação de primitivas com a identidade exata do
momento, e EOS (ρl, ρg entram em SI pelo mesmo `pressure_fun_uj` em BAROTROPIC e
ISOTHERMAL). Os status 60–66 se propagam com códigos compostos sem conversão em `NO_ROOT`.

Limitação já documentada que reaparece na campanha (§7): a correlação é por partes
(min/max). As derivadas Richardson **recusam** atravessar uma quina (`CLOSURE_DERIVATIVE`)
em vez de suavizar.

## 4. Comparação com referências

**Referência externa (artigo e relatório EPRI).** O artigo de 1992 não traz exemplos
numéricos. O EPRI NP-3989-SR Rev. 1 (1986), Apêndice A, Tabela A-1, imprime oito problemas
de vapor–água vertical com C0, Vgj e α. Na p. 260, o artigo de 1992 afirma que a forma
**ascendente** vapor–água "remains the same" desde 1986; a forma descendente foi revista.

Uma transcrição Python independente das equações de 1992
([np3989_sample_check.py](../validation/diagnostics/chexal_riser_20261121/np3989_sample_check.py)),
com o L de vapor–água (3-5/3-6), foi alimentada com estas propriedades:
- ν a partir dos Reynolds impressos na tabela;
- σ pela fórmula IAPWS 1994;
- ρ de tabela de vapor, com banda ±0,5 % em ρ, ±2 % em σ e P_crit 3200–3208 psia.

| Caso | P psia, D ft | Papel | C0 tabela / calc. | Vgj ft/s tabela / calc. | α tabela / calc. | Status |
|---|---|---|---|---|---|---|
| 1 | 14,7; 0,05 | referência | 1,2037 / 1,2038 | 0,5979 / 0,5985 | 0,5361 / 0,5361 | PASS (na banda) |
| 4 | 14,7; 1,0 | referência | 1,1922 / 1,1922 | 0,9054 / 0,9064 | 0,5323 / 0,5322 | PASS (na banda) |
| 2 | 1000; 0,05 | referência | 1,1116 / 1,1107 | 0,1410 / 0,1373 | 0,5947 / 0,5953 | FAIL (não reproduzido) |
| 3 | 1000; 1,0 | referência | 1,1119 / 1,1109 | 0,2234 / 0,2262 | 0,5914 / 0,5921 | FAIL (não reproduzido) |
| 8 | 14,7; 0,05, descendente | informativo | 1,3036 / 1,3036 | 4,6590 / 4,6587 | 0,6714 / 0,6713 | INFO |
| 5–7 | descendente | informativo | difere | difere | difere | INFO (forma revista em 1992) |

Os casos 2–3 a 1000 psia não são reproduzidos (Vgj −2,6 % e +1,3 %, C0 −0,08 %). Os
desvios de Vgj têm **sinais opostos** entre os dois diâmetros, o que exclui um erro de
escala de propriedades. Na p. 3-3 do NP-3989, C2 de 1986 não tem o ramo ρf/ρg ≤ 18, mas a
1000 psia ρf/ρg ≈ 20,6 e os dois ramos dão C2 = 1. **A causa não foi identificada.**
Esses casos ficam fora do regime do PRIMS-F (ar–água a ~2–5 bar, ρl/ρg ≈ 200–800),
que é o dos casos 1 e 4, reproduzidos a 1e-4–1e-3.

**Núcleo Fortran contra a transcrição independente**
([core_crosscheck.py](../validation/diagnostics/chexal_riser_20261121/core_crosscheck.py)):
945 estados ar–água com α 0,02–0,97, ρg 1,2–4,2, θ 0/5/30/51,6/76,47/89/90° e fluxos
ascendentes, descendentes, nulos e contracorrentes.

- 750 aceitos pelos dois: diferença máxima de C0 2,4e-16 e de Vgj 6,9e-14 (relativa).
- 195 recusados pelos dois: contracorrente inclinada/horizontal, horizontal negativo, gás ↓ com líquido ↑.
- Contracorrente vertical: o núcleo exige ramo (AUTO → 62); com HIGH explícito, concorda com a Eq. 3-40 transcrita.
- **Zero** divergências de status.

**Testes do núcleo já existentes**, contra referências mpmath de 60 dígitos, com GNU
Fortran 13.3.0 e `-fcheck=all -ffpe-trap`:

- `test_chexal_1992_core`: 147 verificações (vertical ↑/↓, horizontal, CCFL e ramos, j = 0, fases ausentes, unidades);
- `test_chexal_velocity_reference`: 21;
- `test_chexal_derivative_tail`: 87;
- `test_chexal_integration`: 612;
- `test_chexal_junction`: 120.

Todos PASS.

## 5. Primitivas, fluxos e Jacobianos

- **Primitivas e fluxo pontual:** `test_chexal_integration` cobre ida e volta U ↔ primitivas,
  derivadas de fluxo contra perturbações independentes de U, hiperbolicidade sem reparo
  espectral, contracorrente vertical com ramo e recusa de contracorrente inclinada.
- **Jacobiano global montado, nova evidência**
  ([global_jacobian_check.f90](../validation/diagnostics/chexal_riser_20261121/global_jacobian_check.f90)):
  a configuração completa 51+51 da catenária D foi reconstruída a partir de estados
  **aceitos** pelo solver Chexal (`fields.csv`). O resíduo completo foi montado
  (oleoduto estratificado, junção MATCHED, riser Chexal, Roe, LAGGED), e cada uma das
  306 colunas dos blocos foi comparada com diferenças centradas do resíduo completo em
  h = 1e-4, 1e-5 e 1e-6 (vezes a escala). O erro cai com h², como esperado:

| Estado | Fechamento | Máx. erro (melhor h): oleoduto / junção / riser | Colunas > 3e-4 |
|---|---|---|---:|
| OSC t = 60 s | Chexal | 4,7e-8 / 4,1e-8 / 1,5e-9 | 0 |
| OSC t = 135 s | Chexal | 9,5e-9 / 1,1e-8 / 2,8e-9 | 0 |
| Stable t = 600 s | Chexal | 2,4e-9 / 5,7e-9 / 1,5e-9 | 0 |
| SS1 t = 30 s | Chexal | 3,4e-9 / 4,3e-8 / 3,2e-9 | 0 |
| SS3 t = 30 s | Chexal | 6,4e-9 / 2,3e-8 / 2,0e-8 | 0 |
| OSC t = 60 s (controle) | Bhagwat | 1,7e-8 / 1,4e-8 / 2,2e-9 | 0 |

- **Uso efetivo do Chexal em todo o caminho**
  ([fields_closure_check.py](../validation/diagnostics/chexal_riser_20261121/fields_closure_check.py)):
  em cada instante gravado e em cada célula do riser, o resíduo Chexal foi recalculado pela
  transcrição independente a partir de (α, ρl, ρg, ul, ug, θ). Execuções Chexal (8): máximo
  relativo ≤ 8,3e-14 em todos os estados (8211 no Stable). Execuções Bhagwat (8):
  6,2e-2 a 6,6e-1, o que mostra que o teste discrimina o fechamento. Esse mesmo teste
  detectou um erro do meu lançador (§8) antes que qualquer resultado fosse usado.

## 6. Junção estratificado–Chexal (MATCHED)

A junção continua uma face interna com **um** fluxo mássico por fase compartilhado.

- Em todos os estados da tabela acima: `|Fm − Fp|` das fases = 0 exato e status 0.
- O α* do lado do riser, obtido em θ = 0 com Chexal horizontal, satisfaz a relação Chexal
  com resíduo relativo 2e-16 (OSC t = 60 s: Gl* = 0,4714, Gg* = 1,11e-3 adimensionais,
  α*_oleoduto = 0,551, α*_riser = 0,364).
- `test_chexal_junction` (120) cobre conservação, telescopagem de inventário com larguras
  desiguais, traço horizontal com gás reverso recusado por todos os ramos (sem clipping)
  e contornos.
- E1/E2 foram estudados com Bhagwat e **não** foram requalificados: ficam recusados com
  Chexal (§1).

**Limite encontrado (SS1).** No estado aceito em 52,94 s, a junção está bloqueada:
Gg* = 5,8e-13 adimensional e α*_riser = 3,6e-10. Toda proposta de Newton exige Gg* < 0
e é recusada com 3061 (traço do riser θ = 0 fora do domínio, p. 264–265). É o mesmo
bloqueio de inversão de gás estudado com Bhagwat (1012/1021, E2), agora alcançado mais
cedo (52,9 contra 77,8 s).

## 7. Campanha Stable, OSC, SS1 e SS3

Protocolo:

- decks primários com mudança apenas de `closure_riser` (5 ou 3);
- `tol_res = 1e-10`; N51+51; Δt_max 0,5 s; 2400 s; Roe, LAGGED, MATCHED;
- BAROTROPIC e ISOTHERMAL;
- limites de custo pré-declarados (4 h por execução, estagnação de 1800 s);
- binário congelado `d342750…`.

Tabela completa em [evidence/P1_CAMPAIGN.md](../validation/diagnostics/chexal_riser_20261121/evidence/P1_CAMPAIGN.md).

| Caso / EOS | Chexal (C5) | t s | Bhagwat (C3) | t s |
|---|---|---:|---|---:|
| Stable BARO | COMPLETED | 2400 | COMPLETED | 2400 |
| Stable ISO | COMPLETED | 2400 | COMPLETED | 2400 |
| OSC BARO | BLOCKED_CONSTITUTIVE_DOMAIN (célula 103, 76,1°, j_l = 1,0e-9 m/s) | 145,70 | COMPLETED, mas **decai** ao estacionário | 2400 |
| OSC ISO | BLOCKED_CONSTITUTIVE_DOMAIN (célula 101, 75,5°) | 111,13 | FAIL MAX_NEWTON | 382,98 |
| SS1 BARO | BLOCKED_JUNCTION_TRACE_DOMAIN (3061) | 52,94 | FAIL 3012 (junção) | 77,81 |
| SS1 ISO | BLOCKED_JUNCTION_TRACE_DOMAIN | 60,56 | FAIL MAX_NEWTON | 50,24 |
| SS3 BARO | BLOCKED_CONSTITUTIVE_DOMAIN (célula 85, 68,7°, j_l = 9,2e-10 m/s) | 41,55 | FAIL 1017 (primitiva Bhagwat sem raiz) | 97,97 |
| SS3 ISO | BLOCKED_CONSTITUTIVE_DOMAIN (célula 80, 65,2°) | 50,73 | FAIL 3012 | 72,63 |

O mecanismo dos bloqueios constitutivos foi verificado de forma independente
([domain_boundary_probe.py](../validation/diagnostics/chexal_riser_20261121/domain_boundary_probe.py)):

- O status 1071 = 1000 + 10 + 61 é a recuperação de primitivas com `CHEXAL_DOMAIN`.
- No estado aceito, o j_l recuperado decresce linearmente para o topo (OSC BARO: 8,2;
  5,8; 3,7; 1,7 mm/s nas células 99–102) e é ~1e-9 m/s na célula bloqueada.
- A busca linear recusa até λ = 7,6e-6. O passo seguinte exige líquido descendente com
  gás ascendente a 65–76° da horizontal: contracorrente inclinada.
- No OSC Chexal, a amplitude de p_in **cresce** (±230 Pa nos primeiros 20 s → ±14 kPa em
  146 s) e j_l na saída cai de 0,54 para 0,014 m/s. É uma instabilidade em crescimento
  que leva ao refluxo do líquido no topo.
- O Bhagwat, com a mesma discretização, **decai**: a oscilação relativa no último terço é
  0 em double. Com Δt = 0,5 s, a comparação OSC é dominada pelo amortecimento do Euler
  implícito e não diz nada sobre o período.

**Refino (OSC e SS3 Chexal, BAROTROPIC):**

| Caso | N (por seção) | Δt_max s | Evento | t s | Célula de bloqueio (θ, j_l) |
|---|---:|---:|---|---:|---|
| OSC | 51 | 0,5 | BLOCKED_CONSTITUTIVE_DOMAIN | 145,70 | topo, 76,08°, 1,0e-9 m/s |
| OSC | 51 | 0,25 | BLOCKED_CONSTITUTIVE_DOMAIN | 110,61 | topo, 76,08°, 4,3e-10 |
| OSC | 51 | 0,125 | BLOCKED_CONSTITUTIVE_DOMAIN | 93,10 | topo, 76,08°, 1,2e-9 |
| OSC | 101 | 0,5 | BLOCKED_CONSTITUTIVE_DOMAIN | 163,53 | topo, 76,28°, 8,1e-10 |
| OSC | 201 | 0,5 | BLOCKED_CONSTITUTIVE_DOMAIN | 181,21 | topo, 76,38°, 1,3e-9 |
| OSC | 101 | 0,25 | BLOCKED_CONSTITUTIVE_DOMAIN | 111,16 | topo, 76,28°, 1,1e-9 |
| SS3 | 51 | 0,5 | BLOCKED_CONSTITUTIVE_DOMAIN | 41,55 | 68,72°, 9,2e-10 |
| SS3 | 51 | 0,25 | BLOCKED_CONSTITUTIVE_DOMAIN | 41,95 | 67,42°, 3,0e-9 |
| SS3 | 51 | 0,125 | BLOCKED_CONSTITUTIVE_DOMAIN | 42,18 | 65,96°, 5,6e-9 |
| SS3 | 101 | 0,25 | BLOCKED_CONSTITUTIVE_DOMAIN | 44,65 | 71,23°, 2,5e-10 |
| SS3 | 101 | 0,5 | FAIL 9323: derivada Chexal em quina no caminho de Roe (estado aceito) | 20,45 | — |
| SS3 | 201 | 0,5 | FAIL 9323, idem | 9,95 | — |

Leitura do refino:

- **OSC.** O evento constitutivo aparece em todas as seis discretizações, sempre na
  **célula do topo** (76,1–76,4°) e com j_l → 0. O instante **não** converge:
  - em tempo (N51), cai de 145,7 para 110,6 e 93,1 s, porque um Δt menor amortece menos
    e o crescimento é mais rápido;
  - em espaço (Δt 0,5), sobe de 145,7 para 163,5 e 181,2 s.

  A taxa de crescimento da instabilidade depende fortemente da discretização. A
  conclusão robusta é qualitativa: o OSC com Chexal cresce até exigir refluxo de
  líquido no topo. O instante desse evento não está convergido, e o período e a
  amplitude não são definidos (nenhum ciclo completo).
- **SS3.** Com N51, o instante converge em tempo (41,55 → 41,95 → 42,18 s; diferenças
  0,40 e 0,23 s, ordem aparente ≈ 0,8), e a célula de bloqueio desce de 68,7° para 66,0°.
  Com N101 e Δt 0,25, o bloqueio vem em 44,65 s a 71,2°.
- **Falha 9323 (N101/N201 com Δt 0,5).** Em um estado **aceito**, ao preparar a dissipação
  LAGGED: 9323 = 8000 + 1200 + 100 + 10 + 13 (`CLOSURE_DERIVATIVE`). Um nó de Gauss do
  caminho de Roe numa face do riser caiu sobre uma quina da correlação (min/max), e a
  derivada Chexal recusou atravessá-la, como especificado. É limitação do operador de
  dissipação com uma correlação por partes, não domínio físico. Não foi suavizada.
- **Balanços:** máx. |balanço|/M₀ ≤ 2,5e-11 nas 10 execuções de refino.

Tabelas: [P1_REFINE.md](../validation/diagnostics/chexal_riser_20261121/evidence/P1_REFINE.md),
[DOMAIN_PROBE.json](../validation/diagnostics/chexal_riser_20261121/evidence/DOMAIN_PROBE.json).

**Stable (preservação do estacionário).** Chexal e Bhagwat retornam a um estado
estacionário. A diferença final em p_in em relação ao perfil inicial é +26,4 Pa (5e-5),
contra +22,3 Pa com Bhagwat: é a distância entre o equilíbrio FV e o perfil da EDO, da
mesma ordem nos dois. Os estados estacionários **diferem** por causa do fechamento:

| Stable BARO | p_in Pa | M_l kg | M_g kg | M_g riser kg |
|---|---:|---:|---:|---:|
| Chexal | 490 432 | 80,40 | 0,2960 | 0,0244 |
| Bhagwat | 496 400 | 82,03 | 0,2931 | 0,0182 |

O Chexal prevê 34 % mais gás no riser e p_in 1,2 % menor.

**Eventos e inventários.** E_g(0,01), o primeiro instante em que o fluxo de gás na base cai
abaixo de 1 % do inicial, é um limiar de **fluxo**, não desaparecimento de fase:
- SS1 Chexal 52,45 s (BARO) e 59,99 s (ISO);
- SS1 Bhagwat 77,81 s e 48,95 s.

Nenhuma execução Chexal chega à fase de produção ou blowdown do SS1. Inventários e faixas
de pressão por execução estão em `P1_CAMPAIGN.json`.

Custo: CPU por execução de 87 a 547 s. O Chexal é ~2× mais caro que o Bhagwat no Stable
(536 contra 276 s), por causa das derivadas por Richardson e da enumeração de raízes.

## 8. Regressões e balanços

- **Build:** release (`-O3`, `make`) e debug (`make check_all`, `-fcheck=all -ffpe-trap`) sem erro.
- **Fortran (`scripts/run_unit_tests.py`):** 26/26 programas e 5947 asserções. Esses testes
  fixam o fechamento explicitamente, logo o default novo não os altera em silêncio.
- **Python:** `unittest discover` 123 OK.
- **CLI** (os seis scripts do CI, com o binário de verificação):
  - `test_pipe_regime_cli`: 23 PASS. Atualizado: o caso ACTIVE_SET agora seleciona
    Bhagwat explicitamente e compara com um par Bhagwat, e ACTIVE_SET + Chexal deve ser recusado.
  - `test_accepted_telemetry_cli`, `test_feed_contract_cli` (11) e `test_initial_inventory_cli`: PASS.
  - `test_boundary_telemetry_cli`: PASS com `closure_riser=3` fixado. Esse teste é de
    ledger de contorno; com Chexal em `-O0 -fcheck=all`, um caso de 9 passos levou 197 s
    e estourava o limite de 90 s do script.
  - **Novo** `tests/test_chexal_default_cli.py` (11 execuções, incluído no CI):
    1. deck sem `closure_riser` roda Chexal, e seus estados de riser satisfazem a relação a 3,5e-15;
    2. `closure_riser=3` continua aceito;
    3. E1, E2 e ACTIVE_SET são recusados com IDs 5, 6 e 7 antes de qualquer saída;
    4. os decks primários atribuem 5 uma única vez fora de comentários, e os arquivados, 3.
- **Regressão do Bhagwat contra o binário congelado `364b808`** (mesmo `scr/` de `d2c8101`):
  os decks C3 da campanha rodados com o binário anterior dão `monitor`, `fields`,
  `conservation`, `boundary_conditions`, `rejected_proposals` e `last_rejected_state`
  **idênticos byte a byte**; `steps.csv` e `run_status` diferem só em CPU e tempo de relógio.
  Casos: SS1 BARO e SS3 ISO, até a mesma falha terminal (3012); Stable BARO até 300 s, com as 609 primeiras linhas do monitor idênticas. Só o último passo difere, encurtado para terminar exatamente em 300 s ([REGRESS_BHAGWAT.json](../validation/diagnostics/chexal_riser_20261121/evidence/REGRESS_BHAGWAT.json)).
- **Mudança intencional:** um deck que não declara `closure_riser` passa de Bendiksen para Chexal.
- **Balanços:** seis por execução (duas fases × sistema, oleoduto e riser), pelos fluxos
  aceitos com `audit_ledger`. Máximo |balanço|/M₀ por execução de 2,7e-12 a 3,0e-7
  (critério 1e-6), aritmética do ledger PASS nas 16, descasamento de integral na junção
  = 0 kg (exato) em todas.
- **Erro de protocolo detectado e corrigido.** A primeira tentativa da campanha (09h35)
  editava a primeira ocorrência de `closure_riser=`, que estava no comentário de
  proveniência que eu mesmo havia escrito. As execuções "C3" rodaram Chexal. A verificação
  de campos acusou o erro em minutos. As saídas foram para `p1_INVALID_deck_override_bug/`
  com `INVALID.txt` e não são usadas. Correções: o comentário não contém mais a chave, o
  editor ignora comentários e verifica o valor escrito, e o mesmo foi aplicado ao `assign`
  de `scripts/run_validation.py`.

## 9. Limitações de contracorrente inclinada

O resultado central da campanha é que, com Chexal em todo o riser, **OSC e SS3 exigem
contracorrente inclinada no trecho superior da catenária**, e SS1 exige inversão de gás
no traço horizontal da junção. Nenhum dos dois estados está no domínio de 1992:

- contracorrente só vertical (p. 262);
- horizontal só cocorrente (pp. 264–265);
- a p. 270 registra ainda que, mesmo na vertical, não havia dados fora da linha CCFL.

O código para com status explícito. Não houve:
- extensão para inclinado;
- troca automática para Bhagwat ou Bendiksen;
- suavização de quinas;
- clipping de velocidade;
- reservatório ou fonte artificial na junção.

Uma extensão experimental fundamentada exigiria:

1. dados de C0/Vgj em contracorrente inclinada ar–água (por exemplo, a literatura de
   *liquid loading* e ZNLF em tubos inclinados) na faixa α ≈ 0,5–0,6, θ ≈ 65–76°, D ≈ 50 mm;
2. uma forma de CCFL inclinada e o ramo físico correspondente, com a regra de seleção de
   ramo justificada (p. 268);
3. testes locais contra esses dados antes de qualquer transiente.

A transcrição independente (`chexal92_ref.py`) e as ferramentas de fronteira de domínio
permitem medir quão longe da fronteira cada estado está.

## 10. Matriz de status

| Item | Status | Evidência |
|---|---|---|
| Chexal como default de todo o riser (dev) | PASS | §1, CLI novo, verificação de campos |
| Bhagwat preservado como opção legada | PASS | regressão byte a byte, §8 |
| Equações 1992 contra PDF | PASS | §3 |
| Referência externa NP-3989 (baixa pressão, regime PRIMS-F) | PASS | casos 1 e 4 |
| Referência externa NP-3989 (1000 psia) | FAIL (não reproduzido, causa não identificada) | casos 2 e 3 |
| Núcleo Fortran contra transcrição independente | PASS | 945 estados |
| Primitivas, fluxos e Jacobiano pontual | PASS | testes existentes |
| Jacobiano global acoplado (FD, 3 incrementos) | PASS | 6 estados × 306 colunas |
| Junção MATCHED estratificado–Chexal | PASS (estados admissíveis) | §6 |
| E1/E2 com Chexal | NOT_VERIFIED (recusados até requalificação) | §1 |
| ACTIVE_SET com Chexal | NOT_VERIFIED (recusado) | §1 |
| Stable (BARO, ISO) | PASS (estacionário preservado) | §7 |
| OSC Chexal | BLOCKED_CONSTITUTIVE_DOMAIN | §7, refino |
| SS3 Chexal | BLOCKED_CONSTITUTIVE_DOMAIN | §7, refino |
| SS1 Chexal | BLOCKED_JUNCTION_TRACE_DOMAIN | §6–7 |
| Refino SS3 N101/N201 (Δt 0,5) | FAIL: derivada Chexal em quina no caminho de Roe (9323) | §7 |
| Convergência OSC (espaço, tempo, conjunta) com Chexal | NOT_VERIFIED: sem ciclos; instante do bloqueio não converge (93–181 s) | §7 |
| Balanços (6 por execução) | PASS | §8 |
| Validação experimental (Baliño et al. 2010: SS1 85,6 s; SS3 76 s; OSC 19,3–33,8 s) | NOT_VERIFIED: nenhuma execução Chexal tem ciclos completos admissíveis | §7 |
| Promoção à produção | BLOCKED: lacuna de contracorrente inclinada | §9 |

## 11. Commits, hashes, comandos e custos

Commits locais na branch `implementation/chexal-riser-20261121`, sem push e sem merge:

| SHA | Conteúdo |
|---|---|
| `461b2b5` | feat: Chexal como default de todo o riser; decks arquivados; recusas ACTIVE_SET/E1/E2; `assign` só fora de comentários |
| `8cce58a` | test: CLI novo do contrato Chexal; casos só-Bhagwat adaptados; CI |
| `b74a9a3` | diag: transcrição independente, NP-3989, cruzamento, Jacobiano global, sondas, campanha |
| `e3dfd3f` | evidence: campanha, refino, cruzamentos, Jacobianos, regressão |
| (este) | docs: este relatório, handoff, nota de integração |

Custos medidos (CPU de uma thread):
- campanha principal: 87–547 s por execução, 16 execuções;
- refino: 10 execuções, até 181 s simulados;
- `check_all` com Chexal: ~22 s por passo num caso de 4+4 células.

Binário da campanha: `prims_f_chexal_default`, SHA-256
`d342750229f4c250d87888e0ed0cc37ad29a16293d41faad9004e0509b400ad2`, GNU Fortran 13.3.0,
com o diff exato da fonte em `bin/source_diff_vs_d2c8101.patch`. Protocolo e decks com
hashes: `p1/PROTOCOL_*.json` e `p1_refine/PROTOCOL_*.json`.

```sh
make release && make check_all
PYTHONPATH=scripts python3 scripts/run_unit_tests.py
python3 -m unittest discover -s tests -p 'test_*.py'
python3 tests/test_chexal_default_cli.py --binary bin/prims_f_check_all --output tmp/ci_chexal_default
D=validation/diagnostics/chexal_riser_20261121
python3 $D/run_campaign.py BIN RESULTS/p1 --workers 6
python3 $D/run_campaign.py BIN RESULTS/p1_refine --workers 3 --eos BAROTROPIC --closures 5 \
    --cases D_OSC D_SS3 --grid 51:0.25 51:0.125 101:0.5 201:0.5 101:0.25
python3 $D/summarize_campaign.py OUT.json OUT.md RESULTS/p1/runs/*/
python3 $D/fields_closure_check.py OUT.json RESULTS/p1/runs/*/
python3 $D/np3989_sample_check.py OUT.json OUT.md
python3 $D/core_crosscheck.py tmp/cx/core_crosscheck OUT.json
$D/run_global_jacobian.sh tmp/cx/global_jacobian_check RESULTS/p1/runs EVIDENCE
```
