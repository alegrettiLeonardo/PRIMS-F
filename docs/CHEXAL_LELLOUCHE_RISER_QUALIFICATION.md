# Chexal–Lellouche (1992) como fechamento de todo o riser: qualificação

Data: 21/11/2026. Worktree `PRIMS-F-chexal-riser-20261121`, branch
`implementation/chexal-riser-20261121`, base `d2c8101` (E1/E2 experimentais
preservados, desligados). Resultados brutos: `~/Documentos/PRIMS-F-chexal-riser-results-20261121`.
Ferramentas e evidências: [validation/diagnostics/chexal_riser_20261121](../validation/diagnostics/chexal_riser_20261121/README.md).
Sem push, merge ou alteração da `main`.

## Resumo

> **Rodada 2 (Parte II, ao final):**
> - NP-3989 auditado nas versões 1986 e 1992: a revisão descendente está demonstrada; as linhas a 1000 psia continuam FAIL, com causa não identificada.
> - Ângulos e velocidades certificados.
> - Contrato do fechamento solicitado = executado (incidente de rotulagem coberto por teste negativo).
> - Jacobianos nas fronteiras constitutivas verificados.
> - Validação ar–água NOT_VERIFIED (sem pontos primários).
> - Contracorrente inclinada BLOCKED.
> - Campanha reexecutada com contrato e idêntica byte a byte.
> - Stable estável também com Δt 0,125.

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

---

# Parte II — Qualificação P0/P1/P2 (21/11/2026, segunda rodada)

Base: HEAD `db97d3b` confirmada, worktree limpa e `main` em `f3717dd`. Os resultados da
Parte I foram tratados como evidência histórica. Esta rodada reexecutou a campanha
com o binário novo (contrato de fechamento, §II.3) e comparou as execuções com as
anteriores (§II.8). A Parte II prevalece sobre a Parte I onde as duas tratam do mesmo item.

## II.1 Auditoria EPRI a 1000 psia e diferenças 1986/1992 (P0.1)

**Equações das duas versões.** Fontes: NP-3989 pp. 3-1 a 3-4 (PDF 25–28) e Chexal et al.
(1992) pp. 261–268.

| Termo | 1986 | 1992 | Efeito na Tabela A-1 |
|---|---|---|---|
| C0, L(α,p), K0, r, B1, Re | idênticos | idênticos | nenhum |
| C2 | 1 ou 1/(1−e^{−C5/(1−C5)}), C5 = √(150 ρg/ρf) | idem + ramo 0,4757 [ln(ρf/ρg)]^0,7 se ρf/ρg ≤ 18 | nenhum (ρf/ρg ≥ 20,6 em todas as linhas) |
| C3 ascendente, C4, D2 | idênticos | idênticos | nenhum |
| C3 descendente | C3 = C10 | C3 = 2(C10/2)^B2 | casos 5–8 |
| expoente de (1−α) com Reg<0 | K1 = min(0,65; 0,5 e^{|Reg|/4000}) | C9 = min(0,7; (1−α)^0,65) | casos 5–8 |
| C0 com Reg<0 | só 3-4a | máx(3-4a; Vgj⁰(1−α)^0,2/(|jf|+|jg|)) (3-4b) | casos 5–8 |

**Unidades e propriedades.** A tabela usa P em psia (absoluta), D em ft, j e Vgj em ft/s;
o modelo é numericamente idêntico em SI (p. 3-4). As propriedades vapor–água não são
impressas e foram declaradas: tabelas clássicas de vapor saturado (14,7 psia:
vf 0,01672, vg 26,80 ft³/lb; 1000 psia: vf 0,02159, vg 0,4456 ft³/lb), σ pela IAPWS 1994
em Tsat, viscosidades a partir dos Reynolds **impressos** e P_crit = 3208,2 psia. Nenhuma
propriedade ar–água foi usada.

**Teste por linha** (`tests/test_chexal_np3989.f90`): o núcleo Fortran ganhou o L de
vapor–água (Eq. 3-5/3-6) como opção explícita de verificação. O caminho ar–água continua
bit a bit, conforme o cruzamento de 945 estados repetido. Em todas as linhas, Fortran =
transcrição independente de 1992 com erro ≤ 1e-11.

| Caso | Versão aplicável | Referência C0/Vgj/α | Fortran (1992) | Erro máx. | Causa | Status |
|---|---|---|---|---:|---|---|
| 1 | ambas (ascendente) | 1,2037/0,5979/0,5361 | 1,2038/0,5984/0,5361 | 8,6e-4 | reproduzido na banda de propriedades | PASS |
| 4 | ambas | 1,1922/0,9054/0,5323 | 1,1922/0,9063/0,5322 | 9,8e-4 | idem | PASS |
| 8 | ambas | 1,3036/4,6590/0,6714 | 1,3037/4,6584/0,6713 | 1,5e-4 | idem | PASS |
| 5 | 1986 | 1,0400/5,1172/0,9538 | 1,5869/11,391/0,8056 | — | revisão descendente de 1992; as equações de **1986** reproduzem (1,0399/5,1162/0,9539) | PASS (diferença de versão demonstrada) |
| 3 | ambas | 1,1119/0,2234/0,5914 | 1,1109/0,2263/0,5921 | 1,3e-2 | a linha impressa viola α(C0 j+Vgj)=jg em −4,2e-4, limite de arredondamento 1,3e-4: referência defeituosa | FAIL (referência inconsistente) |
| 2 | ambas | 1,1116/0,1410/0,5947 | 1,1107/0,1373/0,5953 | 2,6e-2 | não explicada | FAIL |
| 6 | ambas | 1,0498/4,8444/0,9171 | 1,0921/4,6374/0,8515 (1986: 1,0514/4,8318/0,9142) | 7,2e-2 | não explicada | FAIL |
| 7 | ambas | 1,2948/1,8775/0,5700 | 1,2934/1,8193/0,5688 (1986: 1,2927/1,8363/0,5697) | 3,1e-2 | não explicada | FAIL |

Hipóteses testadas para a falha a 1000 psia ([np3989_audit.py](../validation/diagnostics/chexal_riser_20261121/np3989_audit.py),
[np3989_property_fit.py](../validation/diagnostics/chexal_riser_20261121/np3989_property_fit.py)):

- **Revisão de versão:** rejeitada. A forma ascendente é idêntica, e os casos 6–7 também
  falham com as equações de 1986.
- **Reynolds efetivo diferente em B1/C3:** rejeitada (erro mínimo ainda 5,7e-3).
- **Um único conjunto de propriedades:** o ajuste é degenerado; só ρg/ρf e σ(ρf−ρg)/ρf²
  são identificáveis. Os casos 2, 6 e 7 são reproduzidos a ≤ 1,1e-3 com ρg/ρf = 0,0468
  (tabela: 0,0485) e σ(ρf−ρg)/ρf² +10 %. Esse par não corresponde a um estado saturado
  coerente (a razão indica ~970 psia, a σ ~273 °C), e o caso 3, fora do ajuste, piora.

A causa continua **não identificada** e os casos 2, 6 e 7 ficam FAIL, como exige o gate.

Relevância: a 1000 psia, C5 ≥ 1 e C2 = 1. Nos estados ar–água do PRIMS-F, ρg/ρl ≈
0,002–0,006 e C5 ≈ 0,5–0,95 < 1, o mesmo ramo das linhas de 14,7 psia, que são reproduzidas.
Nenhum coeficiente foi ajustado.

**Gate P0.1:** parcialmente atendido. Os valores aplicáveis a baixa pressão e a diferença
de versão descendente estão demonstrados; as linhas a 1000 psia continuam **FAIL**.

## II.2 Certificado de ângulos e velocidades (P0.2)

Convenções:
- **PRIMS-F:** θ em rad a partir da horizontal, positivo quando o eixo sobe (riser 0 a 76,47°,
  oleoduto −β), fluxos assinados ao longo do eixo.
- **Chexal:** θ_C = 90 − |θ| graus a partir da vertical e fluxos no sentido físico ascendente,
  f = sign(θ)·j.
- A correlação usa o ângulo só por Fr (3-3); o seno da gravidade fica no termo-fonte, coberto
  por `test_fv_operators`.

`tests/test_chexal_angles.f90`: **468 verificações, PASS**, em 0, ±15, ±45, ±65, ±76 e ±90°
e 12 pares de fluxo (cocorrente ↑/↓, contracorrente nos dois sentidos, nulos, ±1e-12).

- O status de domínio coincide com o previsto pelas pp. 262–265: 81 estados admitidos, 46
  recusas de domínio, 5 contracorrentes verticais com `BRANCH_REQUIRED`. Não há tolerância de
  sinal: 1e-12 m/s conta como sinal.
- A reversão de eixo (θ, j) → (−θ, −j) dá C0 idêntico e Vgj com sinal oposto, **bit a bit**.
- As velocidades recuperadas de (α, ρl j_l + ρg j_g) têm o sinal estrito, com erro ≤ 1e-12 na
  escala do momento.
- C0 e Vgj são contínuos em θ = 0. Há uma quina em Fr (inclinação 1/90 contra 1/10), sem salto.
- Com líquido nulo descendente, a raiz é o limite de gás puro α = 1 (C0(1) = 1, Vgj(1) = 0, pp. 259–260).
- Na **vertical**, o par (α, M) também admite uma solução contracorrente. Sem dica, AUTO recusa
  com 62; com a dica do estado anterior, recupera o ramo cocorrente exato. Isso nunca ocorre na
  catenária (≤ 76,47°).
- Com |j_g| < 1e-4 na vertical, aparece `NUMERIC` (66): é a limitação documentada da busca de
  CCFL (`CHEXAL_REMAINING_LIMITS.md`).
- No topo da catenária (76,47°), a contracorrente é recusada **mesmo com ramo explícito**.

**Gate P0.2: PASS.**

## II.3 Contrato do fechamento efetivo (P0.3)

Causa do incidente de 21/11 (manhã): o gerador de decks editou um `closure_riser=` dentro de
um comentário, e nenhum artefato da execução registrava o fechamento executado de forma
comparável ao rótulo.

Implementação:

- **`main.f90`:** `closure_riser` passa a ter sentinela, o que distingue DEFAULT de valor explícito.
- **`effective_closures`** (procedimento do simulador): devolve o fechamento **despachado** a
  cada célula pelo mesmo `cell_closure` usado nos operadores (estado pontual, faces, fontes,
  Jacobianos e recuperação).
- **`closure_contract.txt`**, gravado antes do transiente, com:
  - `closure_{pipe,riser}_{requested,effective}`, nome do modelo e células despachadas por seção;
  - fechamento do traço do riser na junção e do contorno de saída;
  - `gas_eos`, `junction_model/flux/trace/riser_limit`, globalização, dissipação, `tol_res` e arquivo de entrada.
- Se alguma célula despachada divergir do fechamento resolvido, o solver para (`error stop`).
- O cabeçalho do checkpoint traz os IDs de fechamento, a EOS e as opções da junção. Os leitores
  de checkpoint só usam as linhas seguintes, então continuam compatíveis.
- **`scripts/verify_closure_contract.py`** cruza:
  - deck (só linhas não comentadas) e rótulo `_C<k>`;
  - contrato, `run_metadata` e checkpoint;
  - **dados**: os estados aceitos do riser devem satisfazer o Chexal (transcrição
    independente) se o efetivo for 5–7, e não satisfazê-lo (≤ 1 %) caso contrário.
  
  Grava `CLOSURE_CONTRACT_STATUS.json`. `run_campaign.py` executa o verificador após cada
  execução e registra `source_sha256` (hash de `scr/*.f90`) e o SHA-256 do binário.
- **`tests/test_closure_contract_cli.py`** (CI), 7 casos:
  - deck sem chave → DEFAULT/5, PASS;
  - Bhagwat explícito → PASS;
  - **comentário com `closure_riser=3`, corpo 5, rótulo C3** (o incidente) → FAIL detectado
    (`label_riser == contract_effective`);
  - o mesmo deck rotulado C5 → PASS;
  - `closure_riser=9` → recusado antes de qualquer saída;
  - E1 + Chexal → recusado;
  - contrato adulterado depois da execução → FAIL em 5 verificações (deck, rótulo,
    metadados, checkpoint e dados).

**Gate P0.3: PASS.** Um rótulo incorreto é detectado automaticamente, inclusive pelos dados.

## II.4 Jacobianos nas fronteiras constitutivas (P0.4)

[boundary_jacobian_probe.f90](../validation/diagnostics/chexal_riser_20261121/boundary_jacobian_probe.f90)
([log](../validation/diagnostics/chexal_riser_20261121/evidence/BOUNDARY_JACOBIAN.log)):
- coluna a coluna, A = ∂F/∂U do operador do riser contra diferenças independentes de F(U);
- estados **consistentes** com o Chexal (j_g resolve a relação de deriva);
- h = c·max(|U_m|, 1e-14), c de 1e-2 a 1e-7;
- diferenças centradas, para frente e para trás **só quando admissíveis**; nenhum estêncil
  atravessa domínio inválido;
- autovalores por `dgeev_3x3`.

| Família | Melhor erro por coluna (pior coluna) | Leitura |
|---|---|---|
| α = 1e-2, 1e-4, 1e-6 (45°) | 3e-8 / 4,5e-7 / **2,2e-6** (∂F/∂U₂) | patamar estável em duas décadas de h: degradação real e pequena da derivada Richardson perto de α = 0, abaixo de 3e-4 |
| α = 0,99, 0,9999 | ≤ 8,5e-8 | PASS |
| j_l = 1e-2 … 1e-9 a 76,08° (estado de bloqueio do OSC) e 1e-8 a 68,7° (SS3) | ≤ 8,9e-10 no lado admissível | o estêncil para trás entra em contracorrente inclinada e é recusado; o Jacobiano é a derivada do lado admissível |
| j_g → 0⁺ a θ = 0 (α = 1e-3, 1e-6) | 1,7e-7 / 4,4e-6 | idem α → 0 |
| Re_l = 83 178 (B1 = 0,8 e C3 = 0,5 coincidem) ± 1e-7 | ≤ 1,5e-10 | PASS |
| L = 1 (Eq. 3-7) a α_L(1 ∓ 1e-8) | 7e-10 (para frente, abaixo) / 3,8e-9 (para trás, acima); central 1,1e-2 | **quina detectada**: A é a derivada do pedaço do estado; a central cruza o salto de inclinação |
| L = 1 exato e a 1e-12 | status 53 = 30 + 10 + 13 (`CLOSURE_DERIVATIVE`) | ponto não diferenciável **recusado**, não suavizado |
| θ = 1e-9 rad e 76,47° | ≤ 1,2e-9 | PASS |

Todos os 24 espectros calculados são reais (hiperbólicos). A razão |λ|max/|λ|min chega a
5,6e3 com α → 0 (rigidez acústico–cinemática). O Jacobiano global montado (Parte I, §5)
completa a verificação: 306 colunas, ≤ 4,3e-8.

Descontinuidade descoberta no estudo local (§II.6): o termo |Re_l|^0,03 da Eq. 3-36 dá a C3
inclinação infinita em Re_l = 0. A **contracorrente vertical** tem, portanto, uma fronteira
não diferenciável em j_l = 0⁻ (C3: 2 → ~0,9 já com |Re_l| ≈ 60). Na catenária ela não é
atingida, porque é recusada antes.

**Gate P0.4: PASS.** O erro fica dentro de tolerância justificada, e os pontos não
diferenciáveis são recusados explicitamente.

## II.5 Validação ar–água (P1.1)

Não há, no disco, pontos experimentais primários de fração de vazio ar–água: só o resumo dos
autores. A busca de uma tabulação aberta não encontrou dados anexados (Ghajar & Bhagwat,
OSU, D = 12,7 mm). Portanto **nenhum erro, viés ou RMSE foi calculado**; não se reportam
números que não são nossos.

Cobertura de domínio pela base declarada na Tabela 5-1 (p. 280) e na Tabela 5-2 (p. 282)
([validation_domain_map.py](../validation/diagnostics/chexal_riser_20261121/validation_domain_map.py)):

- **Inclinado: Beggs (1972).** D 25,4–38,1 mm; 2,56–6,81 bar; G 60–5500 kg/(m²·s);
  x 0,0008–0,833; 621 pontos. Os autores reportam erro médio 0,003 e desvio padrão 0,038 em α.
- **Horizontal:** 12,7–127 mm. **Vertical ascendente:** 11–305 mm. **CCFL vertical:**
  Bharathan, 25 e 51 mm.

| Execução Chexal | Estados | Validado | Interpolado | Extrapolado | Não suportado |
|---|---:|---:|---:|---:|---:|
| Stable BARO/ISO | 8211 cada | 0 | 1,000 | 0 | 0 |
| OSC BARO/ISO | 561/459 | 0 | 1,000 | 0 | 0 |
| SS3 BARO/ISO | 204/255 | 0 | 1,000 | 0 | 0 |
| SS1 BARO/ISO | 255/255 | 0 | 0,992 | 0,008 (x < 4e-5) | 0 |

- Nenhum estado é validado diretamente: D = 52,5 mm fica fora dos dados inclinados.
- "Interpolado" quer dizer D dentro das famílias horizontal e vertical que delimitam o ângulo
  (interpolação Fr, Eq. 3-2) e p, G, x cobertos por pelo menos uma família. É uma noção fraca,
  porque combina conjuntos distintos.

**Gate P1.1: NOT_VERIFIED** (dados primários indisponíveis). Concordância algébrica não foi
tratada como validação.

## II.6 Contracorrente inclinada (P1.2)

Fontes locais conferidas:
- **Chexal 1992:** contracorrente só vertical (p. 262); a região fora da CCFL não pôde ser
  verificada com dados (p. 270); CCFL ar–água só em tubos verticais (Tabela 5-2).
- **Ozaki & Hibiki 2015:** feixes de varetas; não remove a restrição de orientação.

Barnea, Ben Yoseph e Taitel, e Ohnuki, Adachi e Murao **não estão disponíveis localmente**. A
limitação fica registrada; nenhuma equação dessas fontes foi transcrita de memória.

Solver constitutivo local nos estados **reais** de bloqueio (11 células de fronteira; α
0,46–0,63; θ 66–76°), com velocidades assinadas
([inclined_cc_local.py](../validation/diagnostics/chexal_riser_20261121/inclined_cc_local.py),
[tabela](../validation/diagnostics/chexal_riser_20261121/evidence/INCLINED_CC_LOCAL.md)):

| j_l (m/s) | Chexal no ângulo real | Chexal vertical HIGH (sensibilidade, inadmissível em θ) | Bendiksen (correlação cocorrente estendida) | Bhagwat |
|---|---|---|---|---|
| 0 (borda admissível) | 1,10–1,82 | 1,07–1,79 | 0,79–1,34 | 0,45–0,51 (só a borda) |
| −1e-3 a −3e-2 | recusado (p. 262) | 0,48–0,92 | 0,72–1,31 | indefinido |

Valores em u_g (m/s). Os candidatos diferem por fator de 2,4–3,5 na borda e 1,4–1,8 em
contracorrente. O Chexal vertical cai de forma quase descontínua ao cruzar j_l = 0 (§II.4).
Não há base experimental, nos dados disponíveis, para escolher um candidato.

**Gate P1.2: BLOCKED.** Nada foi integrado ao solver global. O caminho para uma extensão
fundamentada está na Parte I, §9: dados C0/Vgj ou CCFL inclinados ar–água com D ≈ 50 mm e
θ ≈ 65–76°, e a forma do ramo.

## II.7 Junção estratificado–Chexal (P1.3)

- **Continuidade mássica:** fluxo comum `|Fm − Fp| = 0` exato em todos os estados testados.
  Descasamento da integral na junção = 0 kg em todas as execuções (ledger).
- **Pressão e momento:** os traços MATCHED usam p* e os fluxos de fase comuns. A força pontual
  (convenção T2) é preservada e não é apagada (`test_fv_operators`).
- **Estados e ramos de face:** o traço do riser em θ = 0 satisfaz o Chexal horizontal (2e-16).
  Fluxo de gás nulo é aceito; fluxo negativo é recusado com 3061 (domínio, pp. 264–265), sem
  clipping (`test_chexal_junction`).
- **Distinção dos tipos de falha:** SS1 C5 é domínio constitutivo do traço (3061: proposta com
  G_g* < 0 a partir de G_g* = 5,8e-13). Isso difere de falha de Newton sem causa constitutiva
  (MAX_NEWTON/6000 com status de proposta 0) e de raiz inexistente (Bhagwat 1017/1012).
- **Jacobiano completo:** colunas da junção ≤ 4,3e-8 em 6 estados aceitos.
- **Invariância à partição homogênea:** `test_partition_invariance` foi estendido a Chexal nos
  dois lados: resíduo idêntico a 1e-12 com a junção nas posições 0, 2 e 4 e fluxo de fase comum
  (28 asserções).
- **Condição de transmissão:** a escolha T1 (conservativa, Borda–Carnot) × T2 (pressão
  contínua, atual) foi estudada só com Bhagwat. Com T2 não havia solução de ondas fracas com
  gás no traço. Esse resultado **não foi transferido** para o Chexal.

**Gate P1.3:**
- **PASS** para conservação, admissibilidade, Jacobiano e partição.
- **NOT_VERIFIED/BLOCKED** para a fundamentação física da transmissão: a decisão T1/T2 é
  pendente e muda a física. E1/E2 continuam recusados com Chexal.

## II.8 Replays Stable/OSC/SS1/SS3 (P2.1)

Configuração congelada: binário `prims_f_contract`, SHA-256 `9e54e79f…` (fonte = HEAD
`7ba4d79` em `scr/`, patch em `bin/source_diff_contract_vs_db97d3b.patch`), GNU Fortran 13.3.0,
LAPACK/BLAS do sistema e uma thread. São os mesmos 16 decks primários, com `tol_res = 1e-10`, e
estados estacionários recalculados pelo solver com o fechamento de cada execução.

- **Contrato:** 16/16 `CLOSURE_CONTRACT_STATUS = PASS` (8 efetivo 5, 8 efetivo 3).
  Deck, rótulo, contrato, metadados, checkpoint e dados concordam
  ([P2_CONTRACT_STATUS.txt](../validation/diagnostics/chexal_riser_20261121/evidence/P2_CONTRACT_STATUS.txt)).
- **Reprodução:** as 16 execuções são **idênticas byte a byte** à campanha da Parte I
  ([P2_VS_P1_BYTEWISE.json](../validation/diagnostics/chexal_riser_20261121/evidence/P2_VS_P1_BYTEWISE.json)).
  Isso inclui monitor, campos, conservação, contornos, propostas rejeitadas, estado rejeitado,
  `steps` sem as colunas de tempo de relógio e `run_status`.

  Duas consequências: a opção de fluido e o contrato não alteram a física; e as conclusões da
  Parte I (§7), cujos rótulos agora estão verificados, valem integralmente. Tempos alcançados,
  passos, rejeições, resíduos, motivos de parada, inventários e custos estão em
  [P2_CAMPAIGN.md](../validation/diagnostics/chexal_riser_20261121/evidence/P2_CAMPAIGN.md).

| Caso | Chexal (efetivo 5) | Bhagwat (efetivo 3) |
|---|---|---|
| Stable BARO/ISO | COMPLETED; acomodação: p_in em 13/3 s, M_g riser em 111/152 s; desvio após 600 s ≤ 8e-9 | COMPLETED; 80/98 s |
| OSC BARO/ISO | BLOCKED_CONSTITUTIVE_DOMAIN 145,70/111,13 s (contracorrente inclinada no topo) | 2400 s decaindo ao estacionário / MAX_NEWTON 382,98 s |
| SS1 BARO/ISO | BLOCKED_JUNCTION_TRACE_DOMAIN 52,94/60,56 s | FAIL 3012 77,81 s / MAX_NEWTON 50,24 s |
| SS3 BARO/ISO | BLOCKED_CONSTITUTIVE_DOMAIN 41,55/50,73 s | FAIL 1017 97,97 s / 3012 72,63 s |

- **SS1 e SS3:** nenhuma sequência física completa (acumulação → bloqueio → penetração →
  produção → blowdown) foi alcançada com Chexal. Nada é declarado além do bloqueio.
- **OSC:** período, amplitude e taxa de crescimento não estão definidos, porque não há ciclo
  completo. A taxa de crescimento e o instante do bloqueio dependem fortemente da discretização
  (Parte I, §7: 93–181 s).

**Stable com Δt menor (campanha curta discriminante).** Com Δt = 0,5 s o próprio OSC Bhagwat
decai, então a estabilidade nesse passo não distingue física de dissipação numérica:

| Execução (BARO, Δt 0,125, 2400 s) | Resultado | Semiamplitude de p_in em janelas de 200 s (Pa) | Taxa de decaimento (200–600 s) | M_g riser estável (< 1e-4) | Máx. \|balanço\|/M₀ | CPU (s) |
|---|---|---|---:|---:|---:|---:|
| Chexal N51 | COMPLETED, contrato PASS | 127 → 2,3 → 0,094 → 3,2e-3 → 1,3e-4 → 0 | −0,016 s⁻¹ | 195 s | 8,7e-7 | 1534 |
| Chexal N101 | COMPLETED, contrato PASS | 72 → 2,1 → 0,15 → 9,4e-3 → 6,8e-4 → 0 | −0,013 s⁻¹ | 194 s | 6,9e-8 | 2509 |
| Bhagwat N51 | COMPLETED, contrato PASS | 102 → 0,073 → 1,8e-4 → 0 | −0,030 s⁻¹ | 105 s | 2,4e-7 | 891 |

Com Δt 0,125 o Stable continua decaindo ao estacionário com os dois fechamentos, sem
passos rejeitados. A estabilidade não é artefato do Δt 0,5. O Chexal é **menos amortecido**:
a taxa é metade da do Bhagwat e a acomodação ~2× mais lenta. A taxa varia 20 % entre N51 e
N101, logo ainda não está convergida na malha. A pressão de entrada final difere 42 Pa entre
as malhas (490 433 contra 490 391 Pa).
([P2_STABLE_DT.md](../validation/diagnostics/chexal_riser_20261121/evidence/P2_STABLE_DT.md))

## II.9 Balanços e convergência

- **Seis balanços** (duas fases × sistema, oleoduto e riser), com os fluxos aceitos
  (`audit_ledger`), em todas as execuções:
  - P2: máx. |balanço|/M₀ de 2,7e-12 a 3,0e-7 (critério 1e-6);
  - aritmética do ledger PASS 16/16;
  - integral na junção sem descasamento (0 kg).
- **Inventário inicial nulo:** não ocorre nestes casos; todas as fases começam presentes.
- **Convergência:** OSC e SS3 com Chexal não têm ciclos, então não há grandeza periódica para
  convergir. O que converge em tempo é o **instante de bloqueio do SS3** (N51: 41,55 → 41,95 →
  42,18 s, ordem aparente ≈ 0,8). O do OSC não converge.
- **`tol_res = 1e-10`** em todos os decks de referência.

## II.10 Validação com Baliño et al. (2010) (P2.2)

O PDF de Baliño, Burr e Nemoto (2010) não está nesta máquina; só o registro da conferência
anterior: Tabela 3, p. 653 — SS3 76 s; OSC 19,3–33,8 s; SS1 85,6 s; dissertação 26 s para
SS3, divergência preservada). Além disso, nenhuma trajetória Chexal tem ciclos completos e
admissíveis.

**P2.2: NOT_VERIFIED.** Não se compararam previsões numéricas publicadas com medidas.

## II.11 Matriz de qualificação (rodada 2)

| Item | Status | Evidência |
|---|---|---|
| P0.1 NP-3989: linhas 1, 4, 8 | PASS | §II.1, `test_chexal_np3989` |
| P0.1 NP-3989: linha 5 (diferença 1986/1992 demonstrada) | PASS | §II.1 |
| P0.1 NP-3989: linha 3 (referência inconsistente) | FAIL | identidade impressa −4,2e-4 |
| P0.1 NP-3989: linhas 2, 6, 7 (1000 psia) | FAIL (causa não identificada) | §II.1 |
| P0.2 ângulos e velocidades | PASS | `test_chexal_angles` (468) |
| P0.3 contrato do fechamento efetivo | PASS | `closure_contract.txt`, verificador, `test_closure_contract_cli` |
| P0.4 Jacobianos nas fronteiras | PASS | §II.4 |
| P1.1 validação experimental ar–água | NOT_VERIFIED (dados indisponíveis) | §II.5 |
| P1.2 contracorrente inclinada | BLOCKED | §II.6 |
| P1.3 junção: conservação, admissibilidade, Jacobiano, partição | PASS | §II.7 |
| P1.3 junção: condição de transmissão | BLOCKED (decisão T1/T2) | §II.7 |
| P2.1 replays com contrato (16 + 3) | PASS como execução e reprodução: contrato 19/19 e reprodução byte a byte. Física: Stable estável (Δt 0,5 e 0,125); OSC/SS3 BLOCKED_CONSTITUTIVE_DOMAIN; SS1 BLOCKED_JUNCTION_TRACE_DOMAIN | §II.8 |
| P2.2 Baliño et al. 2010 | NOT_VERIFIED | §II.10 |
| SS1, SS3, OSC validados | não (sem convergência nem confronto) | — |

**Commits locais da rodada 2** (branch `implementation/chexal-riser-20261121`, sem push nem merge):

| SHA | Conteúdo |
|---|---|
| `749a90b` | feat(core): fluido STEAM_WATER explícito para verificação; teste por linha do NP-3989 |
| `1525409` | test: certificado de ângulos e velocidades |
| `f468e67` | feat: contrato do fechamento (solver, checkpoint, verificador, CLI, CI, lançador) |
| `49fbd32` | diag: auditoria NP-3989 1986/1992, ajuste de propriedades, sonda de Jacobiano, mapa de domínio |
| `0a0159a` | test: invariância à partição também com Chexal |
| `7ba4d79` | evidence: auditoria, Jacobianos de fronteira, mapa de domínio, estudo local de contracorrente |
| (seguintes) | evidence P2 e docs |

Binários:
- `prims_f_contract`: `9e54e79f0abac7107b1a2f53108e22a498026399bcead789b5ad010eab26c778`;
- da Parte I: `d342750…`.

Compilador GNU Fortran 13.3.0; LAPACK/BLAS do sistema (`-llapack -lblas`); uma thread.

Testes:
- Fortran 28/28 (5961 asserções);
- Python 123 OK;
- CLI 7/7: `pipe_regime`, `accepted_telemetry`, `feed_contract`, `boundary_telemetry`, `initial_inventory`, `chexal_default` e `closure_contract`.

Custos medidos (CPU de uma thread):
- P2: 87–547 s por execução;
- Stable Δt 0,125: 891–2509 s;
- auditoria NP-3989: 3 min.

```sh
D=validation/diagnostics/chexal_riser_20261121
python3 -I $D/np3989_audit.py OUT.json OUT.md; python3 $D/np3989_property_fit.py OUT.json
./tmp/cx/boundary_jacobian_probe $D/evidence/global_jacobian/gj_D_OSC_C5_BAROTROPIC_N51_dt0.5_t60.params
python3 $D/validation_domain_map.py OUT.json OUT.md RESULTS/p1/runs/*C5*/
python3 $D/inclined_cc_local.py $D/evidence/DOMAIN_PROBE.json OUT.json OUT.md
python3 $D/run_campaign.py BIN RESULTS/p2 --workers 6          # verifies every run (contract)
python3 $D/run_campaign.py BIN RESULTS/p2_stable_dt --workers 3 --eos BAROTROPIC --closures 5 3 --cases D_Stable --grid 51:0.125
python3 scripts/verify_closure_contract.py RUN_DIR --deck DECK
python3 tests/test_closure_contract_cli.py --binary bin/prims_f_check_all --output tmp/ci_closure_contract
```
