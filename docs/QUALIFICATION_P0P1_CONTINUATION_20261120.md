# Continuação P0/P1 — 19–20/11/2026

Branch local `implementation/p0p1-continuation-20261119`, worktree
`~/Documentos/PRIMS-F-p0p1c-20261119`, base `f3717dd` (merge da PR #4; operador
`1603419`). A `main` local foi preservada; nenhum push ou merge. Dados brutos em
`~/Documentos/PRIMS-F-p0p1c-results-20261119` (fora do Git); índice com hashes em
`validation/diagnostics/p0p1c_20261119/evidence/`.

**Conclusão.** Os dois defeitos P0 do modo reduzido estão corrigidos e
reproduzidos antes/depois. A acomodação passou a exigir inventários, fluxos por
fase, holdups e formas de onda. O OSC `ISOTHERMAL`, o SS1 e o SS3 continuam
**BLOCKED**, agora com causa localizada: em SS1 e OSC‑ISO, selo líquido na base
do riser com fluxo gasoso de junção indo a zero; em SS3, refluxo de líquido
(`u_l → 0⁻`) na catenária. Os dois estados ficam fora do domínio constitutivo
declarado. Nenhum regime foi validado experimentalmente.

## Ambiente

| Item | Valor |
|---|---|
| Máquina | `leonardo@leonardo-ubuntu`, Linux 7.0.0-42, 8 núcleos |
| Compilador | GNU Fortran 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1), `-std=f2018` |
| BLAS/LAPACK | alternativas do sistema resolvem para **Intel MKL** (`libmkl_rt.so`); o registro histórico usava LAPACK de referência 3.12.0 |
| Python | 3.12.3, NumPy 1.26.4, SciPy 1.11.4, Matplotlib 3.6.3 |
| `AGENTS.md` | **ausente** no clone e em todos os ramos remotos |
| Skill | `rotodinamica-fem-fortran` carregada |
| Dados históricos | `/home/leonardoab/...` **ausente** nesta máquina; reanálise histórica bruta BLOCKED |

## Commits

| Commit | Conteúdo |
|---|---|
| `5e4882e` | P0‑B: busca mássica no ramo subsônico admissível |
| `1c778bf` | P0: média hidrostática normalizada pelos pesos de Gauss armazenados |
| `aef3774` | P0‑A: Jacobiano reduzido dentro do domínio admissível + regressões |
| `813da75` | critérios de acomodação/convergência v2 **pré-registrados** |
| `9552bab` | analisador com critérios v2 + testes sintéticos |
| `fab19c6` | opção `junction_trace='STEADY'` (opt‑in; default `CENTER` inalterado) |
| `2492670` | diagnósticos de raio/raízes e limites de fluxo na junção |

## P0 — modo reduzido

| Defeito | Causa | Correção | Teste | Evidência | Estado |
|---|---|---|---|---|---|
| A. Jacobiano no limite de refluxo | equilíbrio admissível (G=0, p=Pt, U3=0, alimentação nula, `outlet_reservoir_phase=0`); a perturbação negativa cria refluxo auxiliar → 934 → convertido em **935**; o integrador pede o Jacobiano antes de testar convergência | por coluna: diferença central com incremento controlado (h, h/8, h/64); senão, estêncil **unilateral de 2ª ordem** do lado admissível; senão, o **código original** da recusa; contagem exposta em `boundary_one_sided_count` | `test_reduced_domain_edges` (76 asserções); `test_reduced_equilibrium_cli` | antes: 4/4 equilíbrios (líquido/gás × POINT/HYDROSTATIC) colapsam dt em t=0, `935`, 17 rejeições; depois: 2 s, 57 passos, **0 rejeições**, debug e release. Unilateral vs referência central suave (reservatório da mesma fase declarado): 6,1·10⁻¹⁰ (líquido), 9,3·10⁻⁹ (gás) | **PASS** |
| B1. alimentação negativa | busca em `[nearest(0), 2ρ]` concluía "sem raiz" pelo sinal perto do vácuo | intervalo do ramo admissível: `m(r)=r·v_b(r)` é estritamente crescente no ramo subsônico; extremos sônicos em forma fechada `r₋=ρe^{−1−v/c}` e `√(r₊/ρ)−√(ρ/r₊)=(c−v)/c`; vácuo `p=0` como limite aberto | idem | `ρ=1,c=1,u=−0,1,alvo=−0,1`: antes **934**; depois `ρ_b=1` (≤4 ulp) | **PASS** |
| B2. `log(0)` | `nearest(0,+1)/ρ` arredonda para 0 | `r/ρ ∈ (e⁻², …)`, sem razão subnormal | idem, com `-ffpe-trap` | `ρ=2,u=0,1,alvo=0,2`: antes **SIGFPE** (debug); depois `ρ_b=2` | **PASS** |
| Achado: repouso HYDROSTATIC | `½Σw = 1+2,2·10⁻¹⁶` ⇒ `f_r=1−2,2·10⁻¹⁶`, `v_b=−2,2·10⁻¹⁶c`: o estado de repouso exato era recusado como refluxo | média `Σw·e^x/Σw` (exata para constantes) na reconstrução e na IC | idem | G=0 passa; nenhuma tolerância no sinal de `v_b`, nenhuma composição inferida | **PASS** |

Também cobertos: invariância de escala de densidade (λ=10⁻⁸…10⁶) e de velocidade
(κ=10⁻³, 10³), alimentação positiva/negativa/nula, equivalência `MASS × LOCAL_VOLUME`
(5 velocidades, densidade e fluxo ≤ 10⁻¹³), aproximação dos limites sônicos
(±10⁻¹⁰) e do vácuo líquido (±10⁻⁶), além dos estados explicitamente não suportados
(interior sônico → 934; pressão negativa → 935). O refluxo físico efetivo sem dado
externo continua recusado: coluna sob gravidade com fonte POINT termina com
propostas **934**, nunca 935 nem aceitas.

**Limite remanescente (não alterado):** com G≠0, a IC hidrostática exata ainda
produz `v_b = ±c·O(ε)` na saída; quando negativo, `INITIAL_OUTLET 934` (pré‑existente,
idêntico em `f3717dd`). Tratar esse refluxo de arredondamento como zero exigiria
escolher a composição de entrada. Próximo teste: declarar `outlet_reservoir_phase`
quando o separador contém a fase, ou definir a IC pela própria reconstrução de face.

## P1 — critérios de acomodação e convergência

Auditoria: `analyze_refinement.py` declarava `settled` só com período, amplitude e
média da pressão. Os critérios v2 foram **registrados antes de qualquer análise**
(`813da75`, SHA‑256 `84116c3e…`) e aplicados a 14 canais: pressões base/junção,
holdups base/saída, fluxos mássicos de líquido e gás na junção e na saída, e
inventários líquido/gás do sistema, do oleoduto e do riser.

| Critério | Regra (normalização `max(R̄, piso)`) | Justificativa |
|---|---|---|
| C0 (preservado) | `(máx−mín)/|média| ≤ 0,5%` para T, A e média da pressão | critério histórico |
| C1 deriva da média por ciclo | `máx(ptp m_k, 9·|inclinação|) ≤ 5%` | para a pressão, 5%·R̄ é da mesma ordem que 0,5%·média; para inventários (~70 kg) a média não mede deriva |
| C2 faixa por ciclo | `ptp r_k ≤ 1%` | 2× o limite de amplitude, pois os extremos são amostrados |
| C3 forma | `máx_k RMS(y_k − ȳ) ≤ 2%` | jitter permitido por C0 ≈ 0,6% para senoide; margem ~3 |
| Convergência | período e amplitude ≤ 1% (preservados); nível, faixa e forma **em fase fixa** por canal | o atraso ótimo é só diagnóstico de erro de fase |

Pisos dimensionados: 1 Pa; holdup 10⁻⁴; fluxo 10⁻⁴·Ṁ_alimentação; inventário
10⁻⁴·Ṁ·T. Canal ausente ou não avaliável (faixa e piso nulos) reprova. Onze
marcadores distintos são exigidos para dez ciclos.

Testes sintéticos (13): pressão acomodada com inventário do riser derivando (incluindo
deriva de 0,3% da média, invisível a uma regra normalizada pela média) → reprova;
forma alternando ciclo a ciclo com média e faixa iguais → reprova só em forma;
canal de média quase nula usa a faixa; canal constante usa o piso; alimentação nula
→ não avaliável; deslocamento de 5% do período não é ocultado pelo alinhamento.

**Reanálise:** os brutos históricos estão ausentes. Regenerei a execução N51/dt 0,125
(600 s) a partir de `common_initial_U.txt`: T = 23,108397263824 s, A = 9905,699630720 Pa,
média idêntica a 10⁻¹⁴, os mesmos 4822 passos — idêntica a ~10⁻¹³ relativo, não bit a bit.
Com os critérios v2, os **14 canais passam** (deriva ≤ 10⁻⁴, faixa ≤ 1,4·10⁻³,
forma ≤ 3·10⁻⁴). Regenerei também N51/dt 0,0625 (900 s direto, 14418 passos, 0
rejeições): T = 23,265321 s, A = 13764,28 Pa, contra 13766,38 Pa na história composta
anterior, que teve um reinício. A diferença de 0,015% não foi reconciliada bit a bit.
As duas acomodam com v2 (14/14 canais).

Convergência v2 no eixo temporal (dt 0,125 → 0,0625, N51), em fase fixa:

| Grandeza | Diferença | Alvo | Estado |
|---|---:|---:|---|
| Período | 0,67% | 1% | PASS |
| Semi‑amplitude da pressão | **28,0%** | 1% | **FAIL** |
| Faixa por canal (14) | 27–41% | 1% | FAIL em todos |
| Forma em fase fixa (14) | 10–15% | 2% | FAIL em todos |
| Nível médio (14) | ≤ 1,8% | 1% | 11 PASS / 3 FAIL (holdups base/saída, Mg do riser) |
| Atraso de fase (diagnóstico) | ≈ −1,2% do período | — | não usado para aprovar |

Isso confirma, com critérios mais amplos, a falta de convergência já registrada; não
se recalcularam N101–N401 nem dt 0,03125 (≳ 40 h de CPU nesta máquina).

## P1 — OSC ISOTHERMAL

EOS conferida no código (`c²=RT`, 290,06 m/s com R=287 J/(kg·K) e T=293,15 K;
conversão padrão `P_std/(R·T_std)`; `test_feed_contracts`). Reprodução **bit a bit**
do bloqueio histórico: MASS para em `t=223,64452291713980 s`, com 1829 aceitos e 25
rejeitados (`MAX_NEWTON`); LOCAL_VOLUME em 223,64197 s. O deck LOCAL_VOLUME com `mul`
explícito é recusado (`121`, contrato conflitante), corretamente.

Causa no operador efetivo (replay do estado terminal, Jacobiano local exato
1,4·10⁻¹¹):

1. A oscilação **cresce**: faixa de pressão 272–273 kPa (0–50 s) → 250–292 kPa (200–215 s).
   A fração líquida na base do riser sobe para 0,994 e o fluxo gasoso de junção decai
   linearmente a zero (~1,5·10⁻⁴ kg/s por s). O oleoduto termina estratificado
   (α≈0,55, u_g≈0,057 m/s); a base do riser tem α≈0,0057.
2. O passo de Newton completo pede fluxo gasoso de face **negativo** →
   `3012 = 3000 + 1000 + CLOSURE_DOMAIN(12)` no fechamento cocorrente do riser.
   Na parte admissível do raio o resíduo não cai (razão 0,9998). Isso não é falha
   de busca nem de globalização.
3. O fluxo gasoso Roe de face (5·10⁻¹³) fica **fora** do intervalo dos fluxos dos dois
   traços [0,017; 0,100] kg/m²s; 84% da dissipação gasosa vem da diferença de pressão
   de equilíbrio entre os centros (73 Pa), lida como salto acústico.

Teste discriminante (o instante de bloqueio converge com h/dt?):

| ISOTHERMAL/MASS | t de bloqueio (s) | local |
|---|---:|---|
| N51, dt 0,125 | 223,645 | face 53, 3012 |
| N51, dt 0,0625 | 222,518 | face 53, 3012 |
| N51, dt 0,03125 | 198,667 | face 53, 3012 |
| N101, dt 0,125 | 249,400 | face 103, 3012 |
| N51, dt 0,125, `junction_trace='STEADY'` | 223,650 | face 53, 3012 |

O evento (selo líquido e corte de gás na junção) é **robusto** em todas as
discretizações; o instante **não converge**. A reconstrução estacionária da junção
remove o artefato, mas não altera a trajetória. Não há ciclo limite acomodado:
dez ciclos úteis, refino e convergência de amplitude ISOTHERMAL ficam **BLOCKED**.
Nenhum U de outra EOS foi reutilizado. Resultados BAROTROPIC permanecem apenas
regressão e comparação.

## P1 — SS1/MATCHED e SS3

Mapa ao longo do raio de Newton nos estados terminais (`ray_map.f90`,
`junction_bounds.f90`, `steady_face_check.f90`; operador atual, sem mudança de física):

| Caso | Primeira inadmissibilidade | Mecanismo | Classificação |
|---|---|---|---|
| SS1 (78,58 s) | λ = 9,98·10⁻⁹, face 53, 3012 | Gg Roe de face 5·10⁻¹³, fora de [0,040; 0,071]; o Newton pede Gg<0; o resíduo não cai na parte admissível | ausência de estado no domínio do fechamento cocorrente (não busca, não globalização) |
| SS3 (97,13 s) | λ = 4,15·10⁻⁶, célula 69, 1017 | raiz de Φ(u_g) **na extremidade** `u_g=u₃/u₂` ⇔ **u_l=0** (|Φ|/α = 7,7·10⁻⁸; 0 trocas de sinal internas); o Newton pede u_l<0 | ausência **provada** de raiz no domínio declarado (Bhagwat cocorrente); não é falha de busca |
| OSC‑ISO (223,64 s) | λ = 3,37·10⁻⁸, face 53, 3012 | como SS1 | idem |

Decomposição característica da dissipação gasosa na junção (fluxo gasoso negativo é
fluxo assinado, não massa negativa; proposta recusada não é contrafluxo observado):

| Caso | Δp centro‑a‑centro | dissipação gasosa total | parcela de Δp |
|---|---:|---:|---:|
| SS1 | 57,0 Pa | −0,0557 kg/m²s | −0,0337 (61%) |
| OSC‑ISO | 73,1 Pa | −0,0585 | −0,0490 (84%) |
| SS3 | 19,6 Pa | −0,0851 | −0,0200 (23%) |

Correção implementada para essa parcela: `junction_trace='STEADY'` (opt‑in) reconstrói
cada traço até a face da junção com a mesma EDO estacionária (`fun_dpds`: gravidade,
atrito, inclinação, fechamento do próprio trecho), RK4 com 4 passos fixos, mantendo
um fluxo por fase. Perfil estacionário independente 9+9:

| Fechamento do riser | erro relativo Gg* `CENTER` | `STEADY` | p* `CENTER` → `STEADY` |
|---|---:|---:|---:|
| Bendiksen | 29,7% | 5,4·10⁻⁸ | 1,5·10⁻³ → 1,4·10⁻¹⁰ |
| Bhagwat | 28,7% | 4,9·10⁻⁸ | 1,5·10⁻³ → 1,3·10⁻¹⁰ |

O erro do traço escala como h⁴ (1/2/4/16 passos: 1,5·10⁻⁵ / 8,8·10⁻⁷ / 5,4·10⁻⁸ /
7·10⁻¹¹). A telescopagem é de 10⁻²² e o Jacobiano direcional tem erro de
5·10⁻⁷/1·10⁻⁷. O default é byte a byte idêntico a `f3717dd` (D_OSC e D_Stable, 40 s:
monitor, conservation, boundary, fields, rejected, checkpoint). O custo sobe ~4×
por passo.

Efeito dinâmico: com o protocolo histórico, `CENTER` para em 78,5766 s (706/29;
histórico 78,5768 s, 711/30) e `STEADY` em **78,8378 s**, ambos na face 53/3012. No
novo estado terminal, as duas células ainda transportam gás (0,033 e 0,039 kg/m²s),
mas o Gg de face `STEADY` volta a 2·10⁻⁸: resta a dissipação da onda material (salto
de α entre o traço do oleoduto e o do riser mapeado, quando os fluxos de fase dos dois
lados diferem). **SS1 BLOCKED; SS3 BLOCKED.** A opção fica qualificada apenas para
*well‑balancing* estacionário, sem benefício dinâmico demonstrado.

Próximos testes discriminantes, sem alterar a física para fabricar aprovação:

1. Fluxo de junção que respeite limites de fluxo por fase (positividade da massa
   gasosa transportada entre traços com u_g>0), verificado no mesmo snapshot antes de
   qualquer promoção. HLL/Rusanov **pioram** (Gg de face −3,2 e −1,2 kg/m²s).
2. Modelo de face para gás nulo sob líquido ascendente (traço monofásico líquido na
   face com células bifásicas): exige transição misto–puro conservativa, hoje BLOCKED.
3. SS3: fechamento contracorrente aplicável à catenária inclinada **com fonte
   primária**; Chexal 1992 só cobre contracorrente vertical e não foi estendido.

## Verificação

| Item | Resultado |
|---|---|
| Builds Fortran 2018 release e `check_all` | PASS (warnings preexistentes) |
| Suíte Fortran (`run_unit_tests.py`, debug com `-fcheck=all -ffpe-trap`) | **22/22 PASS** (20 preexistentes + 2 novos) |
| Python `unittest discover` | **123 OK** (110 + 13) |
| CLI | boundary telemetry, feed contract, initial inventory, accepted telemetry, pipe regime, reduced equilibrium: **6/6 PASS** |
| Default byte a byte vs `f3717dd` | PASS (D_OSC/D_Stable 40 s) |
| Balanços (2 fases × sistema/oleoduto/riser, ledger aceito, pesos BE) | máximo/M₀ ≤ 1,09·10⁻⁹ em todas as execuções; cancelamento da junção exato |
| `tol_res` | `1e-10` mantido nos protocolos; defaults não alterados |

Custos (CPU): OSC‑ISO 311–815 s por execução; SS1 `CENTER`/`STEADY` 213/629 s;
barotrópico N51/dt 0,125 600 s: 2268 s; N51/dt 0,0625 900 s: 6482 s. TIMEOUT: nenhum. Três processos foram
encerrados por mim antes do fim (SS1 e OSC‑ISO `STEADY` com 16 passos; SS1 `STEADY`
de 4 passos com deck fora do protocolo); foram substituídos e não são usados como
evidência. A primeira execução SS1 `CENTER` com o deck default (tol 1e-8, dt_max
0,5 s; 77,81 s) também fica fora da comparação, por não seguir o protocolo.

## Estado

| Item | Estado |
|---|---|
| P0‑A Jacobiano reduzido | **PASS** |
| P0‑B busca mássica | **PASS** |
| Critérios v2 implementados e testados | **PASS** |
| Reanálise histórica dos brutos | **BLOCKED** (dados ausentes) |
| Acomodação v2 em execução regenerada | PASS (N51/dt 0,125) |
| Convergência v2 h/dt (amplitude ≤ 1%) | **FAIL** no eixo dt (28%); eixo h NOT_VERIFIED nesta continuação |
| OSC ISOTHERMAL: EOS/contratos | PASS |
| OSC ISOTHERMAL: dez ciclos, refino | **BLOCKED** (selo líquido; instante não convergente) |
| SS1 | **BLOCKED** |
| SS3 | **BLOCKED** |
| Junção well‑balanced (estacionário) | PASS (opt‑in); desbloqueio **FAIL** |
| Validação experimental (Baliño, Burr e Nemoto 2010) | **NOT_VERIFIED** |
