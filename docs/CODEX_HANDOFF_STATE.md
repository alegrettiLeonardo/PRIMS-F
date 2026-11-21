## Atualização: domínio físico Chexal e transmissão na junção, 21/11/2026 (rodada 3)

Relatório: [CHEXAL_PHYSICAL_DOMAIN_AND_JUNCTION_QUALIFICATION.md](CHEXAL_PHYSICAL_DOMAIN_AND_JUNCTION_QUALIFICATION.md).

- **P0.1:**
  - contrato executável ar–água no núcleo (`CHEXAL_UNQUALIFIED` = 68 para ρg/ρl > 8,1e-3
    ou D ∉ [0,010; 0,305] m);
  - EPRI 1000 psia decomposto: as linhas 2, 6 e 7 implicam as mesmas propriedades
    (ρg/ρf = 0,0468, grupo de Vgj +2,5 %), com equações verificadas;
  - FAIL mantido para uso generalizado; subconjunto qualificado = ramo C5 < 1 (todo o PRIMS-F).
- **P0.2: recomendação T1.**
  - Fundamentos: Burr (72)–(75) e discretização nodal; balanço integral; dobra sem força
    tangencial; Riemann local com Chexal (8 estados, multistart).
  - No bloqueio do SS1, a T1 tem solução única com G_g* = 2,94e-3; a T2 não tem solução
    (só G_g* ≤ 0).
  - O operador atual ≈ T2. Ainda **não implementado**: próximo passo, T1 opt-in com
    Jacobiano, partição e replay curto do SS1.
- **P1.1:** registro rastreável de 9 fontes em `validation/experimental_airwater/`; nenhum
  ponto primário de α acessível: NOT_VERIFIED.
- **P1.2:**
  - mecanismo B (líquido descendo a 65–76°): BLOCKED (Ousaka et al. 2006 é limite de
    inundação, 16/26 mm, 30–60°);
  - mecanismo A: não exigido sob T1.
- **Regressão:** 29/29 Fortran, Python OK, 7/7 CLI; 16 réplicas idênticas byte a byte à P2,
  contrato 16/16.

## Atualização: qualificação Chexal, rodada 2 (P0/P1/P2), 21/11/2026

Mesma branch e worktree, sem push nem merge. Relatório: Parte II de
[CHEXAL_LELLOUCHE_RISER_QUALIFICATION.md](CHEXAL_LELLOUCHE_RISER_QUALIFICATION.md).

- **P0.1 NP-3989:**
  - linhas 1, 4 e 8 PASS;
  - linha 5 PASS só com as equações de 1986 (revisão descendente demonstrada);
  - linha 3 é referência inconsistente;
  - linhas 2, 6 e 7 (1000 psia) FAIL, com causa não identificada: não é versão nem Re, e propriedades só a explicam com um par incoerente.

  O núcleo ganhou `fluid=CHEXAL_FLUID_STEAM_WATER` (só para verificação; o caminho ar–água é bit a bit).
- **P0.2:** `test_chexal_angles`, 468 verificações, PASS.
- **P0.3:**
  - `closure_contract.txt` (fechamento solicitado, efetivo e despachado por célula), que para o solver se divergirem;
  - checkpoint com os IDs de fechamento;
  - `scripts/verify_closure_contract.py` e `tests/test_closure_contract_cli.py` (CI);
  - `run_campaign.py` verifica cada execução.
- **P0.4:** `boundary_jacobian_probe`:
  - erro ≤ 9e-10 em j_l → 0⁺ no topo;
  - até 2,2e-6 com α → 0;
  - nas quinas, o Jacobiano é a derivada do pedaço do estado; na quina exata, status 53 (`CLOSURE_DERIVATIVE`).
- **P1.1:** NOT_VERIFIED. Não há pontos experimentais. O mapa de domínio dá 0 % validado (D = 52,5 mm, fora dos 25–38 mm de Beggs) e interpolação por Fr.
- **P1.2:** BLOCKED. Nos 11 estados de bloqueio, os candidatos diferem 1,4–3,5× em u_g. Barnea e Ohnuki não estão disponíveis localmente.
- **P1.3:** a junção passa em conservação, Jacobiano e partição (agora também com Chexal). A transmissão T1/T2 está pendente de decisão.
- **P2:**
  - 16 execuções com o binário `prims_f_contract` (`9e54e79f…`): contrato 16/16 PASS e saídas idênticas byte a byte às da Parte I;
  - Stable a Δt 0,125 também decai (Chexal mais devagar).
- **Pendências científicas:** dados de C0/Vgj/CCFL inclinados ar–água com D ≈ 50 mm; decisão T1/T2; dados primários de fração de vazio; PDF de Baliño et al. 2010.

## Atualização: Chexal–Lellouche em todo o riser, 21/11/2026

Worktree `PRIMS-F-chexal-riser-20261121`, branch `implementation/chexal-riser-20261121`
(base `d2c8101`); sem push/merge; `main` intocada. Relatório:
[CHEXAL_LELLOUCHE_RISER_QUALIFICATION.md](CHEXAL_LELLOUCHE_RISER_QUALIFICATION.md).
Resultados: `~/Documentos/PRIMS-F-chexal-riser-results-20261121` (`p1/`, `p1_refine/`,
`regress_bhagwat/`; `p1_INVALID_deck_override_bug/` não é evidência).

- **Default:** `closure_riser=5` (CHEXAL_1992 AUTO) em `main.f90`, `sim_config_t` e
  decks primários. Bhagwat é o ID 3 explícito, com decks arquivados em
  `input/benchmarks/bhagwat/`, e reproduz byte a byte o binário `364b808`.
  ACTIVE_SET, E1 e E2 são recusados com Chexal.
- **Verificação:**
  - Fortran 26/26 (5947 asserções), Python 123 OK, seis CLI PASS e o CLI novo
    `test_chexal_default_cli.py`;
  - núcleo contra transcrição independente em 945 estados (≤ 6,9e-14);
  - NP-3989 Tabela A-1: casos de baixa pressão PASS, 1000 psia não reproduzidos;
  - Jacobiano global contra diferenças finitas em 6 estados aceitos (≤ 4,3e-8);
  - todos os estados aceitos do riser satisfazem o Chexal (≤ 8,3e-14).
- **Física:**
  - Stable completo (BARO/ISO).
  - OSC e SS3 → `BLOCKED_CONSTITUTIVE_DOMAIN`: contracorrente inclinada exigida no topo
    (65–76°), j_l → 0; persiste sob refino.
  - SS1 → `BLOCKED_JUNCTION_TRACE_DOMAIN`: gás precisa inverter no traço horizontal.
  - SS3 com N101/N201 a Δt 0,5 também falha por quina da derivada Chexal no caminho de
    Roe (9323).
  - Sem ciclos com Chexal: validação experimental NOT_VERIFIED; promoção BLOCKED.
- **Próximo passo científico:** dados e forma constitutiva de contracorrente inclinada
  ar–água (§9 do relatório). Não estender Chexal 1992 sem essa fonte.

## Atualização: auditoria avançada P0–P2, 20/11/2026

Mesma branch/worktree; sem push/merge. Relatório:
[QUALIFICATION_P0P2_ADVANCED_20261120.md](QUALIFICATION_P0P2_ADVANCED_20261120.md).
Operador congelado `007359d` (`frozen_007359d/`, `prims_f` `d1184263…`); default byte a
byte = `1d71f00`; 25/25 Fortran, 123 Python, 6/6 CLI.

- **P0.1:** nos estados aceitos terminais, a junção é classe **C** (gás presente e
  subindo, fluxo de face nulo por cancelamento Roe). Com E1/STEADY, o gás no fim do
  oleoduto chega a 1,3–1,6% de |U_d|, perto da estagnação (A). Ferramentas:
  `junction_state_compare.f90` (argumentos: fixture, c_g, fluxo e traço da execução;
  β = 2° agora definido) e `junction_classification.py`.
- **P0.2:** `junction_flux='PHASE_UPWIND'` (E1, opt‑in, experimental, **não promover**).
  Vetor (G_l\*, G_g\*, p\*) por Riemann acústico de Wood mais separação de fluxo
  material por fase. A dinâmica fica igual à do Roe até o bloqueio; todas as execuções
  concluídas terminam em 3021 (traço do riser com j_g = 0). O próximo passo é o E2
  (passagem misto→puro na face).
- **P0.3:** `initial_u_file` (opt‑in), `make_perturbed_ic.py`, `p03_analysis.py`.
  - Estacionário discreto por Euler implícito com Δt grande (`p03/SS_N*`); com ΔP = 0
    fica bit a bit estacionário.
  - Com ΔP = 250 Pa, a semente independe da malha e os níveis refinados falham no
    mesmo vale (173,3–174,2 s).
  - σ ainda não está na região assintótica.
- **P1.1:** `ss3_closure_map.f90`. O ramo Bhagwat termina em u_l = 0⁺, na fronteira do
  domínio. Bendiksen dá um ramo contracorrente contínuo e hiperbólico, mas por
  extrapolação. Literatura primária só no nível bibliográfico.
- **P1.2:** combinado N101/Δt 0,0625 concluído. Amplitude FAIL (19–22% estimado);
  interação espaço–tempo −619 Pa.
- **P2:** NOT_VERIFIED. O PDF de Baliño 2010 não está nesta máquina; valores só via
  `benchmark_sources.json`.

---

## Atualização: continuação 2 da auditoria físico-numérica, 20/11/2026

Mesma branch/worktree; sem push/merge. Relatório:
[QUALIFICATION_P0P1_CONTINUATION2_20261120.md](QUALIFICATION_P0P1_CONTINUATION2_20261120.md).

- P0.1 (`0960246`): face do separador do modo reduzido *well-balanced* com G≠0
  (forma de desvio do equilíbrio discreto da célula; `log1p` estável). Repouso
  exato preservado; refluxo real sem composição continua 934. `test_hydrostatic_outlet`
  (1760 asserções).
- P0.2 (`1dde473`, `9c31dea`): certificados SS1/OSC-ISO/SS3 em
  `validation/diagnostics/p0p1c_20261119/evidence/certificates/`. Todos em proposta de
  Newton a partir de aceitos admissíveis. SS1/OSC: Gg Roe de face cancelado pela
  dissipação com os dois traços a u_g>0 (`NUMERICAL_FLUX_INDUCED_BOUNDARY`). SS3: raiz
  só em u_l=0 no ramo cocorrente.
- P1.1 BLOCKED por limitação documentada (traço exige j_g>0; 21 em j_g=0, 12 em
  j_g<0; `test_junction_vanishing_gas`). Extensão E1 (fluxo de massa por fase
  *upwind*, tipo AUSMV) e E2 (limite α*→0⁺ só como continuação) **propostas, não
  implementadas**.
- P1.2: SS3 sem fechamento qualificado; Chexal 1992 cobre contracorrente só na
  vertical (PDF local, p. 8). `ss3_local_closures.f90` é o solver local assinado.
- P1.3 ISO (`1d71f00`, `1518f5d`; N51–N401, dt 0,125–0,015625): período converge;
  σ assintótico de ordem ≈ 1 (erro 1,2–1,5%); semente O(h^0,8); A*(f) ±2,3%; o
  instante do evento salta por ciclos (FAIL, não declarado).
- P1.3 BARO (`21ad166`): as 4 execuções concluídas acomodam (v2, 14 canais);
  amplitude FAIL (h 22%; dt 28% → 10%, ordem 1,28, erro estimado 6,8%). O combinado
  N101/dt 0,0625 estava em execução em
  `~/Documentos/PRIMS-F-p0p1c-results-20261119/baro_combined/`; analisar com
  `baro_refinement_v2.py` quando houver `run_status.txt`.
- Regressões no `1d71f00`: 24/24 Fortran, 123 Python, 6/6 CLI; binários de campanha
  byte a byte equivalentes ao HEAD.

---

## Atualização: continuação P0/P1, 19–20/11/2026

Branch local `implementation/p0p1-continuation-20261119`, worktree
`~/Documentos/PRIMS-F-p0p1c-20261119` (máquina `leonardo-ubuntu`), base `f3717dd`;
main local intacta, sem push/merge. Relatório:
[QUALIFICATION_P0P1_CONTINUATION_20261120.md](QUALIFICATION_P0P1_CONTINUATION_20261120.md).
Dados brutos: `~/Documentos/PRIMS-F-p0p1c-results-20261119`. `AGENTS.md` não existe
em nenhum ramo; os resultados `/home/leonardoab/...` não estão nesta máquina.

- P0 corrigidos: Jacobiano reduzido no limite de refluxo (`aef3774`), busca
  mássica no ramo subsônico (`5e4882e`), média hidrostática exata para constantes
  (`1c778bf`). 22/22 Fortran, 123 Python, 6 CLI; default byte a byte = `f3717dd`.
- Critérios v2 pré-registrados (`813da75`) e aplicados (`9552bab`): inventários,
  fluxos por fase, holdups e formas de onda entram em `settled`.
- OSC ISOTHERMAL reproduz bit a bit o bloqueio em 223,6445 s: selo líquido na base
  do riser e Gg de junção → 0 → CLOSURE_DOMAIN. Instante não converge em h/dt
  (198,7–249,4 s). BLOCKED.
- SS1 (face 53) e SS3 (célula 69, raiz em u_l=0) mapeados: ausência de estado no
  domínio constitutivo declarado, não falha de busca/globalização. BLOCKED.
- `junction_trace='STEADY'` (`fab19c6`, opt-in): remove erro de 29% do Gg
  estacionário na junção; não desbloqueia SS1 (78,84 s) nem OSC-ISO (223,65 s).
- BLAS/LAPACK desta máquina = Intel MKL via alternativas.

---

## Atualização: correções P0 e qualificação P1,19/11/2026

Branch `audit/p0p1-20261119`, worktree `PRIMS-F-p0p1-20261119`, base remota
`00df4c2`; main local4572f95 e auditoria6019277 preservadas. Operador final
**1603419806f5cc2db62fd2e326877ee232bd8df0**. Commits posteriores são
ferramentas/evidência/documentação, sem mudança de scr/Makefile.

Leia [QUALIFICATION_P0P1.md](QUALIFICATION_P0P1.md) e
[índice de ferramentas](../validation/diagnostics/p0p1_20261119/README.md).
Implementados: partição homogênea, EOS/feed explícitos, IC local-volume,
traço com fontes e modelo puro optativos, inventário inicial declarado,
composição externa no refluxo, guarda de transição932, telemetria aceita.
Defaults preservados;20 programas Fortran debug/release e110 Python PASS,
mais CLI e comparações byte a byte. Último OSC/MASS1603419 reproduz18aa7f7
até2400 s; manter SHA de cada evidência, não reatribuí-las.

Controle2400 s e seis balanços tol_res1e-10 PASS; tol_res1e-8 FAIL, sem
relaxar1e-6/M0. EOS isotérmica passa no controle, mas OSC falha~223,64 s.
SS1 BLOCKED~78,5768 s FACE53; SS3~97,1319 s CELL69. Traço com fontes antecipa
a falha SS1~51,1994 s: não promovê-lo como cura nem atribuir causa única ao BC.

Campanha comum nova congelada18aa7f7: N51/101/201/401 por trecho, três dt e
caso cruzado. Sete histórias concluídas com dez ciclos acomodados; N401 até900 s
da nova IC, não2400 s recalculados por malha. Últimas diferenças de amplitude
4,0189% em h e10,3455% em dt: **FAIL** ao alvo1%; período0,1702%/0,3372%.
42 balanços/M0≤9,442e-9. TIMEOUTs e caudas brutas preservados; composição
de análise auditada sem status fictício. Campos/frentes geométricos não provam
severe slugging; validação física SS1/SS3/OSC NOT_VERIFIED.

Resultados completos: `/home/leonardoab/Documentos/PRIMS-F-p0p1-results-20261119`.
Análise final: `refinement_composite_final_v2`, `refinement_analysis_final_v2`,
`refinement_spatial_final_v2`; índice de77 execuções e hashes em delivery_run_index.
Não misturar05ad799/490936a. Chexal contrafluxo inclinado, comutação BURR15,
transição misto–puro e promoção MUSCL/SDIRK2 continuam bloqueados/não qualificados.

Usuário autorizou publicar MR ao final: abrir **draft para main**, sem merge.
A restrição anterior de não publicar referia-se à etapa histórica abaixo.

---

## Atualização: auditoria P0 dos contornos, 19/11/2026

Branch local `audit/boundary-contract-20261119`, base `2c3e3af`; worktree
`/home/leonardoab/Documentos/PRIMS-F-boundary-audit-20261119`.
Escopo vigente: inspeção/instrumentação/testes/protótipos; **não alterar solver,
física ou defaults sem aprovação**. Nenhum push/merge nesta etapa; main intacta.

Leia primeiro [BOUNDARY_CONTRACT_AUDIT_20261119.md](BOUNDARY_CONTRACT_AUDIT_20261119.md)
e seus detalhes. Produção e Makefile não mudaram; binário SHA256
02ce51462ec3c954cec9364b17444b7279be72ec150f87e307f1b5fc8dd0a4c3.
Protocolo:08cd4a5; instrumentos:4ffe7b4; controles/manifestos:957eaa4.

Demonstrados: massa líquida fixa difere do contrato Ql com EOS compressível;
controle hidrostático não preservado pela pressão central usada na face;
MATCHED horizontal desnecessário quebra invariância da partição no momento.
Protótipos volume/partition corrigem apenas seus controles; permanecem externos.
Jacobiano completo de contorno passou; contornos externos são admissíveis nas
propostas dos replays que falham em FACE53(SS1) ou CELL69(SS3). Não inferir causa
histórica/global dos contornos a partir dessa localização.

Dados: `/home/leonardoab/Documentos/PRIMS-F-boundary-results-20261119`.
Onze execuções com manifestos; janelas OSC120s v1/v2, checkpoints atômicos,
983 estados exatos de extremidade por variante; seis balanços<=1,28e-8.
Semi-excursão muda−0,01059% no protótipo Q; só quatro intervalos entre picos.
Não atribuir essa evidência a dez ciclos, convergência ou validação física.
SS1/SS3 continuam BLOCKED; OSC não ganhou qualificação experimental.
Reversão de reservatório, NSCBC/Riemann não linear com fontes e evolução
exatamente monofásica continuam abertos. Não extrapolar Chexal à catenária.

Relatórios/campanhas anteriores abaixo são preservados com seus próprios SHAs.

---

# Estado para continuação no Codex


## Atualização da qualificação — 19/11/2026

- Continuação local: `audit/qualification2-20261119`, worktree
  `/home/leonardoab/Documentos/PRIMS-F-qualification2-20261119`.
- Leia [QUALIFICATION2_RESULTS_20261119.md](QUALIFICATION2_RESULTS_20261119.md)
  e o índice `validation/diagnostics/qualification2_20261119/delivery/index.json`.
- Operador novo `a179e51`: globalização ACTIVE_SET opt-in, experimental;
  LINE_SEARCH segue padrão. Release/debug, 15 programas Fortran, 103 testes
  Python dos analisadores, 21 checks CLI e 5 de telemetria aprovados.
- Controle estável novo 2400 s reproduz byte a byte os campos e seis balanços
  do 69813ee/tol 1e-10. SS1/SS3 seguem BLOCKED: com tol 1e-10 a conservação
  passa nas trajetórias registradas, mas SS1 sofre TIMEOUT perto 78,578 s e
  SS3 para perto 97,132 s. Replays específicos não qualificam os regimes.
- OSC usa freeze separado 69813ee: N51/101/201/401 e dt.125/.0625/.03125
  completaram 3600 s, dez intervalos acomodados e balanços PASS. Amplitude
  ainda muda 4,013% com h e10,318% com dt: convergência 1% FAIL.
- Seis sensibilidades estáticas + duas extensões 7200 preservadas. jg+
  com tol 1e-8 falha inventário; estudo separado tol 1e-10 passa (3,383e-8),
  sem período qualificado. AUTO_BURR15 continua apenas diagnóstico.
- Nenhum regime foi validado experimentalmente. Convergência, frentes
  materiais/observação experimental e continuação de domínio seguem abertas.
  Não extrapolar Chexal vertical; não escolher parâmetros para forçar ciclos.
- Todos os solvers/campanhas desta continuação terminaram; dados externos
  em `/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119`.
- PR#2 foi incorporado externamente às 15:38:32UTC, merge e4c6aeb, árvore
  igual a ee 5810e. Não houve push/merge por esta execução; main local 4572f95
  foi preservada. O PR fechado não pode continuar draft.

O estado histórico abaixo permanece preservado.


## Onde continuar e o que foi preservado

Repositório `alegrettiLeonardo/PRIMS-F`; PR #1 incorporado em `4572f95`.
Base desta continuação: `origin/audit/chexal-continuation-20261118`, `05ad799`.
O checkout `/home/leonardoab/Documentos/PRIMS-F` permaneceu limpo em `main`.
As alterações estão no worktree `/home/leonardoab/Documentos/PRIMS-F-chexal`,
branch `audit/chexal-implementation-20261118`. Confira status/refs antes de agir;
não faça reset nem descarte alterações locais.

No Windows, o endereço é `\\wsl.localhost\weg-linux\home\leonardoab\Documentos`.
Comandos usados: `wsl -d weg-linux --user leonardoab --exec sh -c ...`.
O perfil de login `sh -lc` não foi necessário; a distribuição Ubuntu-24.04
não contém este projeto. O fetch WSL falhou na rede; a referência foi trazida
de clone GitHub verificado no Windows, sem alterar `main`.

## Correções e implementação já realizadas

- `47433f8`: alpha-tail integrado **antes** de mudar física. O0/O3: 55 guards,
  398 checks estacionários por versão, 101 registros antigos idênticos; revisão
  O2: 129 contratos. Erro relativo máximo nas raízes novas `3,263e-13`.
  Não amplia o domínio de Bhagwat. Evidência em
  `validation/diagnostics/alpha_tail_integration_20261118`.
- `db8d47c`: `prior_evidence.json` recuperado do ZIP local com SHA original
  `f88afc885178f1c3456ed0d03c9177566939b7ab90edf77ce77f411ae169af89`.
  A ausência inicial da branch permanece registrada como lacuna resolvida.
- `ee4b7c6`: Roe interpola `(U,theta)` no caminho; MATCHED continua no ângulo
  do oleoduto. 315 checks específicos; snapshot SS3 antigo status14 passa0.
  A identidade inclui `integral(F_theta*dtheta)`; não se alega propriedade Roe
  exata para a integral de `F_U` isolada. Evidência em `roe_geometry_20261118`.
- `3096c3b`: núcleo Chexal1992 em SI, ar–água, Reynolds assinados, CCFL e ramos,
  referências independentes mpmath60. Núcleo:139 checks diagnósticos.
- `40a7d44`: integração de residual, primitivas, derivadas/Jacobianos,
  estacionário, contorno e junção; `868eeb9`: fixture adicional do ramo baixo.
  `closure_riser=5` é CHEXAL_1992, `6` HIGH, `7` LOW; Bhagwat permanece3.
- Análise: seleção explícita de coortes OSC e bloqueio de mistura de fontes,
  fechamento/EOS/parâmetros. NumPy1/2 suportados no integrador trapezoidal.

- `b0d33cb`: derivada Chexal em velocidade positiva minúscula, com stencil
  unilateral na mesma peça e guarda bilateral adaptativa. Referência independente
  mpmath60, 87 checks de cauda/guarda e 21 de referência; replays antes/depois
  distinguem erro numérico da proposta fora de domínio.
- `490936a`: singularidade/erro numérico propagados na busca de velocidades;
  CCFL exige interioridade, resíduo e máximo local. Busca mal resolvida retorna
  NUMERIC, sem produzir uma dobra espúria nem fallback cocorrente na singularidade.

Última checagem integrada: **3447 verificações Fortran em 12 programas e 19 testes
Python aprovados**, mais build release. `make test FC=gfortran LIBS='-llapack -lblas'`
e `python3 -m unittest discover -s tests -p 'test_validation_evidence.py' -v`.
Compiler WSL: GNU Fortran 13.3.0. Evidência final de código em
`validation/diagnostics/chexal_ccfl_status_20261118`; a evidência anterior
3331/10 permanece em `chexal_integration_20261118`. CI do commit 490936a:
SUCCESS, run 44390781711. Verificação numérica **não é validação experimental**.
O PR #2 permanece draft; a branch remota continua `audit/chexal-continuation-20261118`.

## Contratos físicos e limites que não podem desaparecer

Leia `CHEXAL_1992_IMPLEMENTATION.md` (fontes/equações),
`CHEXAL_1992_INTEGRATION.md` (operadores/unidades/status) e
`CHEXAL_REMAINING_LIMITS.md` (continuação local, busca finita e CCFL mal resolvida).

Chexal1992 p. 262 admite contracorrente **vertical**, gás ascendente e líquido
descendente. P.264 restringe horizontal a cocorrente. P.268 distingue C3 HIGH
e LOW e exige informação de seleção. Não usar variantes1997 ou autorizar
contracorrente inclinada: a catenária D não tem seção vertical.

O residual é `alpha*(C0*(jl+jg)+Vgj)-jg`, sem divisão pelo fluxo total.
AUTO detecta alternativas CC verticais; IDs 6/7 declaram ramo, e a recuperação
revalida a raiz contra o ramo alfa. Enumeração finita não prova inexistência de
pares de raízes/tangências não amostrados. CCFL retorna status de dobra ou limite
excedido; erros numéricos não são convertidos em NO_ROOT.

Derivadas preservam a assinatura constitutiva e recusam quinas. Os testes
exercitam j=0, orientação invertida, escala, U–primitivas, derivadas no mesmo
ramo, conservação e espectros reais; nenhum espectro foi reparado. Inventário
ausente tem recuperação e fluxo exatos. Fluxo nulo isolado exige informação
de inventário para selecionar alfa0/1; o Jacobiano bifásico completo no endpoint
de fase ausente é recusado, sem modo artificial.

**MATCHED/Roe com gás negativo continua sendo um impedimento separado.**
Replay de cinco propostas Newton preserva as 48 células admissíveis, mas falha
na junção; 15 replays dos três seletores Chexal recusam o traço horizontal.
Não foi aplicada barreira/projeção de gás nem ampliação de domínio. Consulte
`validation/diagnostics/chexal_junction_20261118` e os diagnósticos históricos.

SS3: dissertação Tabela 6.7 p. 136 informa 26 s; Baliño 2010 Tabela 3 p. 653, caso 11,
informa 76 s. A divergência permanece explícita, sem ajuste de parâmetros.
OSC/SS3 podem manter penetração gasosa: bloqueio completo não é requisito
universal, e oscilações isoladas de pressão não estabelecem slugging validado.

## Campanhas por versão

A campanha histórica de 05ad799 permanece preservada. As duas repetições antes
interrompidas completaram 2400 s no worktree `PRIMS-F-baseline-20261118`, sem
alpha-tail/Chexal: N201+201/dtmax 0,125 e N51+51/dtmax 0,03125, zero passos
rejeitados. Períodos 23,345470 e 23,344185 s. Dez intervalos finais acomodados
nos 13 canais; convergência das amplitudes ainda não demonstrada: +9,725% no último refino h e
+11,505% no último refino dt. `validation/baseline_repeat_20261118` guarda
métricas, hashes e manifesto do ZIP local. Históricos usam compilador 14.2;
repetições 13.3, diferença explicitada. Não há restart por checkpoint validado.

A primeira comparação usa o freeze `7321e67` em `PRIMS-F-comparison-20261118`:
13 tentativas, 6 COMPLETED e 7 FAILED, nenhum timeout. Inclui ambos os controles
estáveis até 2400, OSC Bhagwat até 2400 com refinos h/dt e extensão grossa até 3600.
Na extensão, dez ciclos finais: T=23,108641 s; semi-amplitude 9909,563 Pa;
variação 0,006414%. Há dependência importante de malha e passo. Falhas Chexal
SS3/OSC nessa versão incluíam piso numérico da derivada, documentado e corrigido
em b0d33cb; não atribuir esses dados ao código corrigido. Relatório e ZIP:
`validation/comparisons/frozen7321e67_20261118`.

A comparação homogênea final foi concluída no worktree separado
`PRIMS-F-comparison-final-20261118`, detached `490936a`: 13 execuções, seis completas,
sete falhas, nenhum timeout. Ambos os controles estáveis chegam a 2400 s;
OSC Bhagwat conclui 2400 com refinos separados e 3600 na extensão. Nos dez ciclos
finais da extensão, T=23,108641 s, semi-amplitude=9909,563 Pa;
todos os 13 canais passam o critério 1%, mas a convergência da amplitude com h/dt não foi demonstrada.
SS1 falha no traço gasoso negativo MATCHED. SS3 Bhagwat retorna NO_ROOT local;
Chexal SS3/OSC/refinos recusam propostas fora do domínio inclinado, já sem o
piso artificial de derivada. Não interpretar propostas recusadas como reversão
física aceita. Fontes e parâmetros permaneceram congelados; consulte
`validation/comparisons/frozen490936a_20261118` e o relatório consolidado
`docs/PRIMSF_CONTINUATION_VALIDATION_20261118.md`.

## Limites e trabalho científico restante

As análises separadas de malha/passo, dez intervalos, 13 canais e balanços estão
arquivadas. A convergência de amplitude e a validação experimental continuam
abertas. Raízes/CCFL usam busca finita; continuação local não define transição
física através da dobra, e o Jacobiano bifásico em fase ausente permanece recusado.
Falhas de domínio na catenária não autorizam estender Chexal vertical nem
ajustar entradas para forçar slugging. Manter PR draft até qualificação física.

Os PDFs e imagens de páginas permanecem fora do Git. Fontes conferidas, hashes,
páginas e a desambiguação pontual de C10 pelo manual primário RELAP 1995 estão
na nota de implementação. O candidato alpha-tail histórico pode ser reproduzido
no worktree `PRIMS-F-roe-review-20261118`/commit47433f8; seu runner exige o hash
exato anterior/candidato e não deve ser aplicado à árvore já modificada por Chexal.
