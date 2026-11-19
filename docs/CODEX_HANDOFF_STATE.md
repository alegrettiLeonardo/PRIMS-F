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
