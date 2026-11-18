# Estado para continuação no Codex

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
- `3096c3b`: núcleo Chexal1992 SI, ar–água, Reynolds assinados, CCFL e ramos,
  referências independentes mpmath60. Núcleo:139 checks diagnósticos.
- `40a7d44`: integração de residual, primitivas, derivadas/Jacobianos,
  estacionário, contorno e junção; `868eeb9`: fixture adicional do ramo baixo.
  `closure_riser=5` é CHEXAL_1992, `6` HIGH, `7` LOW; Bhagwat permanece3.
- Análise: seleção explícita de coortes OSC e bloqueio de mistura de fontes,
  fechamento/EOS/parâmetros. NumPy1/2 suportados no integrador trapezoidal.

Última checagem integrada: **3331 verificações Fortran em10 programas e19 testes
Python aprovados**, mais build release. `make test FC=gfortran LIBS='-llapack -lblas'`
e `python3 -m unittest discover -s tests -p 'test_validation_evidence.py' -v`.
Compiler WSL: GNU Fortran13.3.0. Logs, flags e hashes estão em
`validation/diagnostics/chexal_integration_20261118`; essas aprovações são
verificação numérica, **não validação experimental**.

## Contratos físicos e limites que não podem desaparecer

Leia `CHEXAL_1992_IMPLEMENTATION.md` (fontes/equações) e
`CHEXAL_1992_INTEGRATION.md` (operadores/unidades/status).

Chexal1992 p.262 admite contracorrente **vertical**, gás ascendente e líquido
descendente. P.264 restringe horizontal a cocorrente. P.268 distingue C3 HIGH
e LOW e exige informação de seleção. Não usar variantes1997 ou autorizar
contracorrente inclinada: a catenária D não tem seção vertical.

O residual é `alpha*(C0*(jl+jg)+Vgj)-jg`, sem divisão pelo fluxo total.
AUTO detecta alternativas CC verticais; IDs6/7 declaram ramo, e a recuperação
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
Replay de cinco propostas Newton preserva as48 células admissíveis, mas falha
na junção;15 replays dos três seletores Chexal recusam o traço horizontal.
Não foi aplicada barreira/projeção de gás nem ampliação de domínio. Consulte
`validation/diagnostics/chexal_junction_20261118` e os diagnósticos históricos.

SS3: dissertação Tabela6.7 p.136 informa26s; Baliño2010 Tabela3 p.653, caso11,
informa76s. A divergência permanece explícita, sem ajuste de parâmetros.
OSC/SS3 podem manter penetração gasosa: bloqueio completo não é requisito
universal, e oscilações isoladas de pressão não estabelecem slugging validado.

## Campanhas por versão

A campanha histórica de05ad799 permanece preservada:28 configurações,13 completas,
13 falhas do solver e2 interrupções de sessão, além de2 tentativas de timeout.
Os refinamentos interrompidos foram reiniciados do zero no worktree separado
`PRIMS-F-baseline-20261118`, detached05ad799, sem alpha-tail/Chexal:
201+201/dtmax0,125s e51+51/dtmax0,03125s, ambos até2400s.
Label `repeat_20261118_frozen05ad799`. Não há restart por checkpoint validado;
não concatenar séries reiniciadas. Compilador histórico14.2 e atual13.3 são
identificados como diferença de reprodução.

A comparação nova usa worktree congelado `PRIMS-F-comparison-20261118`,
Bhagwat3 e Chexal5, BAROTROPIC/MATCHED/LAGGED/ROE,51+51,dtmax0,125s,2400s.
O script `scripts/run_closure_comparison.py` preserva EOS, vazões e geometria
e grava plano, hashes, comandos e status. Consulte manifests e resumos de
`validation/comparisons/` para o resultado de cada execução; não confunda
RUNNING/TIMEOUT/FAILED com COMPLETED.

## Próximo trabalho condicionado à evidência

Avaliar separadamente refino de malha e de passo, dez intervalos finais,
pressão, holdup, fluxos por fase, inventários e balanço integral. Examinar
acumulação/compressão, produção, penetração e blowdown conforme o regime.
Falhas de domínio na catenária não autorizam estender Chexal vertical nem
ajustar entradas para forçar slugging. Manter PR draft até qualificação física.

Os PDFs e imagens de páginas permanecem fora do Git. Fontes conferidas, hashes,
páginas e a desambiguação pontual de C10 pelo manual primário RELAP1995 estão
na nota de implementação. O candidato alpha-tail histórico pode ser reproduzido
no worktree `PRIMS-F-roe-review-20261118`/commit47433f8; seu runner exige o hash
exato anterior/candidato e não deve ser aplicado à árvore já modificada porChexal.
