# Candidato: busca de raiz abaixo da cauda original de alfa

**Estado: CANDIDATO REVISADO, NÃO APLICADO À PRODUÇÃO.** O patch não foi usado
nos resultados transientes da campanha. Este pacote contém somente a diferença,
testes, referências numéricas e evidências da revisão isolada; não contém cópias
integrais dos módulos, executáveis ou objetos.

## O que corrige

A busca original cobre `|logit(alpha)| <= 28`, cujo limite inferior é
aproximadamente `6,9144e-13`. Uma raiz positiva menor pode existir mesmo quando
essa busca retorna `STEADY_NO_ROOT`. O patch adiciona uma busca geométrica
seguida de bisseção relativa **somente depois de zero raízes detectadas**, para
Bhagwat ID3, `theta >= 0`, `jl >= 0`, `jg > 0`, com uma condição suficiente de
positividade de C0 e unicidade na cauda inferior. Não substitui uma raiz por um
piso e não converte pequeno fluxo em ausência de fase.

O domínio de contracorrente não é ampliado. Caminhos anteriormente bem-sucedidos,
múltiplas raízes e limites monofásicos permanecem intactos. O patch não resolve
por si só `ul = 0`, fluxo de junção negativo, trajetórias de Roe fora do domínio
ou a validação do transiente. Não é uma implementação de Chexal.

## Arquivos e evidência

- `alpha_tail.patch`: diferença para `scr/module_steady_state.f90`.
- `run_tests.py`: runner portátil; aplica a diferença apenas em cópias temporárias.
- `test_alpha_tail_guard.f90` e `tail_fixtures.dat/json`: 55 verificações e sete
  referências independentes de 80 dígitos, incluindo representabilidade/subnormais.
- `probe_alpha_tail_review.f90`: revisão independente de escalas, tolerâncias e hints.
- `REVIEW.md` e `prior_evidence.json`: revisão já realizada, com seus resultados e limites.
- `candidate_manifest.json`: hashes do conteúdo e dos fontes anterior/candidato.

Evidência anterior: O0/O3 tiveram 55 guards e 398 testes estacionários em cada
versão, com 101 registros anteriores idênticos byte a byte. A revisão O2 passou
129 contratos: 108 novas raízes, 18 sucessos e três exclusões preservados. Há seis
casos com tolerância frouxa `1e-6` em que uma limitação anterior de precisão
permanece deliberadamente inalterada; eles não são apresentados como aprovação
de precisão. A prova de unicidade é local à cauda qualificada, não global.

## Conferir e testar sem modificar o solver

Na raiz do repositório, com Python 3, Git e GNU Fortran disponíveis:

```bash
python3 validation/diagnostics/alpha_tail_candidate/run_tests.py --check-only
FC=gfortran python3 validation/diagnostics/alpha_tail_candidate/run_tests.py
FC=gfortran python3 validation/diagnostics/alpha_tail_candidate/run_tests.py --review
```

O primeiro comando verifica hashes e aplicação do patch sem compilar. O segundo
reproduz O0/O3; o terceiro também executa a revisão independente O2. Para uma
checagem mais curta, use `--optimization O0`. Os logs e `summary.json` são escritos
em `tmp/alpha_tail_candidate_tests`; `--output` permite outro diretório. Nenhum
comando do runner grava em `scr/` ou `bin/`. Não requer NumPy, mpmath ou LAPACK:
as referências de alta precisão já estão incluídas. Não habilitar fast-math,
flush-to-zero ou trap de underflow nestes testes de subnormais.

## Aplicar posteriormente em uma branch própria

Depois de revisar o diff e os resultados, na raiz do clone:

```bash
git switch -c fix/alpha-lower-tail
git apply --check validation/diagnostics/alpha_tail_candidate/alpha_tail.patch
git apply validation/diagnostics/alpha_tail_candidate/alpha_tail.patch
git diff --check
FC=gfortran python3 validation/diagnostics/alpha_tail_candidate/run_tests.py --review
```

O runner reconhece tanto o hash anterior quanto o candidato exato e reconstruirá
as duas versões temporárias para comparação. Se o fonte tiver outra alteração,
ele exige rebase/revisão explícitos. Depois da integração, execute a suíte geral
do repositório e novos transientes identificados pelo novo hash; preserve a
campanha anterior e seus motivos de interrupção. Aprovação local do helper não
constitui confirmação de que os casos SS1/SS3/OSC agora completam ou reproduzem slugging.
