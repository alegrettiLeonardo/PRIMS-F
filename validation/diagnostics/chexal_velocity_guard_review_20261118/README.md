# Revisão independente da guarda de derivadas — 18/11/2026

Este pacote guarda seis arquivos curtos de revisão. São testes locais de um
operador; não são novas execuções transientes nem validação experimental.

| Arquivo | Versão efetivamente exercitada | Resultado |
|---|---|---|
| `probe_before.log` | Adaptador congelado `7321e67` | Ponto suave status 0; quina exata status 13. |
| `probe_fixed_guard_regression.log` | Stencil unilateral intermediário, não commitado, com probe bilateral fixo de 12,5% | Ponto suave também retornava 13: regressão detectada antes do commit final. Não atribuir este resultado a `b0d33cb`. |
| `probe_after.log` | Guarda bilateral adaptativa incorporada em `b0d33cb` | Recupera status 0 no ponto suave; mantém 13 na quina exata. Coincide byte a byte com o probe anterior nesses dois casos. |
| `tail_test.log` | Fonte final do adaptador e teste ampliado | 87 verificações PASS, incluindo três escalas, expressão analítica no ponto suave, zero exato e quinas. |
| `reference_test.log` | Fonte final, referências independentes com 60 dígitos | 21 verificações PASS. |
| `review_hashes.txt` | Hashes preservados durante a revisão | Identifica fontes anterior/final, teste e logs. Os caminhos relativos apontam para o diretório temporário original. |

Ponto do probe: alpha=0,4, rho_l=998 kg/m³, rho_g=2,3 kg/m³, ug=0,005 m/s,
mu_l=0,001 Pa·s, mu_g=0,000018 Pa·s, D=0,05254 m, g=9,81 m/s²,
sigma=0,072 N/m, theta=1,151159 rad desde a horizontal.
Rel=Reg em `ulcrit=4,267794848957174e-4 m/s`; o ponto suave é
`ul=0,95*ulcrit=4,054405106509315e-4 m/s`. No ponto suave, as assinaturas
nos probes ul±1% e no centro coincidem (270369). Na quina exata, as assinaturas
dos dois lados diferem (270369/262176).

A guarda final reduz o raio bilateral até ambos os probes terem a assinatura
do estado, com limite de 16 tentativas. Recusa pontos não representáveis e
mudanças persistentes. O stencil unilateral e sua extrapolação continuam
verificando assinaturas e convergência; nenhuma variável física é recortada.
Essa verificação finita não constitui prova global de ausência de mudanças
de trecho. Status 13 é uma recusa de derivada e, sozinho, não prova violação
do domínio físico.

Compilação executada com GNU Fortran 13.3.0, `-std=f2018 -O0 -fcheck=all
-ffpe-trap=invalid,zero,overflow -ffree-line-length-none`, em
`/home/leonardoab/Documentos/PRIMS-F-chexal/tmp/derivative_guard_review`.
O probe foi compilado separadamente contra os três adaptadores; não se alterou
a baseline. O primeiro comando de link dos testes adicionais procurou objetos
temporários ausentes e falhou; após compilar explicitamente os módulos com `-c`,
os dois programas passaram. Não houve correção física motivada por esse erro
de comando.

Hash final de `module_closure_law.f90`:
`2c7e9b7d81651faa6580be044fcf04e8430fc54670afd8b10dd614d2e955d233`.
As referências analíticas e os testes reproduzíveis estão em
[CHEXAL_VELOCITY_DERIVATIVE_REFERENCE.md](../../../docs/CHEXAL_VELOCITY_DERIVATIVE_REFERENCE.md)
e [test_chexal_derivative_tail.f90](../../../tests/test_chexal_derivative_tail.f90).
