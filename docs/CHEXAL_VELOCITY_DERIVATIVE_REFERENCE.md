# Referência independente: derivada perto de ul=0 positivo

Esta nota documenta o defeito numérico observado após congelar `7321e67` e a
verificação independente da correção do adaptador de derivadas. O núcleo
Chexal, seus coeficientes, domínio físico e parâmetros da campanha não mudam.

## Referência e derivação

Mantêm-se alpha=0,4; rho_l=998 kg/m³; rho_g=2,3 kg/m³; ug=0,5 m/s;
theta=1,151159 rad desde a horizontal; D=0,05254 m; g=9,81 m/s²;
mu_l=0,001 Pa.s; mu_g=0,000018 Pa.s; sigma=0,072 N/m.
Calcula-se a derivada parcial de
`R=alpha*(C0*((1-alpha)*ul+alpha*ug)+Vgj)-alpha*ug`
com respeito a ul SI, mantendo as outras grandezas constantes.

As Eqs.3-2/3-3, 3-15, 3-21, 3-24 e 3-34 de Chexal1992
(pp.261–269) dão, nessa peça positiva:

- Reg=1342,68888888888889 domina Rel; portanto B1 e C0 não dependem de ul.
- C3=`2*exp(-Rel/60000)` está acima do piso 0,5.
- A dependência de Vgj em ul é apenas a exponencial de C3.

Logo, a expressão analítica local independente é

`dR/dul = alpha*(1-alpha)*(C0 - Vgj*rho_l*D/(mu_l*60000))`.

`scripts/generate_chexal_derivative_reference.py` avalia a correlação em
mpmath com 60 dígitos e confronta diferenciação arbitrariamente precisa com
essa expressão analítica. Os resultados concordaram com erro absoluto menor
que 4e-62. `tests/fixtures/chexal_velocity_derivatives.json` preserva os números,
parâmetros e versão mpmath 1.2.1 utilizada. Não foi chamado o código Fortran
para gerar esses valores.

| ul SI (m/s) | dR/dul de referência |
|---|---|
| 2,2586e-8 | 0,1836908994576642193379204514612842432495967349039477632 |
| 1e-14 | 0,1836908982971393978002752516130279724507857393470429714 |

A próxima mudança de expressão positiva, Rel=Reg, ocorre em
ul=0,0426779484895717361 m/s; o piso de C3 ocorre em 2,64383803984041669 m/s.
Portanto não há quina constitutiva local nos dois pontos da tabela. ul=0 é,
contudo, a fronteira do domínio cocorrente para essa inclinação. As avaliações
da tabela são pontos de teste do residual e não satisfazem R=0; elas não são
estados transientes medidos ou uma validação experimental.

## Testes executados

GNU Fortran 13.3.0, `-std=f2018 -O0 -fcheck=all
-ffpe-trap=invalid,zero,overflow`, em diretórios temporários separados.
O programa `tests/test_chexal_velocity_reference.f90` passou **21 verificações**
contra o adaptador corrigido: as duas referências em três escalas; recusa da
derivada bilateral em ul=0; recusa de contracorrente inclinada e recusa da
quina real Rel=Reg. As escalas `(wu,wr)` foram `(1,1)`, `(3,1000)` e `(0.2,10)`.
A derivada `d(R/wu)/d(ul/wu)` deve conservar exatamente a mesma referência SI.

O mesmo programa foi compilado contra uma cópia temporária do adaptador de
`7321e67`, sem modificar a baseline, e falhou. Um probe separado registrou:

| ul SI | Antes, 7321e67 | Depois, adaptador corrigido |
|---|---|---|
| 2,2586e-8 | status 0; derivada 0,18369090253895967; erro 3,0813e-9 | status 0; derivada 0,18369089945825487; erro 5,9064e-13 |
| 1e-14 | status 13, derivada indefinida | status 0; derivada 0,18369089829680174; erro 3,3765e-13 |

Os três pares de escalas deram os mesmos resultados. O decimal arredondado
2,2586e-8 não reproduz exatamente o limiar do estado da campanha: aqui a versão
antiga retornou status 0, porém com perda de precisão. Para 1e-14, o falso bloqueio
por status 13 foi reproduzido diretamente. Os CSVs antes/depois e hashes ficam
em `validation/diagnostics/chexal_velocity_reference_20261118/`.

## Avaliação da correção

Um stencil unilateral de quarta ordem, apontando para o interior do domínio,
é consistente nesses pontos não nulos. A checagem adaptativa da assinatura em ambos os
lados próximos ao estado, seguida de assinaturas iguais em todos os pontos
do stencil e convergência entre h e h/2, evita usar esse recurso para esconder
as quinas exercitadas. O teste independente confirma a derivada; a simples
concordância entre duas diferenças finitas não seria referência suficiente.

Essa verificação é local. Assinaturas em uma malha finita não provam ausência
de toda mudança de peça entre pontos arbitrários; o controle de convergência
e as recusas de domínio continuam necessários. A correção não autoriza
contracorrente na catenária, não transforma ul exatamente nulo em uma derivada
bilateral válida e não garante que a campanha transiente conclua ou reproduza
slugging físico.


## Guarda perto de uma mudança de peça

A revisão independente identificou que um probe bilateral fixo de 12,5% da
velocidade também recusava pontos suaves próximos a Rel=Reg. Com os mesmos
parâmetros acima, mas ug=0,005 m/s, o ponto ul=0,95*ulcrit é suave:
ulcrit=4,267794848957174e-4 m/s e ul=4,054405106509315e-4 m/s.
A guarda passa a reduzir seu raio até encontrar a mesma assinatura de ambos
os lados; não aceita pontos que deixem de ser representáveis nem uma quina
persistente. A regressão cobre esse ponto e a quina exata. O status13 continua
sendo uma recusa de derivada numérica e exige diagnóstico: não prova, sozinho,
que um estado tenha ultrapassado o domínio físico.
