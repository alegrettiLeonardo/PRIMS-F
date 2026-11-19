# Burr Eq. 15: fonte, conven??es e limites

Confer?ncia em 2026-11-19. Esta nota fixa a leitura da equa??o dimensional
solicitada; n?o declara valida??o experimental nem corrige silenciosamente a fonte.
Nenhum PDF, imagem de p?gina ou extra??o integral deve ser inclu?do no Git.

## Fontes conferidas

- Karl Peter Burr, *Modelagem Unidimensional de Sistema Oleoduto-Riser*,
  documento local sem data impressa, 15 p?ginas. Arquivo recebido:
  `C:/Users/leonardoab/Downloads/sistema_oleoduto_riser.pdf`, 210856 bytes,
  SHA256 `9f115779783f87a58719e730f1c3a35c331d221361bff1c78ffed9306f7feb1d`.
  P?ginas PDF e impressas coincidem. Eqs. 15/p.3 e 54/p.9 conferidas tamb?m
  visualmente, para distinguir subscritos de erros de extra??o.
- Y. Taitel, N. Lee e A. E. Dukler, *Gas-liquid transient flow in horizontal
  pipes: a model for predicting flow-pattern transitions*, relat?rio final,
  setembro de 1977, DOE/ID/01571--T1. [Texto prim?rio completo no OSTI](https://www.osti.gov/servlets/purl/5902200).
  SHA256 local `c953722d3abb7e22369b4c89b6fd15c6a28482d896a3bbdfc8f654f41f9854fa`.
  62 p?ginas PDF; p.6 impressa = PDF10, p.8?9 = PDF12?13, p.23 = PDF27.
- Taitel e Dukler (1976), *A model for predicting flow regime transitions in
  horizontal and near horizontal gas-liquid flow*, AIChE J. 22, 47?55,
  [DOI e metadados oficiais](https://doi.org/10.1002/aic.690220105).
  O texto integral desse artigo n?o foi obtido nesta confer?ncia; as equa??es
  abaixo atribu?das ao relat?rio de 1977 foram lidas no pr?prio relat?rio,
  que resume explicitamente o trabalho de 1976 nas pp.3?6.

O arquivo antigo `multiPhaseFlow/figs/sistema_oleoduto_riser.txt` n?o ? uma
transcri??o confi?vel para substituir este PDF: tem equa??es e numera??o
incompat?veis. N?o foi usado para fixar a Eq.15.

## Leitura literal e discrep?ncia interna

A Eq.15, p.3, imprime a transi??o estratificado ? n?o estratificado quando

```text
ug >= Ucrit15
Ucrit15 = (1 - hL/D) * sqrt[(rhol-rhog)*g*cos(beta)*A*alpha_p/(rhol*Sl)]
```

O denominador ? **rhol**, confirmado na imagem. `alpha_p` ? a fra??o de g?s,
pois a Eq.6/p.2 define `jg = alpha_p*ug`; `ug` ? a velocidade de fase, n?o `jg`.
A passagem imediatamente posterior define `Sl` como per?metro l?quido e `hL`
como altura l?quida na se??o. N?o identifica `Sl` como largura da interface.
As Eqs.8/p.2 e 12?13/p.3 separam parede l?quida `gamma*pi*D` e interface
`gamma_i*pi*D`; por essa geometria, a interpreta??o consistente ? o per?metro
molhado da parede, `Sl=gamma*pi*D`. Essa identidade ? uma interpreta??o geom?trica
expl?cita, n?o uma igualdade adicional escrita ao lado da Eq.15.

H? uma inconsist?ncia material no pr?prio documento. As Eqs.53?57, pp.8?9,
com `hL*=hL/D`, `ug*=A*ug/Ql0`, `Sl*=Sl/D` e
`F=(Ql0/A)*sqrt[rhog/(rhol-rhog)]*sqrt[D/(g*A*cos(beta))]`, d?o algebricamente:

```text
Ucrit53 = (1 - hL/D) * sqrt[(rhol-rhog)*g*cos(beta)*A*alpha_p/(rhog*Sl)]
```

Logo `Ucrit53/Ucrit15=sqrt(rhol/rhog)`. A forma adimensional n?o ? equivalente
? Eq.15 impressa no PDF local. Implementar a Eq.15 literal exige congelar essa
escolha e manter a diverg?ncia vis?vel; substituir `rhol` por `rhog` muda o modelo.

## Geometria, ?ngulo e sinais

Para se??o circular e interface plana, seja `phi=pi*gamma`, `0<phi<pi`:

```text
hL/D = (1-cos(phi))/2
alpha_p = 1 - phi/pi + sin(2*phi)/(2*pi)
Sl = D*phi                       [parede molhada pelo l?quido]
Si = D*sin(phi) = dAl/dhL         [largura da interface plana]
A = pi*D**2/4
```

Essas identidades s?o geometria circular exata. A rela??o inversa `hL(alpha_p)`
n?o ? fornecida no PDF de Burr. O relat?rio prim?rio de 1977 distingue `Sl` e
`Si` no texto/p.8 e Eqs.6?7; sua Eq.11/p.9 confirma `dAl/dhL=Si`.
N?o substituir `Sl` por `Si`, nem pelo per?metro total `Sl+Si`, sem mudar a lei.

Burr chama `beta` de inclina??o do oleoduto (p.2); sua Eq.3 usa `+g*sin(beta)`
na dire??o positiva do oleoduto. Isso ? consistente com `beta` positivo para
inclina??o descendente, em contraste com o ?ngulo horizontal positivo ascendente
que acompanha `-g*sin(theta)` no riser. Para o discriminante, `cos(beta)=cos(theta)`
quando `beta=-theta`. N?o usar graus desde a vertical de Chexal nesse termo.

A Eq.15 usa `ug` assinado, sem m?dulo. Os m?dulos nos Reynolds/Eqs.10?11 n?o
autorizam inserir `abs(ug)` na transi??o. O texto n?o documenta contracorrente,
faixa angular num?rica, mudan?a de sentido do eixo, nem limite monof?sico desse
crit?rio. `0<alpha_p<1`, `rhol>rhog>0`, `D>0`, `g>0`, `Sl>0` e `cos(beta)>0`
s?o condi??es de avalia??o da forma estratificada adotada, n?o uma prova de
validade f?sica em qualquer inclina??o que as satisfa?a. Fases ausentes devem
seguir o contrato monof?sico pr?prio, sem extrapolar a Eq.15 ou fabricar alturas.
Crit?rio falso com g?s negativo tamb?m n?o prova que o fechamento estratificado
seja admiss?vel para esse estado.

## N?o estratificado e retorno

Nas pp.2?3, o oleoduto tem duas leis: Eq.8, equil?brio estratificado local sem
in?rcia, e Eq.14, drift-flux n?o estratificado, atribu?do no texto a Swananda e
Gajar. A Eq.15 n?o distingue bolhas, golfadas ou anular. O riser ? tratado como
sempre n?o estratificado, com Bendiksen na Eq.23/p.5. Usar outro fechamento ? uma
op??o expl?cita da implementa??o; n?o ? uma consequ?ncia da Eq.15.

A busca no PDF completo n?o encontrou regra n?o estratificado ? estratificado,
banda de histerese, tempo de relaxa??o, estado de mem?ria nem crit?rio espec?fico
de revers?o de fase. Aplicar a desigualdade complementar ao estado atual para
retornar ? uma pol?tica alg?brica adicional; n?o deve ser atribu?da a Burr como
lei de transi??o reversa validada. Tamb?m ? necess?rio esclarecer qual estado
ou solu??o candidata fornece `alpha_p`, `hL` e `ug` quando as duas leis admitem
estados diferentes para as mesmas vari?veis conservadas.

O relat?rio prim?rio de 1977 fornece contexto, mas n?o fecha essa lacuna:
Eq.1/p.6 usa `rhog*A_l_prime`, com `A_l_prime=Si`, e portanto difere de Burr15
em densidade **e** geometria. Trata tubo horizontal, com evolu??o da altura
l?quida e hip?teses de g?s mais r?pido que o l?quido (p.9); pp.3?4 distinguem
mudan?as quase est?ticas e r?pidas. Na p.23, o caminho
r?pido reverso cruza transi??es em vaz?es distintas do caminho direto.
Isso n?o fornece coeficientes de histerese ou retorno alg?brico transfer?veis ?
caten?ria PRIMS-F. N?o foi localizada uma regra reversa prim?ria diretamente
aplic?vel ? Eq.15 literal. O modelo revers?vel f?sico permanece sem qualifica??o.

## Refer?ncias num?ricas independentes

C?lculo por geometria circular fechada, com aritm?tica `mpmath` de 60 d?gitos;
sem chamar a implementa??o Fortran. Propriedades SI: `rhol=998`, `rhog=2.3`,
`g=9.81`, `D=0.05254`, `beta=0`. Valores abaixo s?o velocidades de fase em m/s.

| `gamma` | `alpha_p` | Eq.15 literal | Eqs.53?57 | Taitel1977 Eq.1 |
|---|---:|---:|---:|---:|
| 1/2 | 0.5 | 0.179274503957725615 | 3.73439181205966400 | 4.68036605232962554 |
| 1/3 | 0.804498890522114679 | 0.417766046911252909 | 8.70230886434247842 | 9.56937307448993684 |
| 2/3 | 0.195501109477885321 | 0.0485410064280379163 | 1.01113729477531016 | 1.57244013476571494 |

No caso meio cheio, `ug=0.2` separa Eq.15 das outras duas; `ug=4` separa Burr53
de Taitel1977. `D` multiplicado por quatro duplica cada limite. Para `beta=+2`
ou `-2` graus, Burr15 d? `0.179219891047293187` m/s; isso verifica a convers?o
angular e a simetria do cosseno, n?o amplia o dom?nio experimental.
A escala comum de ambas as densidades cancela, e converter as mesmas unidades
consistentes mant?m o resultado SI. Casos de igualdade, aus?ncia de fase,
?ngulo fora da geometria estratificada e g?s negativo devem verificar o contrato
de status/sele??o, sem `abs`, clipping ou troca autom?tica de correla??o.

Esta ? verifica??o independente da f?rmula e das conven??es. N?o ? valida??o
experimental da transi??o nem demonstra??o da adequa??o da Eq.15 em transientes.

## Confer?ncia complementar: disserta??o e Bali?o2010

- `Dissertacao_Leonardo.pdf`, vers?o local de agosto/2025, 160 p?ginas PDF,
  SHA256 `dcc6c069d89099950fe737112b331a678727cc3c6929deb4f105072f0f0e8e81`.
  P.61/PDF63 explicita `Sl=gamma*pi*D` e `Si=gamma_i*pi*D` ap?s a Eq.3.14.
  Ap?ndice E, pp.157?158/PDF159?160, Eqs.12.1?12.4 e 12.8, deduz a ?rea do
  segmento circular e `dAl/dh=Si`. H? discrep?ncia interna: Eq.3.16/p.62 imprime
  `sin(pi*gamma)`, enquanto Eq.12.3/p.157 deduz corretamente `sin(2*pi*gamma)`.
  A refer?ncia independente de altura usa esta ?ltima geometria, verific?vel
  diretamente pela ?rea circular, sem copiar a express?o inconsistente.
  O resumo/PDF15 e o Ap?ndice E/p.157 assumem oleoduto estratificado; as se??es
  3.3.1?3.3.2, pp.61?62, apresentam leis distintas de fechamento. Nessas passagens
  conferidas n?o h? regra para retorno n?o estratificado ? estratificado.
- Bali?o et al. (2010), *Modeling and simulation of severe slugging in air?water
  pipeline?riser systems*, [DOI](https://doi.org/10.1016/j.ijmultiphaseflow.2010.04.003),
  SHA256 local `ce04a2d1b79d3a10d442fbaeb3863e4639c8508bd34ba9e1a414ae03e5479e2b`.
  P.645/PDF3 assume escoamento estratificado no oleoduto e despreza in?rcia no
  riser (NPW). A Eq.7 ? o balan?o estratificado avaliado no estado estacion?rio;
  distingue os per?metros l?quido, g?s e interface. Essa formula??o n?o fornece
  um discriminante din?mico Burr15 nem uma lei reversa de troca de fechamento.

Essas confer?ncias refor?am a conven??o geom?trica e delimitam o papel dos
modelos consultados. N?o constituem uma afirma??o de inexist?ncia universal de
modelos de retorno na literatura.
