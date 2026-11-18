# Bhagwat conforme Azevedo, Baliño e Burr (2015)

Atualização: 18/11/2026. **`BHAGWAT=3` está implementado** a partir da transcrição completa das eqs.116–128 em Azevedo, Baliño e Burr (2015), p.121. A identificação precisa é **“Bhagwat–Ghajar conforme Azevedo, Baliño e Burr (2015), tubo circular, domínio cocorrente”**. A conferência integral do artigo original de 2014 continua pendente; a implementação não afirma essa equivalência documental mais forte.

## Fonte conferida e integridade

**G. R. Azevedo, J. L. Baliño e K. P. Burr. _Linear stability analysis for severe slugging: sensitivity to void fraction correlations_. IMAACA 2015, pp.115–124. ISBN 978-88-97999-63-8.** É pesquisa original de aplicação e comparação de fechamentos, com a transcrição utilizada na seção 6.3; não é a publicação original da correlação de 2014.

- [Anais no endereço institucional do evento](https://www.msc-les.org/proceedings/imaaca/2015/IMAACA2015.pdf).
- Localização: p.121 impressa, página 135 do PDF de 166 páginas.
- SHA-256: `51891ab67e550804e732603faf37bfac53504ee2d1e36b4f589a4921aefbd171`.
- Página 121 impressa / 135 do PDF conferida visualmente; consultar a seção 6.3 no link da fonte. O PDF e o recorte de terceiros não são redistribuídos nesta revisão.
- [Cálculo independente](../scripts/generate_bhagwat_fixtures.py) e [valores de referência](../tests/fixtures/bhagwat_azevedo2015.json).

O original identificado é **S. M. Bhagwat e A. J. Ghajar, _A flow pattern independent drift flux model based void fraction correlation for a wide range of gas–liquid two phase flow_, International Journal of Multiphase Flow 59 (2014), 186–205**, [DOI 10.1016/j.ijmultiphaseflow.2013.11.001](https://doi.org/10.1016/j.ijmultiphaseflow.2013.11.001). Seus metadados foram confrontados no [registro da editora](https://www.sciencedirect.com/science/article/abs/pii/S0301932213001705) e no [catálogo institucional dos autores](https://ceat.okstate.edu/mae/research/single-and-two-phase-flow-studies/publications). A busca inicial recuperou apenas resumo e metadados desse original; a obtenção posterior do PDF de 2015 permitiu implementar uma versão integralmente identificada.

## Especificação matemática

Todas as grandezas abaixo são dimensionais em SI. Para seção circular, $D_h=D$. Com $j_l=(1-\alpha)u_l$, $j_g=\alpha u_g$ e $j=j_l+j_g$,

$$
r=\rho_g/\rho_l,\quad Re_{tp}=\rho_ljD/\mu_l,\quad
x=\frac{\rho_gj_g}{\rho_gj_g+\rho_lj_l},\quad \beta=\frac{j_g}{j_l+j_g}.
$$

A eq.116, escrita de forma numericamente equivalente, é

$$
C_0=(2-r^2)(1-w)+(B+C_{0,1})w,\quad
w=\frac{Re_{tp}^2}{10^6+Re_{tp}^2},\quad
B=\left[\frac{1+r^2\cos\theta}{1+\cos\theta}\right]^{(1-\alpha)/5},
$$

$$
C_{0,1}=0.2(1-\sqrt r)\left[(2.6-\beta)^{0.15}-\sqrt{f_{tp}}\right](1-x)^{1.5}.
$$

O fator da eq.118 satisfaz

$$
\frac1{\sqrt{f_{tp}}}=-4\log_{10}\left(\frac{\varepsilon}{3.7D}
+\frac{1.256}{Re_{tp}\sqrt{f_{tp}}}\right).
$$

Esse $f_{tp}$ é um parâmetro interno da correlação, distinto do cálculo por ramos do atrito de parede. A implementação utiliza a equação acima literalmente; não insere $16/Re$ nem transição em 2000/2300 nessa expressão. Sua extensão formal a Reynolds baixos não é apresentada como verificação experimental adicional.

A velocidade de deriva das eqs.122–128 é

$$
U_d=(0.35\sin\theta+0.45\cos\theta)
\sqrt{gD(1-r)}\sqrt{1-\alpha}\,C_2C_3C_4,
$$

$$
C_2=\begin{cases}(0.434/\log_{10}\mu_l^*)^{0.15},&\mu_l^*>10,\\1,&\mu_l^*\le10,\end{cases}
\quad\mu_l^*=\frac{\mu_l}{0.001\,\mathrm{Pa\,s}},
$$

$$
La=\frac1D\sqrt{\frac{\sigma}{g(\rho_l-\rho_g)}},\qquad
C_3=\begin{cases}(La/0.025)^{0.9},&La<0.025,\\1,&La\ge0.025,\end{cases}
$$

$$
C_4=\begin{cases}-1,&-50^\circ\le\theta<0^\circ\ \text{e}\ Fr_{sg}\le0.1,\\
1,&\text{demais casos},\end{cases}
\quad Fr_{sg}=\sqrt{\frac{\rho_g}{\rho_l-\rho_g}}\frac{j_g}{\sqrt{gD\cos\theta}}.
$$

Em $\theta=90^\circ$, $C_4=1$ e o Froude descendente não precisa ser calculado. Para baixo Froude, o ramo descendente também produz salto de $C_4$ ao atravessar $\theta=0$. Uma derivada espacial na base do riser deve respeitar o domínio angular ascendente e usar o limite pela direita, sem amostrar uma correlação descendente diferente. A constante 0.434 e o logaritmo decimal de $C_2$ foram conferidos visualmente e preservados, mesmo que essa expressão apresente salto no limiar de viscosidade. Nenhum ajuste para tornar o modelo suave foi introduzido.

O residual usado por primitivas, estacionário e Jacobiana é

$$
\Phi=\frac{\alpha[-u_g+C_0j+U_d]}{w_u}.
$$

`CdUd_drift(...,closure_id=3)` e `CdUd_bhagwat` devolvem $C_0$ adimensional e $U_d$ em m/s. Entradas $\rho,u$ são restauradas a SI antes do cálculo. Os nomes legados contendo `swananda` encaminham para essa versão; Swanand é o primeiro nome de Bhagwat.

## Domínio qualificado e limites

Exigem-se $0<\alpha<1$, $\rho_l>\rho_g>0$, $D,g,\mu_l,\sigma>0$, $\varepsilon\ge0$, seção circular coerente, $-50^\circ\le\theta\le90^\circ$, $j_l,j_g\ge0$ e $j>0$. Na aplicação ao riser, $\theta\ge0$. A rugosidade deve admitir raiz positiva da eq.118: $\varepsilon/D\ge3.7$ é rejeitado.

Uma fase pode ter velocidade zero e massa positiva. Portanto, $u_l=0$ ou $u_g=0$ são estados bifásicos avaliáveis se $j>0$. Ausência de fase é um limite distinto, declarado pela EOS e pelos chamadores. Fluxo total nulo não recebe valor artificial.

Velocidades negativas produzem falha, sem módulo ou recorte. Em contracorrente, $x$ e $\beta$ deixam de representar frações usuais; por exemplo, $(1-x)^{1.5}$ pode deixar de ser real. A correção para inclinação negativa não estabelece validade para fases de sinais opostos. Uma extensão exige especificação física e testes próprios. Esse limite pode encerrar um transiente SS1/SS3; o diagnóstico deve permanecer visível.

## Solução numérica e derivadas

Com $y=1/\sqrt{f_{tp}}$, resolve-se

$$
H(y)=y+4\log_{10}(q)=0,\qquad
q=\frac{\varepsilon}{3.7D}+\frac{1.256y}{Re_{tp}}.
$$

$H'(y)>0$ para $y>0$. Um intervalo é construído por expansão e Newton permanece dentro dele; passos que sairiam tornam-se bisseções. Há limite de iterações e falha explícita. A derivada vem da diferenciação implícita de $H$, sem perturbação numérica da raiz.

`Dlei` deriva analiticamente $j,r,x,\beta,Re_{tp},f_{tp},C_0,U_d$ e aplica a cadeia das escalas. Dentro do ramo $La<0.025$, $dC_3/C_3=-0.45\,d(\rho_l-\rho_g)/(\rho_l-\rho_g)$; no outro ramo essa parcela é zero. $C_3$ é contínuo e pode ter mudança de derivada. O salto de $C_4$ em $Fr_{sg}=0.1$ retorna `CLOSURE_DERIVATIVE` quando se solicita a derivada.

Na inversão do estado, $C_0$ e $U_d$ não são constantes: cada tentativa de raiz deve recomputar a correlação e conferir o mesmo residual. `DCdUd_swananda` retorna derivadas dos coeficientes em relação às entradas adimensionais, mantendo $U_d$ dimensional; `Dlei` retorna derivadas do residual adimensional.

## Comparação com os códigos fornecidos

| Item | Python do ZIP, `closureLaw.py` | Fortran do GitHub original, `module_closure_law.f90` |
|---|---|---|
| $j_l,j_g,Re_{tp}$ | Equivalentes às eqs.119–121 após restauração de SI. | Fatores de densidade nas grandezas chamadas superficiais e outra expressão de Reynolds. |
| $x,\beta,La$ | Correspondem às definições no domínio cocorrente. | Definições e recortes diferentes das eqs.120,121,128. |
| $C_0,U_d$ | Próximos à transcrição, incluindo potências e coeficiente 0.45. | Fórmulas distintas da seção 6.3. |
| Descida | Acrescenta $C_{0,1}=0$ para alguns estados e recorta $x$. | Outras regularizações. |
| Colebrook | Ponto fixo sem limite de iterações no operador principal. | Ramos e recortes não correspondentes à eq.118 literal. |

A regra adicional do ZIP que zera $C_{0,1}$ em descida não aparece na eq.117 nem no texto adjacente conferido. Ela não foi importada. Recortes de qualidade, fluxo ou argumentos também não foram mantidos. A proximidade do Python em cocorrente é evidência de intenção; não certifica suas derivadas ou seus tratamentos de domínio.

## Verificação independente

```sh
python3 scripts/generate_bhagwat_fixtures.py
FC=gfortran LIBS='-llapack -lblas' python3 scripts/run_unit_tests.py
```

O primeiro comando exige `mpmath` apenas para regenerar os valores. Usa 80 dígitos e bisseção de Colebrook, independentemente do Newton do Fortran. Os sete casos incluem $C_2$ viscoso, $C_3$ capilar, ambos os sinais de $C_4$, inclinação vertical e fase estagnada.

O teste Fortran confere coeficientes, residual e três escalas; compara derivadas com diferenças finitas; fecha uma raiz de fração de gás e verifica falhas de domínio. Estacionário e primitivas têm testes separados de integração do ID3. Essas evidências qualificam a implementação matemática; não constituem validação experimental do período ou amplitude de slugging.
