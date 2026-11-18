# Física implementada e critérios de verificação

Este documento descreve as equações efetivamente usadas nesta revisão do PRIMS-F. A base conservativa segue Burr [B] e a extensão barotrópica da dissertação de Leonardo Alegretti Belaroza [D]. O fechamento do riser é uma escolha explícita: a variante de Bendiksen publicada por Baliño, Burr e Nemoto [BNB] tem código `2`; Bhagwat–Ghajar, conforme a transcrição completa de Azevedo, Baliño e Burr (2015) [ABB], tem código `3`. A conferência integral do artigo original de Bhagwat–Ghajar (2014) continua pendente; o código `3` identifica precisamente a versão e o domínio verificados em [ABB].

Uma execução com os parâmetros da dissertação e `closure_riser=2` deve ser descrita como **parâmetros da dissertação com fechamento Bendiksen–Baliño**; com `closure_riser=3`, como **Bhagwat conforme Azevedo–Baliño–Burr (2015), domínio cocorrente**. A aprovação de testes unitários verifica operações matemáticas e contratos do programa; a reprodução quantitativa de severe slugging exige ensaios transientes, balanços, refinamentos e comparação com dados experimentais compatíveis.

## 1. Coordenada, variáveis e equações de conservação

A coordenada $s$ cresce da entrada do oleoduto para o separador. A inclinação $\theta$ é medida em relação à horizontal nesse sentido: $\theta=-\beta<0$ no oleoduto descendente e $\theta\geq0$ no riser. Pressões são absolutas. A fração de gás é $\alpha$, com

$$
j_l=(1-\alpha)u_l,\qquad j_g=\alpha u_g,\qquad j=j_l+j_g,
$$

$$
\rho_m=(1-\alpha)\rho_l+\alpha\rho_g,\qquad
\mu_m=(1-\alpha)\mu_l+\alpha\mu_g.
$$

O sistema dimensional resolvido é

$$
\frac{\partial U}{\partial t}+\frac{\partial F(U)}{\partial s}=S(U,s),
\quad
U=\begin{bmatrix}(1-\alpha)\rho_l\\\alpha\rho_g\\\rho_lj_l+\rho_gj_g\end{bmatrix},
$$

$$
F=\begin{bmatrix}
\rho_lj_l\\
\rho_gj_g\\
(1-\alpha)\rho_lu_l^2+\alpha\rho_gu_g^2+p
\end{bmatrix},\qquad
S=\begin{bmatrix}0\\0\\-\rho_mg\sin\theta-\dfrac{2f_F}{D}\rho_mj|j|\end{bmatrix}.
$$

As fases não trocam massa. O atrito se opõe à velocidade superficial total com sinal, inclusive quando $j<0$. Esta forma reúne as equações de massa e momento de [B, eqs. 1–6 e 16–21], com a compressibilidade líquida de [D, eqs. 3.10–3.12 e 3.31–3.33]. O líquido incompressível assumido por Burr é uma hipótese diferente do líquido com velocidade acústica finita adotado aqui.

## 2. Unidades, EOS e interpretação das vazões

### 2.1 Unidades e escalas

Diâmetro, comprimentos e rugosidade absoluta entram em metros; $g$ em m/s²; viscosidades dinâmicas em Pa·s; tensão superficial em N/m; ângulos em radianos. Os operadores constitutivos recebem densidades e velocidades adimensionais, restauram SI internamente e reconvertem apenas o resultado. Isso impede que a escolha arbitrária de `omega_u` ou `omega_rho` altere o fechamento físico.

Com escalas positivas $w_u,w_\rho,w_p$,

$$
s^*=s/L_r,\quad t^*=tw_u/L_r,\quad
u^*=u/w_u,\quad \rho^*=\rho/w_\rho,\quad p^*=p/w_p.
$$

As escalas de $U_1,U_2,U_3$ são, respectivamente, $w_\rho,w_\rho,w_\rho w_u$; as de $F_1,F_2,F_3$ são $w_\rho w_u,w_\rho w_u,w_\rho w_u^2$. O termo de pressão no fluxo de momento é

$$
\Pi_Pp^*,\qquad \Pi_P=\frac{w_p}{w_\rho w_u^2}.
$$

Os coeficientes `nCg` e `nCl` da EOS são $C_g\sqrt{w_\rho/w_p}$ e $C_l\sqrt{w_\rho/w_p}$, respectivamente. Não são, em geral, $C_g/w_u$ e $C_l/w_u$. O CFL usa os autovalores da Jacobiana do fluxo, que já incorporam $\Pi_P$.

### 2.2 Equações de estado

$$
\rho_l=\rho_{l0}+\frac{p-p_{l0}}{C_l^2},\qquad
\rho_g=\frac{p}{C_g^2}.
$$

`gas_eos='ISOTHERMAL'` define $C_g^2=R_gT$. `gas_eos='BAROTROPIC'` usa `Cg_fixed` explicitamente, permitindo reproduzir o valor $343\,\mathrm{m/s}$ da tabela 6.1 da dissertação. A densidade do gás ideal isotérmico é $p/(R_gT)$; inserir $\gamma$ nessa expressão altera a EOS. Os dois perfis devem ser identificados nos resultados.

A inversão das conservativas resolve a restrição de volumes

$$
\frac{U_1}{\rho_l(p)}+\frac{U_2}{\rho_g(p)}=1,
$$

com avaliação estável da raiz física da equação quadrática. Densidades ou pressão não positivas são falhas. Não há reparação silenciosa por recorte das massas conservadas.

### 2.3 Vazões de entrada

Os dados `mul` e `mug` representam **fluxos de massa por área**, em kg/(m²·s) antes da adimensionalização. Multiplicar por $A=\pi D^2/4$ fornece kg/s. Quando informados explicitamente, substituem a conversão de velocidades superficiais.

Por convenção declarada para o líquido barotrópico, $G_l=\rho_{l0}j_{l0}$; portanto, a velocidade superficial local é $j_l=G_l/\rho_l(p)$. Para gás com `gas_flow_basis='STANDARD'`,

$$
G_g=\frac{P_{\rm std}}{R_gT_{\rm std}}j_{g0}.
$$

As condições padrão de [BNB, eqs. 78–79] são $P_{\rm std}=1.013\,\mathrm{bar}$ absoluto e $T_{\rm std}=20\,^{\circ}\mathrm C$. `gas_flow_basis='SEPARATOR'` representa outra definição: $G_g=\rho_g(P_s)j_{g0}$. Escolher essa opção muda a massa de entrada; não é uma simples troca de unidade. O valor padrão de $G_g$ deve permanecer o mesmo ao comparar EOS diferentes.

## 3. Atrito de Fanning e limite de velocidade nula

O código usa o fator de **Fanning** em todas as tensões e fontes:

$$
f_F=16/Re\quad (Re<2000),\qquad
\tau_w=\tfrac12f_F\rho u|u|.
$$

O fator de Darcy é $4f_F$. Usar $64/Re$ com a fonte $2f_F\rho_mj|j|/D$, ou dividir essa fonte por quatro, produz outro modelo. A convenção acima é compatível com os balanços e tensões de [B, eqs. 3, 8 e 18]. A ocorrência de $64/Re$ na dissertação deve ser reconciliada com a convenção de tensões antes de uma comparação literal.

Entre $Re=2000$ e $2300$, há interpolação linear entre $16/2000$ e o valor turbulento em 2300. Para $Re>2300$, foi preservada a aproximação explícita turbulenta do código fornecido:

$$
f_F=\left[-4\log_{10}\left(\frac{\varepsilon/D}{3.7065}
-\frac{5.0452}{Re}\log_{10}\left(\frac{5.8506}{Re^{0.8981}}+
\frac{(\varepsilon/D)^{1.1098}}{2.8257}\right)\right)\right]^{-2}.
$$

Os limites de transição e essa aproximação são escolhas complementares do programa; Burr não os especifica. `Dffan` deriva a mesma interpolação e conserva sua dependência da rugosidade. A função de fator isolado exige $Re>0$. Em velocidade nula, o operador físico avalia o limite laminar diretamente:

$$
\tau_w=8\mu u/D_h,\qquad
S_{3,\rm atrito}=-32\mu_mj/D^2.
$$

Essas expressões são finitas e têm derivada definida em zero, sem formar $0\times\infty$.

## 4. Modelos de fechamento disponíveis

| ID público | Nome | Uso e fonte |
|---|---|---|
| `1` | `STRATIFIED` | Equilíbrio de paredes, interface e gravidade de Burr, eqs. 8–13. |
| `2` | `BENDIKSEN` | Variante publicada por Baliño, Burr e Nemoto (2010), eqs. 32–36. |
| `3` | `BHAGWAT` | Transcrição de Azevedo, Baliño e Burr (2015), p.121, eqs.116–128, seção circular e domínio cocorrente. |
| `4` | `MONTGOMERY` | Variante impressa na tese de Montgomery (2002), p.151, eqs. 5.43–5.48. |

O identificador é propagado à recuperação primitiva, à lei e às derivadas. Para chamadas antigas que o omitem, $\theta<0$ seleciona `1` e $\theta\geq0$ seleciona `2`. Essa hipótese fixa de escoamento por trecho não constitui implementação validada do critério dinâmico de transição de regime proposto em outras partes de Burr.

### 4.1 Estratificado: `STRATIFIED=1`

A geometria circular satisfaz

$$
\alpha=1-\gamma+\frac{\sin(2\pi\gamma)}{2\pi},\qquad
\gamma_i=\frac{\sin(\pi\gamma)}{\pi},\qquad
D_l=\frac{(1-\alpha)D}{\gamma},\quad D_g=\frac{\alpha D}{1-\gamma+\gamma_i}.
$$

A raiz geométrica é obtida por bisseção em $0\leq\gamma\leq1$. O residual dimensional de [B, eq.8], dividido por $\pi D$, é

$$
\Phi_{\rm SI}=
\frac{\tau_g(1-\gamma)}{\alpha}-\frac{\tau_l\gamma}{1-\alpha}
+\tau_i\gamma_i\left(\frac{1}{1-\alpha}+\frac{1}{\alpha}\right)
-\frac{(\rho_l-\rho_g)gA\sin\theta}{\pi D}=0.
$$

Aqui, $\tau_l$ e $\tau_g$ usam Fanning e as velocidades com sinal; a tensão interfacial é

$$
\tau_i=\tfrac12 f_i\rho_g(u_g-u_i)|u_g-u_i|,\qquad f_i=0.0142.
$$

O valor constante de $f_i$ é o de [D, eq.3.14, p.63] e está identificado por `FI_FANNING` no módulo. Com $Re_l=\rho_l|u_l|D_l/\mu_l$, usa-se $u_i=1.8u_l$ até 2000, interpolação linear até 2200 e $u_i=u_l$ a partir de 2200. Isso corrige o intervalo incompatível impresso na eq.9 de Burr, mantendo a transição de [D, eq.3.18].

A função pública retorna $\Phi_{\rm SI}/(w_\rho w_u^2)$. O termo gravitacional tem o mesmo fator de normalização das tensões; a conversão $\beta=-\theta$ ocorre uma única vez. São exigidos $0<\alpha<1$, propriedades positivas, seção circular coerente e $|\theta|<90^\circ$.

### 4.2 Bendiksen conforme Baliño: `BENDIKSEN=2`

Definem-se $Fr_j=|j_l+j_g|/\sqrt{gD}$ e $u_g=C_0j+U_d$. As expressões de [BNB, eqs.32–36, p.647] são

$$
\begin{array}{lll}
Fr_j<3.5:& C_0=1.05+0.15\sin\theta,&
U_d=\sqrt{gD}(0.35\sin\theta+0.54\cos\theta),\\
Fr_j\geq3.5:& C_0=1.2,&U_d=0.35\sqrt{gD}\sin\theta.
\end{array}
$$

A igualdade em 3.5 foi atribuída ao ramo superior como decisão numérica explícita. A derivada exatamente na mudança de ramo é sinalizada como indefinida. O residual adimensional retornado é

$$
\Phi=\alpha\left[-u_g^*+C_0\big((1-\alpha)u_l^*+\alpha u_g^*\big)+U_d/w_u\right].
$$

O sinal global desse residual não altera sua raiz; suas derivadas usam exatamente a mesma convenção. `CdUd_drift` e `CdUd_bendiksen` devolvem $U_d$ **dimensional**, em m/s. Para $g=0$, o código define explicitamente o limite $C_0=1.2, U_d=0$; essa extensão não constitui validação para microgravidade.

### 4.3 Variante Montgomery: `MONTGOMERY=4`

A variante de [M, p.151] usa $Fr_l=|j_l|/\sqrt{gD}$, $C_0=1.05+0.15\sin^2\theta$ abaixo de 3.5, $C_0=1.2$ acima e $U_d=\sqrt{gD}(0.35\sin\theta+0.54\cos\theta)$ nos dois ramos. A fonte apresenta $j_l$ positivo para escoamento ascendente; o módulo usa seu módulo apenas no número de Froude, retendo o sinal das velocidades no residual. Essa extensão algébrica a velocidades negativas não demonstra validade experimental em contracorrente.

As variantes `2` e `4` são distintas: trocar $\sin\theta$ por $\sin^2\theta$, a velocidade no Froude ou o ramo de $U_d$ muda o modelo. O intervalo angular computacional admitido para esses operadores é $[-30^\circ,90^\circ]$; admissibilidade computacional e validação física não são equivalentes. A tensão superficial é preservada na interface das rotinas, mas não participa dessas duas correlações nem do equilíbrio estratificado aqui implementado.

### 4.4 Bhagwat conforme Azevedo, Baliño e Burr: `BHAGWAT=3`

A dissertação escolhe Bhagwat (2014) nas pp.62 e 65. A implementação usa a transcrição completa de [ABB, seção 6.3], conferida visualmente, cujo PDF e imagem têm proveniência registrada em [BHAGWAT_SOURCE_STATUS.md](BHAGWAT_SOURCE_STATUS.md). Esse documento contém todas as equações, limites e diferenças em relação aos códigos fornecidos. O fechamento tem o mesmo residual de deriva da seção 4.2, mas seus coeficientes dependem de $\alpha$, densidades, velocidades, viscosidade e tensão superficial.

O fator $f_{tp}$ da eq.118 de [ABB] resolve Colebrook implicitamente e é um parâmetro interno da correlação; ele não substitui nem reutiliza automaticamente a regra de atrito de parede da seção 3. Coeficientes, residual e derivadas são calculados em SI. As funções `CdUd_bhagwat`, `CdUd_drift(...,closure_id=3)` e os nomes legados contendo `swananda` encaminham para a mesma versão identificada e devolvem $U_d$ dimensional.

O domínio qualificado admite $0<\alpha<1$, $\rho_l>\rho_g>0$, propriedades positivas, $-50^\circ\le\theta\le90^\circ$, $j_l,j_g\ge0$ e $j_l+j_g>0$. Uma fase estagnada é admitida se a outra tem fluxo positivo. Não se usam módulos nem recortes para tornar estados contracorrentes avaliáveis. Fora desse domínio há falha explícita; a abrangência angular da correlação não comprova validade para filme líquido descendente.

### 4.5 Seleção de raiz e contracorrente

As velocidades podem ter sinais opostos quando existe solução admissível da lei escolhida, como no estratificado; o modelo `3` tem a restrição específica acima. A seleção da raiz precisa ser explícita: `solve_alpha_from_superficials` retorna `STEADY_MULTIPLE_ROOTS` quando detecta mais de uma raiz e não recebeu indicação de ramo. O opcional `alpha_hint` seleciona a raiz mais próxima em coordenada logit, como continuação do estado indicado pelo chamador; um empate continua sendo falha. `nroots` solicita a busca diagnóstica completa. A busca finita por mudanças de sinal não prova ausência de raízes tangentes. O opcional `phase_hint` distingue resolução bifásica (`PHASE_AUTO=0`) de ausência explicitamente declarada de gás (`PHASE_LIQUID_ONLY=1`) ou líquido (`PHASE_GAS_ONLY=2`). Fluxo de uma fase igual a zero não significa que sua massa ou fração volumétrica sejam zero: a fase pode estar estagnada.

Esse tratamento evita descartar contracorrente por um recorte de velocidades, mas **não qualifica os modelos para filme líquido descendente durante SS1/SS2**. Essa limitação já é discutida em [D, pp.139–140]. Fechamento constitutivo, seleção de ramo e validade experimental são verificações diferentes.

## 5. Derivadas, recuperação primitiva e falhas

As derivadas de `Dlei_fechamento_or_ndim_simp` têm ordem pública

```fortran
dfda, dfdrhog, dfdrhol, dfdug, dfdul
```

As leis de deriva usam derivadas analíticas no ramo selecionado, incluindo diferenciação implícita de Colebrook e das frações mássica e volumétrica para o modelo `3`. O estratificado usa derivadas analíticas da geometria circular, dos diâmetros hidráulicos, das tensões de Fanning, da velocidade e tensão interfaciais e da gravidade, com regra da cadeia na reconversão. Isso evita misturar inclinações de ramos de Reynolds. A interpolação linear original de Fanning permanece inalterada: é contínua em valor e pode não possuir derivada bilateral nos limiares. Os testes próximos a 2000, 2200 e 2300 usam diferenças unilaterais de segunda ordem no mesmo ramo. Saltos de valor, como o Froude de Bendiksen ou $C_4$ descendente de Bhagwat, têm diagnóstico específico de derivada indefinida. `module_jacobian` deriva implicitamente EOS, momento e fechamento; as matrizes são verificadas contra perturbações das conservativas.

| Código de fechamento | Significado |
|---|---|
| `0` | Avaliação concluída. |
| `10` | Entrada inválida: escala, geometria, propriedade ou estado inadmissível. |
| `11` | Modelo não disponível ou identificador não suportado. |
| `12` | Estado fora do domínio admitido. |
| `13` | Derivada não definida no limiar de troca de ramo. |

O módulo estacionário acrescenta `20` para entrada inválida, `21` para ausência de raiz encontrada, `22` para raízes múltiplas sem seleção inequívoca, `23` para singularidade e `24` para falha de integração. Primitivas, contornos e álgebra linear também têm códigos próprios; o contexto da chamada faz parte do diagnóstico. O chamador deve sempre verificar `status` ou `INFO` antes de consumir saídas. Valores indefinidos são sinalizados com quiet NaN, sem convertê-los em zero e declarar sucesso.

Os limites monofásicos são implementados na EOS, na recuperação primitiva e no inicializador estacionário. Eles não tornam automaticamente diferenciável o fechamento bifásico em $\alpha=0$ ou $1$. Operadores que exigem essa derivada precisam rejeitar o estado ou adotar uma formulação de limite explicitamente testada.

## 6. Catenária, inicialização e junção

Para a geometria de [D, eqs.3.26–3.30, pp.66–67],

$$
z=a[\cosh(x/a)-1],\quad s_r=a\sinh(x/a),\quad
\theta(s_r)=\arctan(s_r/a),\quad \theta'(s_r)=\frac{a}{a^2+s_r^2}.
$$

A malha possui uma face exatamente na junção e outra no topo; centros e faces têm papéis distintos. A inicialização fixa $p=P_s$ na **face** do topo e integra o riser até a junção; depois integra o oleoduto com a mesma pressão de junção e seu próprio fechamento.

No estacionário, $G_l=\rho_lj_l$ e $G_g=\rho_gj_g$ são constantes. Definindo o fluxo de momento

$$
H(p,s)=G_lu_l(p,s)+G_gu_g(p,s)+p,
$$

a conservação fornece

$$
\frac{dp}{ds}=\frac{S_3-\left.\partial H/\partial s\right|_p}
{\left.\partial H/\partial p\right|_s}.
$$

As derivadas recompõem o fechamento com o ângulo local. O termo $H_s$ inclui a mudança de inclinação na catenária; omiti-lo não satisfaz o balanço total de momento. As duas fases conservam massa mesmo com líquido compressível, em vez de impor $j_l$ constante enquanto $\rho_l$ varia. A integração falha explicitamente diante de pressão inválida, fechamento insolúvel ou denominador singular.

As condições de compatibilidade de [B, eqs.27–29] igualam pressão e velocidades superficiais, permitindo salto de $\alpha$. A opção `MATCHED` impõe dois traços de fechamento na face de junção, pressão e fluxos de fase comuns, mas permite dois fluxos de momento. A diferença representa uma força localizada implícita na hipótese de acoplamento. Portanto, o momento global deve incluir esse salto de junção; a soma das massas continua telescópica. Não se deve interpretar o acoplamento como uma face homogênea com um único fluxo de momento.

## 7. Contornos, conservação e integração temporal

`boundary_flux` aplica condições na face, sem modificar células físicas. Na entrada, recupera a pressão interior e impõe $F_1=G_l$, $F_2=G_g$ independentemente. Reconstrói o fechamento a essa pressão e calcula $F_3=G_lu_l+G_gu_g+p$. A raiz usa a fração interior como indicação explícita de continuação de ramo. Um limite monofásico é declarado somente quando a massa interior da fase e seu fluxo imposto são ambos zero. Não existe cache global para fixar massas de uma execução anterior.

No topo, impõe $P_s$, conserva a fração de gás e o momento extrapolados do interior, recompõe densidades e velocidades pelo fechamento e calcula o fluxo físico desse estado. `apply_ghost_cells` apenas copia as células físicas adjacentes para os fantasmas. Não aplica pisos positivos, não zera momento e não altera a primeira ou a última célula física.

O transiente usa volumes finitos conservativos e Euler implícito. Todas as células físicas pertencem ao sistema de Newton. O residual de uma fase soma, com os volumes, exatamente a variação de seu inventário e os fluxos externos; fontes de fase são zero. As Jacobianas de face podem incluir a variação completa da dissipação numérica, verificada por perturbação do fluxo completo. A opção de fluxo e de junção deve acompanhar os resultados, pois altera a aproximação numérica.

O CFL usa o raio espectral calculado por LAPACK `DGEEV`, sem substituir autovalores pelo maior elemento da matriz nem reduzir arbitrariamente a velocidade acústica. A estrutura característica de Roe verifica espectro real, invertibilidade e resíduo dos autovetores. `DGESV` resolve sistemas densos pequenos; `DGBSV` resolve o sistema global de blocos 3×3 com pivotamento, $KL=KU=5$, `LDAB=16`. Não se regularizam pivôs somando um epsilon ao sistema físico.

## 8. Reprodução dos testes

O script exige os seis programas de teste; um arquivo ausente produz falha, nunca um teste omitido silenciosamente:

```sh
FC=gfortran LIBS='-llapack -lblas' python3 scripts/run_unit_tests.py
```

Em instalações Linux que oferecem apenas as bibliotecas de execução, pode-se usar:

```sh
FC=gfortran LIBS='-l:liblapack.so.3 -l:libblas.so.3' python3 scripts/run_unit_tests.py
```

O compilador pode ser um caminho absoluto ou um comando em `FC`. A qualificação desta revisão usa GNU Fortran 14, LAPACK e pilha não executável. O código de produção não exige trampolins de procedimentos; `CALLBACK_FLAGS` é vazio por padrão. `LDFLAGS` pode ser ajustado para outro compilador previamente qualificado. A compilação usa `-std=f2018`, `-fcheck=all`, rastreamento de falhas e exceções para operações inválidas, divisão por zero e overflow; não inclui `main.f90`.

| Teste obrigatório | Evidência principal |
|---|---|
| `test_geometry` | Parâmetro e comprimento da catenária, coordenadas, inclinações e topologia da malha. |
| `test_closure_physics` | Fanning, limites com sinal, gravidade, escalas, deriva, Bhagwat contra referências de 80 dígitos e derivadas do mesmo operador. |
| `test_primitives_numerics` | Inversão EOS/fechamento, Jacobiana, espectro, dissipação e sistemas densos. |
| `test_steady` | Hidrostática analítica, pressão de face, geometria, massas constantes, dois traços e termo espacial da catenária. |
| `test_boundaries_time` | Fluxos de fase independentes, ausência de cache/alteração física, pressão de saída, CFL espectral e banda contra sistema denso. |
| `test_fv_operators` | Fontes e fluxos completos por diferenças finitas, junção compatível e soma conservativa nas células físicas. |

O diretório padrão `tmp/unit_tests` recebe executáveis, logs individuais e `summary.json`, contendo compilador, argumentos, duração, contagem de verificações, revisão Git e SHA-256 das fontes. Fontes alteradas durante a execução invalidam a qualificação dessa rodada. `--build-dir` e `--timeout` permitem alterar o diretório e o limite por teste.

Esses testes não aprovam, por si, períodos ou amplitudes de slugging. A aceitação física exige condição de vazão correta, janelas transientes suficientes, inventários e produção de cada fase, refinamento espacial e temporal independente, verificação de sensibilidade às tolerâncias e comparação com a mesma geometria/EOS/lei de fechamento. A tese de Montgomery usa outra instalação; seus dados não podem ser transferidos silenciosamente à catenária de [D]/[BNB].

## 9. Referências e limites da documentação disponível

- **[B]** Karl Peter Burr. *Modelagem Unidimensional de Sistema Oleoduto-Riser*. Anexo `03-sistema_oleoduto_riser.pdf`, 15 páginas. Conservação: eqs.1–6 e 16–21; estratificado: eqs.8–13, pp.2–3; fechamento de deriva: eq.23; contornos e junção: eqs.24–29, p.5; adimensionalização: eqs.30–72; estacionário: eqs.73–108, pp.11–15. Foram identificados erros de impressão nas desigualdades da velocidade interfacial e em partes da derivação estacionária; o código parte do balanço conservativo completo.
- **[D]** Leonardo Alegretti Belaroza. *Solução numérica para escoamentos multifásicos unidimensionais em sistemas oleoduto-riser*. Dissertação, Universidade Federal do ABC, agosto de 2025, versão revisada fornecida. EOS: eqs.3.10–3.12, p.62; estratificado e fator interfacial: eqs.3.14–3.21, pp.63–64; escolha Bhagwat e contornos: p.65; catenária: eqs.3.26–3.30, pp.66–67; conservativas: eqs.3.31–3.33, pp.68–69; parâmetros e casos: tabelas 6.1, 6.2 e 6.7, pp.108,109,138. Páginas citadas são impressas; no PDF fornecido há deslocamento de duas páginas.
- **[BNB]** J. L. Baliño, K. P. Burr e R. H. Nemoto (2010). *Modeling and simulation of severe slugging in air–water pipeline–riser systems*. International Journal of Multiphase Flow **36**, 643–660. [DOI:10.1016/j.ijmultiphaseflow.2010.04.003](https://doi.org/10.1016/j.ijmultiphaseflow.2010.04.003). Fechamento de Bendiksen: eqs.32–36, p.647; vazão padrão: eqs.78–79, p.649; casos experimentais: tabelas 3 e 5, pp.653–654. O conteúdo foi confrontado por transcrições indexadas do artigo; as tabelas numéricas com divergências de transcrição exigem conferência visual antes de aprovação quantitativa. O modelo NPW do artigo e o modelo conservativo com inércia desta implementação não são idênticos.
- **[M]** J. A. Montgomery (2002). *Severe Slugging and Unstable Flows in an S-Shaped Riser*. PhD thesis, Cranfield University. Variante de Bendiksen: seção 5.3.3, p.151 impressa / p.169 do PDF, eqs.5.43–5.48, verificadas visualmente no anexo.
- **[BG]** Swanand M. Bhagwat e Afshin J. Ghajar (2014). *A flow pattern independent drift flux model based void fraction correlation for a wide range of gas–liquid two phase flow*. International Journal of Multiphase Flow **59**, 186–205. [DOI:10.1016/j.ijmultiphaseflow.2013.11.001](https://doi.org/10.1016/j.ijmultiphaseflow.2013.11.001). Identificação e resumo confirmados; texto completo original não auditado. A implementação usa a transcrição precisamente identificada de [ABB].
- **[ABB]** G. R. Azevedo, J. L. Baliño e K. P. Burr (2015). *Linear stability analysis for severe slugging: sensitivity to void fraction correlations*. IMAACA 2015, pp.115–124, ISBN 978-88-97999-63-8. [Anais oficiais](https://www.msc-les.org/proceedings/imaaca/2015/IMAACA2015.pdf). Bhagwat: seção 6.3, eqs.116–128, p.121 impressa / p.135 do PDF, conferidas visualmente. É pesquisa primária de aplicação/comparação; não se confunde com o original da correlação de 2014.

Persistem limitações documentais que não devem ser preenchidas por ajuste de coeficientes: o texto original de Bhagwat de 2014 não foi confrontado integralmente, faltam séries experimentais brutas de produção no topo e há divergência no período SS3 entre a dissertação e transcrições do artigo de 2010. Os relatórios devem registrar quais perfis e versões de fechamento foram testados.

A busca complementar da fonte completa de Bhagwat está registrada em [BHAGWAT_SOURCE_STATUS.md](BHAGWAT_SOURCE_STATUS.md).
