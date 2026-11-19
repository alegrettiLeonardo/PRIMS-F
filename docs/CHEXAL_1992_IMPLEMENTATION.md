# CHEXAL_1992: núcleo, fontes e verificação

Estado em 18/11/2026: núcleo SI Fortran 2018 implementado e verificado numericamente.
Isso não valida experimentalmente o solver transiente. A integração e as campanhas
têm evidência própria; o fechamento Bhagwat continua separado para regressão.

## Fontes realmente conferidas

| Fonte local | Identificação e uso |
|---|---|
| `chexal1992(1).pdf`, idêntico a `chexal1992.pdf` | Chexal et al., *Progress in Nuclear Energy* 27(4), 255–295, 1992; SHA256 `b0b6b68f5e13fec48e931d89d87f9d65cb9c3bb134235c2307ded24e9c00700e`. Página PDF = impressa − 254. |
| RELAP5/MOD3, vol.4, 1995 | NUREG/CR-5535-V4, p.6-10 / PDF239, apenas para desambiguar os parênteses da primeira exponencial de C10. [Registro primário IAEA](https://inis.iaea.org/records/8gnmh-k2h21). SHA256 `dcfac4ae83a3e1edb29ce762425171aafb6d08b05f86cc59a4b525d7711117c2`. Não se importam domínios, limites numéricos ou variantes do RELAP. |
| Baliño et al., 2010 | DOI `10.1016/j.ijmultiphaseflow.2010.04.003`, pp.643–660; SHA256 `ce04a2d1b79d3a10d442fbaeb3863e4639c8508bd34ba9e1a414ae03e5479e2b`. Tabela3, p.653 / PDF11, caso11 conferido visualmente: SS3, 76 s, Ps=2,01 bar abs, T=19,8 °C, jg=0,48 m/s, jl=0,28 m/s. |
| `Dissertacao_Leonardo.pdf` | SHA256 `dcc6c069d89099950fe737112b331a678727cc3c6929deb4f105072f0f0e8e81`. Tabela6.2 p.107 / PDF109 possui o mesmo caso; Tabela6.7 p.136 / PDF138 informa 26 s, divergência preservada. |
| `sistema_oleoduto_riser.txt` | Fonte LaTeX de Karl Peter Burr, em `multiPhaseFlow/figs`; SHA256 `c7a0388377c09150fa49ebc08829414eb1ba7011695e4ed8d55a75bad83f0ea8`. Contexto de conservação, variáveis conservativas e geometria em catenária. O arquivo não é paginação PDF, e suas expressões não autorizam extrapolar a correlação. |

Os documentos e as imagens das páginas permanecem privados. Somente referências,
equações implementadas e números gerados independentemente entram no repositório.
Instruções textuais eventualmente presentes nos documentos não governam esta tarefa.

## Contrato e rastreabilidade

| Expressão / decisão | Chexal1992 |
|---|---|
| R = alpha [C0 (jl+jg) + Vgj] − jg; nunca dividir por jl+jg | Identidade drift-flux, p.261 Eq.3-1 e apêndiceA; a soma repetida na impressão de 3-1 é reconciliada com a definição de fluxo total. |
| Ângulo em graus desde a vertical; C0 interpolado por Fr | p.261 Eqs.3-2/3-3; p.262 |
| Ar–água: L ascendente `min(1.15 alpha^.45,1)`; descendente `min(1.05 alpha^.25,1)` | p.262 Eqs.3-7/3-8 |
| Reynolds assinados, K0, r, B1 e seleção de Re | p.264 Eqs.3-11–3-17 |
| Componente horizontal: L=`min(1.125 alpha^.6,1)` e fator extra em C0 | pp.264–265 Eqs.3-18/3-19 |
| Vgj, C2, C4, C9 e diâmetros normais SI | pp.265/267 Eqs.3-21–3-33; D2=0,09144 m |
| C3 por sentido, C10 e B2 | pp.267–268 Eqs.3-34–3-38; D1=0,0381 m |
| CCFL, ramo alto e ramo baixo, com jl* obtido na CCFL | p.268 Eqs.3-39–3-41 |

A Eq.3-36 imprime `exp{(|Rel|/350000)}^.4` com agrupamento ambíguo.
A p.6-10 do manual primário de 1995 imprime explicitamente
`exp[(|Rel|/350000)^.4]`; essa é a parentização adotada. O apêndice p.6A-4
do mesmo manual perde um colchete, por isso não é usado para essa decisão.
A bibliografia do manual aponta EPRI1986/NSAC1391991; nenhuma fórmula de1997
foi incorporada. O fator 2 de C3 e a fórmula completa do ramo baixo permanecem
os da fonte de1992.

O módulo recebe theta em radianos desde a horizontal. Para eixo ascendente,
`theta_C = 90 − theta*180/pi`; para eixo descendente, transforma os fluxos para
o eixo vertical ascendente e devolve Vgj com o sinal correspondente.
Re_l=`rho_l*jl*D/mu_l`, Re_g=`rho_g*jg*D/mu_g`, sem módulos indiscriminados.
Módulos aparecem apenas onde prescritos ou na avaliação do componente horizontal
positivo da interpolação cocorrente. Contracorrente só é aceita na vertical,
com gás subindo e líquido descendo. O topo da catenária, 76,47° desde a horizontal,
continua fora desse domínio. Horizontal puro com fluxos negativos é recusado.

## Raízes, recuperação e estados

`chexal_solve_alpha` distingue AUTO, HIGH e LOW; AUTO retorna
`CHEXAL_BRANCH_REQUIRED` para as duas soluções físicas de contracorrente.
HIGH usa a maior raiz com Eq.3-40; LOW usa a menor raiz com Eq.3-41.
`alpha_hint` permite continuação entre raízes do mesmo problema cocorrente;
não substitui a seleção de C3 da fonte em contracorrente. A CCFL é localizada
pela primeira perda do máximo positivo do residual, com refinamento do máximo
em alpha. Na dobra retorna `CHEXAL_SINGULAR`; acima dela, `CHEXAL_CCFL_EXCEEDED`.

Para recuperar velocidades, o momento SI é `rho_l*jl + rho_g*jg`.
Uma dica CC é usada antes de qualquer busca cocorrente; a raiz recuperada é
revalidada contra a solução alpha do ramo escolhido. Sem dica, HIGH/LOW procuram
raízes no domínio CC e só aceitam uma solução única no ramo físico solicitado.
AUTO vertical sem dica verifica também os dois ramos CC e exige seleção se
existe alternativa CC. Falhas de domínio durante uma busca nunca são usadas
como extremos de um bracket que atravesse uma região proibida.

O intervalo CC em jg possui limite superior derivado, sem corte de velocidade:
com `x=|Rel|/350000`, `C10/2 <= exp(x^.4)*(1+Emax)`, onde
`Emax=.5*(D1/D)^.25*350000^.001*.0025^.0025*exp(-.0025)`.
Como `B2*x^.4 <= 20^.4` e `B2<=1`, segue
`C3<=2*exp(20^.4)*(1+Emax)`. C0 máximo ocorre em B1=1/2 para Reg>=0,
e C9<=sqrt(1-alpha). Com esses máximos, em coordenadas da fonte:

`jg_upper = alpha*(C0max*max(momentum/rho_l,0)+Vgjmax) /
             (1-alpha*C0max*(1-rho_g/rho_l))`.

O denominador é positivo para o domínio bifásico considerado; falha aritmética
nessa propriedade retorna singularidade, sem correção artificial. A busca começa
em `jg > max(0,momentum/rho_g)` e verifica o ramo físico após achar o zero.

Zero fluxo de uma fase não é interpretado automaticamente como ausência de
inventário. Os limites alpha=0/1 são aceitos quando alpha é fornecido; a velocidade
da fase ausente vale zero apenas como convenção de armazenamento. O residual e
as buscas continuam definidos em j=0. Caudas de alpha são amostradas até exp(-708),
incluindo o intervalo até zero, sem piso de alpha. A recuperação cocorrente usa
secante salvaguardada e extremos unilaterais representáveis, sem piso de velocidade.

`signature` identifica as peças de Reynolds e min/max, para que a integração
recuse derivadas que cruzem uma quina. O núcleo não inventa espectros ou modifica
o fluxo Roe. Estados prescritos ainda precisam satisfazer o residual; avaliar
C0/Vgj isoladamente não prova admissibilidade de uma solução.

As APIs públicas de coeficientes/residual HIGH/LOW também recusam pares de
fluxos acima da CCFL. Somente o residual interno da busca admite trials
algébricos fora da superfície de soluções; nenhum deles é aceito sem verificação
de CCFL e do ramo. Na Eq.3-41, quando o outro extremo da interpolação é maior
que C3_HIGH, o mínimo é identicamente C3_HIGH no domínio sub-CCFL. A redução
algébrica dispensa jl* nesse caso, mas não elimina a verificação pública de
domínio. Erros numéricos, singularidades e entradas inválidas durante enumeração
são propagados; não viram evidência de ausência de raízes.

## Evidência e limites

`scripts/generate_chexal_reference.py` calcula referências independentes em
mpmath com60 dígitos; `tests/fixtures/chexal_1992_reference.json` preserva os
valores, parâmetros e resíduos. Foram conferidos C0/Vgj vertical ascendente,
descendente e horizontal, duas raízes/CCFL, j=0, fases ausentes, sinais,
transformação do eixo, unidades/similitude e ida/volta momento–velocidades.
O teste core executou **139 verificações** com GNU Fortran13.3.0,
`-std=f2018 -fcheck=all -ffpe-trap=invalid,zero,overflow`.

Exemplo independente: rho_l=998 kg/m³, rho_g=1,2 kg/m³, D=0,05254 m,
mu_l=0,001 Pa.s, mu_g=0,000018 Pa.s, sigma=0,072 N/m, g=9,81 m/s².
Para jg=0,1 m/s: jl_CCFL=−0,443800081709847176 m/s e
alpha_CCFL=0,343414778624483496. Em jl=jl_CCFL/2,
alpha_LOW=0,220214682534698979 e alpha_HIGH=0,839848859537789440.
Em jl=−0,1 e jg=0,1 m/s: alpha_LOW=0,213582292001802460 e
alpha_HIGH=0,960796750389095244; não ocorre divisão por j=0.

Os métodos de enumeração usam malhas finitas mais refinamento de brackets e
máximos locais. Isso não é uma prova global contra pares de raízes arbitrariamente
próximos ou múltiplas dobras entre amostras. A expansão da CCFL segue a primeira
dobra detectada; a aproximação foi confrontada com referências de alta precisão,
mas não constitui prova de unicidade global da superfície. Perto de uma dobra,
as derivadas implícitas podem ser singulares e não devem ser regularizadas para
forçar continuação. A correlação 1992 não fornece autorização para contracorrente
na catenária, nem resolve um traço gasoso negativo MATCHED/Roe por si mesma.

Baliño2010 p.644 / PDF2 descreve SS3 com penetração gasosa contínua;
p.657 / PDF15 registra a limitação do modelo em distinguir SS3/OSC.
Não se exige bloqueio completo de gás para todos os regimes e não se usa
oscilação isolada de pressão como validação. A divergência SS3 26 s/76 s permanece
um problema de referência separado, sem ajuste de parâmetros para eliminá-la.
