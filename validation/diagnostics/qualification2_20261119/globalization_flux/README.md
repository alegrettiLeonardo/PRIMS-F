> Registro histórico da versão3940bb7 (L2 tardio). Para o gatilho antecipado,
> minimax e limites dos smokes v1/v2, consulte [o adendo](ACTIVE_SET_MINIMAX.md).

# Direção ativa nos mesmos snapshots SS1/SS3

O operador `69813ee` foi reproduzido sem alterar suas fontes: os dois
vetores de resíduo coincidem **bit a bit** com os snapshots terminais.
A opção `newton_globalization='ACTIVE_SET'`, congelada em `3940bb7`,
modifica apenas a direção após esgotar a busca de Newton original.
O padrão permanece `LINE_SEARCH`. Não há projeção de fluxo, clipping de
massa, alteração de fechamento, tolerância ou coeficiente de dissipação.

## Problema quadrático e solução

Sejam `s` as escalas originais do resíduo, fixadas por `Uold`:
`[max(ml,1e-3), max(mg,1e-7), max(abs(M),1)]` em cada célula.
Usam-se `S=diag(tol_res*s)`, `W=S^-1`, `A=W J S`, `b=W R`.
Para uma margem `g(U)>=0`, `c=S grad(g)/gscale` e
`gbar=g/gscale`, resolve-se

\[
 \min_z \tfrac12\|Az+b\|_2^2,\qquad gbar+c^Tz\ge0.
\]

O Newton já calculado fornece `z0=S^-1 d_Newton`. Se sua margem
linearizada é negativa, as condições KKT dão

\[
 q=A^{-T}c,\quad \mu={-gbar-c^Tz_0\over q^Tq},\quad
 z=z_0+A^{-1}q\mu,\qquad d=Sz.
\]

Isso é verificado substituindo na restrição e em
`A^T(Az+b)-mu*c=0`; `mu>=0`. O código resolve duas vezes o sistema
em banda adicionalmente ao Newton existente. Não forma `A^T A`.
Normaliza `q` antes do produto interno para evitar overflow de `q^Tq`.
Singularidade e gradiente não resolvido retornam status, sem regularização
oculta. `gscale` apenas normaliza a mesma desigualdade; `S` e `W` ficam
fixos durante todo o passo temporal.

O tratamento por subproblema quadrático com restrições linearizadas se
relaciona a Byrd, Hribar e Nocedal, *An Interior Point Algorithm for
Large-Scale Nonlinear Programming*, SIAM J. Optimization 9(4),
877–900 (1999), [DOI](https://doi.org/10.1137/S1052623497325107),
[manuscrito primário, §2](https://users.iems.northwestern.edu/~nocedal/PDFfiles/nitro.pdf).
A fórmula acima é derivada aqui para **uma** restrição ativa. Não se
implementou o algoritmo completo de barreira/trust-region daquele artigo,
nem se invoca sua garantia global para o PRIMS-F.

## Margens e critérios de aceite

A opção atual é específica de `MATCHED`, pipe estratificado (ID1),
riser Bhagwat (ID3), dissipação `LAGGED`.

- SS1: `g=FM(2)`, o fluxo mássico gasoso bruto comum da junção. A rotina
  de margem admite status físico 0 ou 1012 somente, pois 1012 ocorre
  depois do cálculo desse fluxo e recusa seu traço Bhagwat. Qualquer outro
  status deixa a margem indisponível. `g>=0` restringe o traço cocorrente;
  **não** equivale a positividade das massas médias.
- SS3: no riser ascendente, `g=ug_max-C0(alpha,ul=0,ug_max)*alpha*ug_max-Ud`
  em SI, `ug_max=M/mg`. É a margem no endpoint cocorrente da mesma lei.
  Trata a aproximação à fronteira da recuperação; não estende a lei a
  contrafluxo. Os coeficientes vêm da implementação já existente. O
  wrapper recebe velocidades normalizadas e devolve `Ud` dimensional.

O gradiente usa diferenças centrais de quarta ordem em seis variáveis da
junção ou três da célula. Os stencils `h` e `h/2` devem concordar dentro
de `1e-5*max(abs(dc),abs(df),1e-8)`, até oito reduções; a discordância
recusa a alternativa. `h=1e-5*max(abs(Uk),floor_k)`, com pisos `1e-5`
para massa gasosa e `1e-2` para as demais componentes.

Depois das 18 propostas do raio original, a alternativa tem até 18
amortecimentos próprios. A fração de segurança de massa original (0,90)
permanece. Cada proposta monta o **resíduo conservativo original**,
recusa estados físicos inválidos e exige também a margem não linear
selecionada não negativa. Usa a mesma condição de redução e só conclui o
passo quando a norma original é `<=tol_res`. Uma margem disponível não
substitui a recuperação de primitivas nem o teste do residual.

Seleciona-se primeiro a junção, se violada, senão a primeira célula do
riser cuja proposta cruza o endpoint. É uma alternativa com uma
restrição ativa por direção; não há garantia de resolver múltiplas
fronteiras simultâneas. O ledger `active_set.csv` registra tipo, índice,
status, margens, amortecimento e resultado, separado dos status físicos.
Sua coluna `accepted` indica aceitação da **proposta Newton por descida**;
somente `steps.csv` e `conservation.csv` registram um passo temporal
concluído com a tolerância original.

## Replays executados

Snapshots originais:

- SS1 SHA256 `1e9629e5293c8ee55ae324b6094b0a73d46aee47b6673dda2e5ee456917f27dc`,
  `t=78.57178847673931 s`, `dt=1.3286025e-7 s`.
- SS3 SHA256 `e2cafe7245ac084fe8b8db764d4d19bbc6d4dc6e0f891ec4a58ad68f054fac3c`,
  `t=97.1318523387001 s`, `dt=1.10716875e-7 s`.

| Mesmo estado e dt | Raio Newton original | Direção ativa Fortran |
|---|---|---|
| SS1 | proposta integral 3012, FACE53; norma inicial `2.976406330921032e-7` | 1 direção, lambda1; norma `1.6709772876859093e-12`; `jg*=+5.79677615946e-13 m/s` |
| SS3 | proposta integral 1017, CELL69; norma inicial `5.069257660924518e-8` | 2 direções, lambda0.5 e1; norma `4.0275959284536995e-9`; margem `+1.1102230246251565e-16 m/s` |

`tol_res=1e-8`. O resultado do replay é `ORIGINAL_RESIDUAL_CONVERGED`,
não apenas um critério de estacionariedade de otimizador. Não foi usado
SciPy neste caminho e não há `OptimizeResult` a atribuir a estes testes.
As massas médias de ambas as fases permanecem positivas.
O certificado anterior de ausência de raiz nas propostas SS3 originais
continua válido: a direção nova produz outros estados conservativos.
Não era uma prova de inexistência para todas as direções possíveis.

No passo SS1, os balanços de fase são `[-8.65959777666e-13,
+1.02922562338e-16] kg`; em SS3, `[+5.86371767560e-13,
-4.30475956770e-15] kg`. São os resíduos integrados reais, não zero
forçado. A identidade `delta inventario - dt*(entrada-saida) = sum(R*dx)`
concorda até `4.30e-23 kg`. Fluxos comuns cancelam entre as duas seções;
o momento retém as duas traças prescritas pelo modelo MATCHED.

Nos dois snapshots, multiplicar o passo do gradiente offline por
0,25, 0,5, 2 ou 4 preservou os mesmos valores finais acima. O teste
Fortran verifica também KKT com blocos não simétricos, direção já
admissível, sinais, escala inválida, gradiente não finito, dual finito
da ordem `1e200` e sistema singular: 11 verificações. A suíte completa
passou 14/14 programas.

## Comparação de fluxo

Harten, Lax e van Leer (1983), p.50, Eqs.(3.13)–(3.15),
[DOI](https://doi.org/10.1137/1025002),
[artigo primário na SIAM](https://epubs.siam.org/doi/10.1137/1025002),
define o estado intermediário e fluxo de duas ondas a partir da forma
integral. Aplicou-se essa álgebra ao **mesmo par mapeado** usado pelo
MATCHED. As velocidades usadas foram o envelope dos autovalores e
velocidades de fase dos endpoints, sem multiplicador espectral.
Não se demonstrou que esse envelope limita todas as ondas não lineares.

No SS1 aceito, HLL dá massas intermediárias positivas
`[433.2570016,1.35404615] kg/m3`. A atualização explícita das equações
de massa, usando o dt congelado e os demais fluxos originais, também
mantém massas positivas (mínimo gasoso `0.04470513036 kg/m3`). Porém
o fluxo gasoso comum é `-0.00700837780 kg/s`, `jg*=-1.339686934 m/s`:
o traço riser Bhagwat continua inadmissível. A equação de momento desse
traço não foi inventada nem considerada definida. No estado de
linearização ocorre a mesma distinção. O Rusanov já implementado, sem
fator extra, recusa o preparo com status1012 no mesmo snapshot.
Trocando somente a matriz de dissipação da junção por `rho(|A_Roe|) I`,
a montagem localiza a recusa em FACE53/status3012, com
`Gg*=-0.00701223655 kg/s`, enquanto as massas médias seguem positivas.

Além disso, a massa gasosa do estado direito mapeado é cerca de 30,585
vezes a massa da célula direita original. Um resultado de positividade
para o problema homogêneo mapeado não transfere automaticamente uma
região invariante ao acoplamento físico original. Não se promoveu HLL
ou Rusanov a solução dessa falha.

## Limites e reprodução

Convergir estes passos minúsculos não demonstra avanço temporal
sustentado, ciclos acomodados ou validade experimental. A margem SS3
final está na precisão de máquina; dados longos e sensibilidade de
tolerância são necessários para distinguir avanço e estagnação junto à
fronteira. O método não fornece a física de uma frente livre nem autoriza
contracorrente inclinada. Os smokes SS1/100s e SS3/110s são registrados
separadamente, com o mesmo binário, entradas e hashes.

`build_manifest.json`, `fortran_replay_manifest.json`,
`flux_comparison.json` e `gradient_sensitivity.json` identificam os
artefatos. Os scripts e fontes do probe estão nesta pasta; a biblioteca
e os executáveis ficam somente no diretório externo de resultados.
Nenhum PDF bibliográfico é incluído.
