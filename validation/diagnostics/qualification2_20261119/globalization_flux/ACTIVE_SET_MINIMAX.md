# Gatilho antecipado e direção minimax: qualificação numérica

Este adendo sucede o README histórico do operador `3940bb7`. As fontes
avaliadas nos replays v3 têm SHA256 `67267f05590043d07c129e2357bd93eba99aaf355c5f7e16c257894f0fa2168a`
(`module_active_set.f90`) e `46b1fecef59c980a9cbaf748b15600f6a808e159039a8cbbc84a7a2cf4ef796f`
(`module_simulate_pipe_riser.f90`). O commit final é `a179e51d9d63e8aeed91a94d829b5b483d107864`, com
build fresco cujo executável tem SHA256
`02ce51462ec3c954cec9364b17444b7279be72ec150f87e307f1b5fc8dd0a4c3`. Os replays aqui não são trajetórias v3 completas.

## Alteração restrita à globalização

`LINE_SEARCH` continua padrão. Em `ACTIVE_SET`, a primeira proposta Newton
recusada com status **3012 ou 1017** antecipa a tentativa ativa, desde que
a margem correspondente esteja disponível e sua linearização viole a
restrição. Um `NO_ROOT` genérico ou margem indisponível não autoriza a
alternativa. Se ela não produz descida admissível, restaura-se exatamente
a direção Newton e seu amortecimento para continuar a busca original.

O critério anterior ativava L2 somente depois de esgotar o raio Newton.
No terminal SS1 v1, pequenos passos de descida desse raio consumiam o
limite de 25 iterações, embora outra direção já resolvesse o mesmo passo.
O novo gatilho resolve esse replay em uma iteração. Não foram alterados
o resíduo, o Jacobiano, os fluxos, as leis físicas nem as tolerâncias.

## Derivação da segunda direção

Conservam-se as escalas fixas `S`, `W`, `A=WJS`, `b=WR`, `z0=S^-1 dNewton`
e a única desigualdade linearizada `gbar+c^T z>=0` descritas no README.
Definindo `r=Az+b`, `q=A^-T c` e `delta=-gbar-c^T z0>0`, a restrição
equivale a `q^T r>=delta`. A desigualdade de Hölder implica

\[
 \|r\|_\infty\ge\delta/\|q\|_1.
\]

O limite é atingido por `r_i=delta*sign(q_i)/||q||1` quando `q_i!=0`,
com `r_i=0` quando `q_i=0`; então `d=S*(z0+A^-1*r)`. Se `delta<=0`,
Newton já é admissível para esse subproblema. Se `q=0` com `delta>0`,
o subproblema é inviável. Não se substitui uma norma nula por piso.
A derivação é explícita para uma restrição linear; não é uma prova de
viabilidade do problema não linear nem um método de múltiplas restrições.
A referência primária para subproblemas com restrições linearizadas
permanece Byrd, Hribar e Nocedal (1999), §2,
[DOI](https://doi.org/10.1137/S1052623497325107),
[manuscrito](https://users.iems.northwestern.edu/~nocedal/PDFfiles/nitro.pdf).
Não se implementou o algoritmo completo de barreira daquele artigo.

L2 permanece a primeira direção. L∞ é tentada somente se a direção L2
válida esgotar sua busca sem descida admissível. Cada alternativa tem
até 18 amortecimentos; com as 18 propostas originais, o máximo é 54 por
iteração Newton. Cada direção ativa requer até duas soluções em banda;
o contador registra as soluções realmente tentadas. Os 25 Newton e a
fração de segurança de massa 0,90 permanecem inalterados. Esse custo
adicional é explícito e não implica melhora geral do tempo de execução.

`q` é normalizado antes de sua norma; os coeficientes escalados usam
expoentes separados para evitar overflow intermediário. NaN, singularidade
e resultado não representável retornam status. Comparações de faixa
são feitas depois dos testes de finitude, sem pressupor curto-circuito
Fortran. Toda proposta volta ao resíduo conservativo original, ao domínio
físico e à margem não linear **g>=0**, sem tolerância negativa. A coluna
`direction_objective` em `active_set.csv` distingue `L2` e `LINF`.
`accepted=T` significa descida Newton, não conclusão do passo temporal.

## Experimentos causais no mesmo estado e dt

| Replay a partir de Uold | LINE_SEARCH | ACTIVE_SET com gatilho/L∞ |
|---|---|---|
| SS1 terminal v1; dt=1,77978515625e-7 s | status6000, 9 Newton | status0, 1 passo, 1 Newton, 3 montagens de resíduo |
| SS3 terminal v1; dt=1,3286025e-7 s | status6000, 1 Newton | status0, 1 passo, 4 Newton, 42 montagens, 14 soluções lineares |
| SS3 terminal v2; dt=1,9131876e-7 s | status6000, 1 Newton | status6000, 16 Newton; passo não concluído |

SS3 v1 snapshot SHA256
`683f5c44e2acad1de051e130040613156b0b666dddb18a2007b80430d6d2b1fe`.
No seu ponto de linearização, a norma inicial é `4,49997595858e-8`.
L2 com lambda1 dá norma `7,48628874208e-9`, mas margem
`-5,55111512313e-17 m/s`: é recusada. L∞ com lambda1 dá norma
`7,46514076553e-9`, margem exatamente zero e status físico0. A tolerância
é `1e-8`. O teste integrado parte do mesmo Uold/dt e verifica também
massas positivas, igualdade dos fluxos de fase e identidade entre balanço
integral e resíduo integrado. Não se aceita resíduo alterado ou zerado.

Na célula69 desse snapshot, theta=`0,9012386672268461 rad` desde a
horizontal: ela não é vertical. No estado aceito anterior, em SI,
`[alpha,rhol,rhog,p,ul,ug]=[0.6382716659196366,997.10664660672,
1.94170976434467,228440.212065386,0,0.624559711373577]`.
Na linearização, os valores são
`[0.6382716662903444,997.106646607824,1.94170978540585,
228440.214543213,0,0.624559711664570]`.
A condição2 do Jacobiano escalado é `1,00058516378`; ±1 ulp de massa
gasosa ou momento (`2,16840434497e-19` adimensional) muda a margem em
aproximadamente `5,55111512313e-17 m/s`. A sensibilidade do sinal nesse
ponto é de arredondamento; não se transformou a margem negativa em zero.

SS3 v2 snapshot SHA256
`a92b9f0652254d6366697e73e3bdf4773cfc6d402b190a601ee12f8b268ea7c1`.
Norma inicial `1,077897799327e-8`; margem `8,32667268469e-17 m/s`;
condição2 `1,00084271214`. O limite ótimo **linear local** em norma∞ é
`1,07497553507e-8`, acima da tolerância. A proposta L∞ integral produz
`1,07497703859e-8` e margem negativa de uma ordem de ulp. Os replays
L2/L∞ e o passo completo não convergem. Isso limita as alternativas
examinadas nesse estado/dt; não prova inexistência global nem necessidade
inevitável de contrafluxo. Não se alterou Bhagwat ou a física de frente.

## Smokes anteriores, separados da regressão

| Versão | SS1 solicitado100s | SS3 solicitado110s |
|---|---|---|
| v1, commit3940bb7, L2 tardio | status6100/MAX_NEWTON em78,5717914962s; CPU99,8564s; wall711,586s | status6000 em97,1318529144s; CPU42,5804s; wall43,6881s |
| v2, gatilho cedo/L2 | TIMEOUT externo900s; último registro completo78,5804787653s; CPU registrada497,598s (limite inferior) | status6000 em97,1318527417s; CPU32,1117s; wall32,2177s |

No SS1 v2, a última linha completa tem 2195 passos aceitos, 277 rejeições
e dt=`3,14008230439e-6 s`. O avanço de milissegundos com dt pequeno não
qualifica a trajetória. O arquivo de rejeição existente no TIMEOUT é uma
rejeição **anterior**, não um estado terminal capturado no instante da
interrupção. Não existe `run_status.txt` final desse caso.

Os seis balanços finais registrados, em kg, na ordem
`[l,g,l_pipe,g_pipe,l_riser,g_riser]`, são:

- v1 SS1: `[-4.0473e-11,4.3195e-15,7.7023e-12,3.2064e-12,-4.8175e-11,-3.2020e-12]`.
- v1 SS3: `[-1.7656e-8,2.6379e-11,1.1914e-10,-1.5058e-11,-1.7775e-8,4.1436e-11]`.
- v2 SS1, último registro: `[-6.0262e-6,9.3305e-10,3.6293e-8,5.7996e-10,-6.0625e-6,3.5309e-10]`.
- v2 SS3: `[3.0269e-12,-2.3204e-14,1.3385e-10,-5.9221e-12,-1.3087e-10,5.8989e-12]`.

As discrepâncias entre fluxos mássicos esquerdo/direito da junção são
exatamente zero nos registros aceitos. Os balanços são resíduos reais
acumulados, não valores impostos; os JSONs preservam magnitudes completas,
inventários e erro da identidade integral. Nenhum desses smokes atingiu
o tempo pedido ou demonstrou ciclos acomodados/validação experimental.

## Verificação e proveniência

`make test FC=gfortran LIBS='-llapack -lblas'`: 15/15 programas,
incluindo 27 checks matemáticos e 22 checks de replay, com traps de
depuração. Python da CI: 19/19. As verificações cobrem KKT, norma dual1,
singularidade, escala, sinais, NaN e overflow/representabilidade. Os
gradientes offline dos snapshots originais foram repetidos com fatores
0,25/0,5/2/4, mantendo os mesmos resultados finais (README/JSON anterior).

O diretório `minimax/` contém manifestos de fontes, replays v1/v2, resumos
dos smokes históricos e certificado local de geometria. Binários e PDFs
bibliográficos permanecem fora do Git. Os novos smokes v3 deverão usar
build fresco de um commit exato, manifesto por caso, CPU, wall e dt real;
não foram substituídos por estes replays.


Arquivo local dos dados históricos: `C:/Users/leonardoab/Documents/rotordin_frontend/artifacts/primsf-qualification2-20261119/active_set_v1_v2_replays_20261119.zip`,
16.582.095 bytes, SHA256 `0bff87a8a2f0dfab759530a0fa1b9019bb25c408d1e704eea6bf946c83a10f58`.
CRC de todos os membros verificado; zero membros ELF/PE ou PDFs
bibliográficos. O manifesto de membros e exclusões está em
`minimax/historical_archive.json`. O pacote preserva séries v1/v2 e replays;
não contém uma trajetória v3 completada.


No freeze final a179e51, a verificação fresca completa passou15 programas
Fortran,76 testes Python,21 verificações CLI e5 de telemetria. Os19
testes Python mencionados acima eram apenas o subconjunto da CI histórica.
O terminal do smoke SS3 v3 e seu replay independente80d estão em
[minimax/v3_ss3/README.md](minimax/v3_ss3/README.md); esse caso não atingiu110s.
