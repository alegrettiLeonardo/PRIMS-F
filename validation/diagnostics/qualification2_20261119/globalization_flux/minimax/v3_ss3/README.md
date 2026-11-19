# SS3 v3: terminal no operador a179e51

O smoke solicitou110s com `ACTIVE_SET`, usando build fresco do commit
`a179e51d9d63e8aeed91a94d829b5b483d107864`, binário SHA256
`02ce51462ec3c954cec9364b17444b7279be72ec150f87e307f1b5fc8dd0a4c3`.
Terminou com status6000/`NEWTON_GLOBALIZATION_FAILURE` em
`97.131852741709523s`, após812 passos aceitos e24 rejeições. CPU do
processo44.711895s, wall44.967827s. Não houve TIMEOUT. O tempo atingido
é igual ao v2, mas o snapshot terminal é diferente; nenhum dado anterior
é reutilizado como se fosse o novo terminal.

Snapshot SHA256 `9a85ce3c5769d1dd11a9e8860b79a2936db00e42c6922771727f9e6102778139`.
O sidecar registra Newton11, backtrack18, lambda=`2^-17`, proposta1017,
célula69. A norma dessa proposta é indisponível porque a montagem foi
recusada. O `huge()` do CSV bruto é um sentinela, não um resíduo calculado.

## Mesmo estado e direções alternativas

O probe ligado aos objetos a179e51 reproduziu bit a bit o vetor de resíduo
salvo na linearização. Norma original=`1.075081558190227e-8`, acima de
`tol_res=1e-8`. Margem double no endpoint cocorrente=`0`. A célula69
tem theta=`0.9012386672268461rad` desde a horizontal e usa Bhagwat ID3;
não é uma seção vertical nem um teste de fechamento contracorrente.

| Direção, lambda1 | Status físico | Margem double (m/s) | Norma conservativa original |
|---|---|---:|---:|
| Newton |1017/CELL69|−1.60219951029e-6|indisponível|
| L2 |0|+1.38777878078e-16|1.07824340961e-8|
| L∞ |0|+2.77555756156e-17|1.07497558552e-8|

A direção L∞ integral é fisicamente admitida pelo operador e reduz a
norma, mas permanece acima da tolerância e não satisfaz a redução Armijo
exigida. L2 integral aumenta a norma. As18 propostas de cada direção
estão no CSV; não se relaxou tolerância ou margem para aceitá-las.

A condição2 do Jacobiano escalado é`1.000842712204`. O mínimo linear
local em norma∞ para a única desigualdade é`1.074975535314e-8`, também
acima da tolerância. Esse limite se refere ao subproblema linear nesse
estado/dt. Ele não prova inexistência global de soluções nem necessidade
inevitável de contrafluxo. A geometria double independente usou a API
física congelada69813ee, cujo resíduo foi novamente verificado bit a bit
contra este snapshot; o replay das direções de produção usa objetos a179.

## Verificação independente com80 dígitos

Reutilizou-se `scripts/replay_root_bounds.py`, com a transcrição primária
de Azevedo, Baliño e Burr(2015), p.121, Eqs.(116)–(128), e solução de
Colebrook por bisseção. As referências completas estão em
`docs/BHAGWAT_SOURCE_STATUS.md` e no README anterior. Os bits de U e
theta foram exportados diretamente do probe Fortran. Cada estado foi
avaliado usando seus literais serializados e os valores binary64 exatos.
As propriedades SI nominais são decimais nas duas avaliações; não se
afirma uma prova por aritmética intervalar de todas as operações do solver.

O limite suficiente foi testado nas propostas lambda1,1/2 e2^-17 de
Newton/L2/L∞ e nos dois estados de referência. São11 estados e22
avaliações, não22 estados independentes.

- Nas três propostas Newton, o limite inferior de `F=R/alpha` é positivo
  com folga: `[1.60219951031e-6,8.01099755082e-7,1.22237774523e-11]m/s`.
  As hipóteses do limite global cocorrente são satisfeitas e a folga
  excede`128eps*max(omega_u,Ud,ug_max)`. A exclusão refere-se a essas
  propostas e a essa lei.
- No estado aceito anterior e na linearização, o avaliador binary64
  exato com propriedades nominais dá `F(1)=-1.94919607e-17` e
  `-3.66006763e-17m/s`. As raízes nominais têm ul da ordem de1e-19m/s;
  isso fica abaixo da resolução prática do sinal double e é classificado
  como marginal/continuação não qualificada.
- As propostas L2/L∞ examinadas também ficam nessa escala marginal.
  Variações de sinal double de um ULP não são usadas como prova matemática
  de exclusão. O bloqueio local relatado permanece mesmo na proposta L∞
  com margem double positiva, pois sua norma supera a tolerância.

Em SI, a referência nominal da linearização fornece
`[alpha,rhol,rhog,p,ug_max]=[0.6382716624417565,997.1066466071015,
1.941709771628229,228440.2129222895,0.6245597085796658]`.
O estado aceito anterior fornece
`[0.6382716651810471,997.1066466045193,1.941709722374810,
228440.2071276740,0.6245597107938116]`.

## Conservação e limite da conclusão

Os seis balanços finais reais em kg, na ordem
`[l,g,l_pipe,g_pipe,l_riser,g_riser]`, são
`[3.02691205e-12,-2.32036612e-14,1.33844935e-10,-5.92212390e-12,
-1.30853550e-10,5.89894800e-12]`.
O maior erro entre balanço e resíduo integrado é`9.50e-14kg`.
As diferenças entre fluxos esquerdo/direito da junção são exatamente
zero para ambas as fases nos registros aceitos; os inventários permanecem
positivos. Isso verifica a conservação dos passos aceitos, não valida a
trajetória incompleta nem experimentalmente o regime SS3.

O resultado documentado é uma tentativa numérica local limitada, com
propostas Newton robustamente fora do ramo e alternativas ativas
marginais que não atingem a norma requerida. O operador, a física e o
controlador temporal congelados não foram alterados por esse diagnóstico.
