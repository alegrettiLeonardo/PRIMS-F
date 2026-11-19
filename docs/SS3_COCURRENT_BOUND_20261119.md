# Exclusão de raiz cocorrente em SS3

Veja [a transcrição de Bhagwat–Ghajar](BHAGWAT_SOURCE_STATUS.md) para a
correlação efetivamente implementada. A referência é Azevedo, Baliño e
Burr (2015), p. 121, Eqs. (116)–(128). O certificado usa o ramo ascendente
dessa versão, com \(0\le\theta\le\pi/2\), \(0<\alpha<1\) e
\(\rho_l>\rho_g>0\). Ângulo descendente pode cruzar uma troca em \(C_4\).

Fixadas as massas de fase \(m_l,m_g\) e o momento \(M>0\), o intervalo
cocorrente completo é parametrizado pela contribuição gasosa
\(q=m_g u_g/M\). Então

\[
u_l=(1-q)M/m_l,\quad u_g=qM/m_g,\quad
j(q)=(1-q)(1-\alpha)M/m_l+q\alpha M/m_g,\qquad q\in[0,1].
\]

O resíduo dimensional dividido por \(\alpha\) é
\(F(q)=C_0(q)j(q)+U_d-u_g(q)\). \(U_d\) é constante neste intervalo
porque \(C_4=1\) no ramo ascendente e as propriedades são fixas. Seja
\(A=2-(\rho_g/\rho_l)^2\) e
\(B=[(1+(\rho_g/\rho_l)^2\cos\theta)/(1+\cos\theta)]^{(1-\alpha)/5}\).
A Eq. (116) é a média ponderada de \(A\) e \(B+C_{0,1}\) com peso
crescente \(w=Re^2/(10^6+Re^2)\). Vale \(A>1\ge B\ge 2^{-1/5}\).

- Na Eq. (118), \(H(y)=y+4\log_{10}(\varepsilon/(3.7D)+1.256y/Re)\)
  cresce com \(y\) e decresce com \(Re\). Se \(H(1/3;Re(0))<0\)
  e \(j\) cresce, então \(y>1/3\) e \(C_{0,1}>-0,6\) em todo o
  intervalo. Isso dá \(C_0>0\).
- Escolha \(q_c=\min(1/2,U_d/(2u_{g,\max}))\). Em \(q\le q_c\),
  \(F\ge U_d-q_c u_{g,\max}>0\).
- Se \(H(1;Re(q_c))<0\), então \(C_{0,1}\ge0\) para \(q\ge q_c\).
  Como \(C_{0,1}(1)=0\) e \(A\ge B\), resulta
  \(C_0(q)\ge C_0(1)\), logo
  \(F(q)\ge U_d+\min(0,\alpha C_0(1)-1)u_{g,\max}\).

O operador declara `NO_PHYSICAL_ROOT` quando as hipóteses são satisfeitas
e o limite é positivo acima da margem de arredondamento. Um momento
não positivo também exclui um estado bifásico Bhagwat com \(j>0\).
Falhas de enquadramento ou convergência sem certificado ficam com
`ROOT_SEARCH_FAILURE` (4, alias de compatibilidade `PRIMITIVE_NO_ROOT`).
Uma recuperação final inválida usa `INADMISSIBLE_ROOT` (8). O status
`NEWTON_GLOBALIZATION_FAILURE` (6000) pertence ao solver de tempo e
carrega separadamente o status da proposta e seu local.

## Célula 69 do snapshot Burr15

O [replay do caso](../validation/diagnostics/burr15_rejection_replay_20261119/)
reproduziu exatamente o resíduo salvo. O novo probe serializa os blocos do
Jacobiano, o resíduo por célula e os estados conservativos para cada
\(\lambda\). A equação de estado e o limite acima foram reavaliados com
80 dígitos, usando bisseção independente da Eq. (118).

Na linearização, \(p=228440,194157\) Pa, \(\alpha=0,638271663241\) e
\(F(1)=-2,46867\times10^{-12}\) m/s. Há raiz cocorrente com
\(u_l=2,19740\times10^{-14}\) m/s. Em \(\lambda=1\), o limite
inferior global é \(+5,98498\times10^{-7}\) m/s. Em \(2^{-17}\),
é \(+2,09752\times10^{-12}\) m/s. As 18 propostas nos dois modos de
Jacobiano foram excluídas. A inclinação era \(0,901239\) rad e
\(\mathrm{d}F/\mathrm{d}q\) à esquerda de 1 era aproximadamente
\(-0,241098\) m/s. O resíduo conservativo da proposta fica indisponível
quando a recuperação de primitivas impede sua montagem.

O certificado exclui as raízes para esses estados e esse fechamento.
Um Newton restrito poderia tentar outra direção, mas isso por si só
não especifica a evolução física após a fronteira cocorrente. Retorno
e contrafluxo inclinado exigem um fechamento qualificado ou um modelo
diferente. O limite evita trocar silenciosamente a lei ou atribuir uma
exclusão demonstrada ao limite de iterações.

Reprodução com a árvore congelada, apontando para os objetos do build
informado no manifesto:

```sh
python3 scripts/replay_root_bounds.py --objects DIR_OBJ --runs BURR15_V1_RUNS \
  --output REPLAY_NOVO
```

O conjunto desta execução está em
`PRIMS-F-numerics-results-20261119/replay_burr15_base`, com SHA de cada
snapshot, probe, objeto compilado, log, Jacobiano e proposta.

Na nova campanha, ambos os snapshots terminais foram reproduzidos
byte a byte. Com os objetos do release `69813ee`, o replay em
`PRIMS-F-numerics-results-20261119/replay_69813ee` usa os diretórios de
cenário sem sufixo (`--suffix=`). Confirme a procedência no manifesto
daquele replay antes de comparar os status das propostas.
