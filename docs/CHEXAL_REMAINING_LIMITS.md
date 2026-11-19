# CHEXAL_1992: limites de continuação e busca

Revisão pontual de 18/11/2026, sobre `b0d33cb25ab2480b3eb0f5d37c18edc340b0f122`
mais a correção do núcleo descrita abaixo. Complementa
`CHEXAL_1992_IMPLEMENTATION.md`; não substitui os relatórios das campanhas.

## Correções verificadas nesta revisão

A enumeração de velocidades podia descartar um status `SINGULAR` ou `NUMERIC`
obtido no refinamento ou na revalidação do ramo e terminar em `NO_ROOT`. Isso
permitia tentar uma raiz cocorrente depois de uma falha da solução CC. Esses
status agora são propagados. Rejeições esperadas de domínio, ausência de raiz ou
CCFL continuam distinguíveis de falhas numéricas.

Na referência ar–água do teste core, `jg=0,1 m/s`, a recuperação do estado na
dobra, sem dicas de velocidade, devolvia `OK` em HIGH e LOW, mas com
`jl=-0,443604697759832` e `jg=-0,0624943184289925 m/s`: uma solução cocorrente
descendente. Agora ambos retornam `SINGULAR (65)`, sem trocar o ramo.

A localização de CCFL também podia confundir o desaparecimento de um máximo
entre os nós da malha com uma dobra. O resultado agora precisa satisfazer:

- `0 < alpha_CCFL < 1`;
- resíduo absoluto menor que `2048 epsilon` vezes a maior magnitude dos três
  termos `alpha C0(jl+jg)`, `alpha Vgj` e `jg`;
- valores vizinhos, a distância `epsilon^(1/3) min(alpha,1-alpha)`, que não
  excedam o máximo além da mesma tolerância.

As desigualdades unilaterais permitem máximo em quina da correlação. Elas são
checagens locais de aceitação, sem constituir prova global da dobra ou de sua
unicidade. A malha, as fórmulas e o domínio físico não foram ampliados.

Com `rho_l=998`, `rho_g=1,2 kg/m³`, `D=0,05254 m`, `g=9,81 m/s²`,
`mu_l=0,001`, `mu_g=0,000018 Pa.s` e `sigma=0,072 N/m`:

| jg (m/s) | Antes | Depois |
|---|---|---|
| 0,1 | CCFL com resíduo relativo `6,01e-16` | Mesmo resultado |
| 0,01 | CCFL com resíduo relativo `-7,63e-16` | Mesmo resultado |
| 1e-4 | `OK`, resíduo relativo `0,22523` | `NUMERIC (66)` |
| 1e-6, 1e-8, 1e-12 | `OK`, alpha=0, resíduo relativo -1 | `NUMERIC (66)` |

O teste isolado do núcleo passou **147 verificações**, incluindo oito novas
regressões, com GNU Fortran 13.3.0, `-std=f2018 -O0 -fcheck=all
-ffpe-trap=invalid,zero,overflow`. A evidência antes/depois está em
`validation/diagnostics/chexal_ccfl_status_20261118/`. `NUMERIC` nesses exemplos
identifica insuficiência da busca atual; não afirma inexistência física de CCFL.

A suíte integrada executou 3447 verificações em 12 programas Fortran, 19 testes
Python da CI e build release; comandos e hashes constam na mesma pasta.

## Limites materiais remanescentes

**Completude das raízes.** As malhas finitas, os refinamentos de brackets e os
máximos locais detectam as raízes exercitadas pelas referências. Não certificam
que todo par muito próximo, tangência entre amostras ou primeira dobra foi
encontrado. O limite superior analítico do intervalo de velocidades não transforma
a amostragem desse intervalo em prova de unicidade. Uma única raiz detectada pode
ser aceita; essa qualificação vale também para AUTO após procurar alternativas
CC. Ampliar a validade da busca requer algoritmo adicional e referências próprias.

**Continuação.** HIGH/LOW selecionam os ramos físicos definidos na p.268 da fonte;
dicas permitem continuação local e são revalidadas. Não há identidade global de
raiz persistida em toda a evolução, nem regra física implementada para atravessar
uma dobra e escolher o estado depois dela. Na singularidade, o contrato é parar
com status. Isso não equivale a um modelo dinâmico de transição entre ramos.

**Fases ausentes.** Inventário conservativo exatamente nulo e alpha=0/1 têm
tratamento explícito. Zero fluxo superficial, sozinho, não determina inventário
nulo. A velocidade armazenada da fase ausente é convenção; não é uma raiz
bifásica adicional. O Jacobiano bifásico completo não é regularizado nas fases
ausentes. A verificação desses limites não demonstra uma evolução transiente
com aparecimento/desaparecimento de fases e troca automática de modelo.

**Derivadas e condicionamento.** As derivadas mantêm a peça da correlação e
verificam convergência local, inclusive perto de velocidade não nula minúscula.
Nas quinas e na dobra podem não existir derivadas bilaterais ou inversão regular.
Status de falha são preservados; não há suavização de coeficientes, clipping ou
espectro artificial para forçar a continuação.

**Domínio e validação.** A contracorrente autorizada permanece vertical, com gás
subindo e líquido descendo, conforme Chexal1992 p.262; horizontal é cocorrente
(p.264). A catenária atual não satisfaz essa condição. Traços gasosos negativos
MATCHED/Roe e a discretização de inclinação variável exigem análise própria.
Referências de alta precisão e testes de operadores verificam implementação;
não demonstram reprodução experimental de OSC/SS3, ciclos acomodados ou
convergência da campanha. Os resultados dessas tarefas são registrados
separadamente, preservando a divergência documental SS3 de 26 s versus 76 s.
