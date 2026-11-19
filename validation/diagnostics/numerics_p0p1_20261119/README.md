# P0/P1: contornos e recusa SS1 na base 9f2ed21

Auditoria somente de leitura do solver, em 2026-11-19. O probe usa os objetos
release de `8a85e37`; `git diff 8a85e37 9f2ed21 -- scr` é vazio. Fontes,
binários, snapshots e comandos estão em `provenance.json`. As alterações
concorrentes de status, histórico de Roe e instrumentação no worktree não
entraram neste replay. Não houve clipping, projeção de gás, troca de coeficientes
ou reparo de espectros.

## P0: características e dados de contorno

Foram calculados 36 espectros: estados interiores, estados efetivamente
construídos por `boundary_flux`, vizinhos/mapeamento/traços da junção, tanto no
último estado aceito como na linearização salva, em SS1 e SS3. Todos têm parte
imaginária exatamente zero na saída do LAPACK; maior resíduo relativo
`A R - R Lambda` = 3,722e-18. Os 16 espectros de contorno mantêm duas
características entrantes na entrada (`lambda>0`) e uma na saída (`lambda<0`).

Autovalores SI dos estados de face impostos na linearização salva:

| Caso | Face | Autovalores (m/s) | Entrantes |
|---|---|---|---:|
| SS1 | Entrada | −246,681985; 1,615399; 247,888413 | 2 |
| SS1 | Separador | −33,392533; 0,931228; 34,369022 | 1 |
| SS3 | Entrada | −312,232972; 1,411841; 313,365066 | 2 |
| SS3 | Separador | −34,851999; 0,693482; 35,244741 | 1 |

As duas condições de massa na entrada são independentes sobre os dois modos
entrantes: os valores singulares da matriz projetada, com linhas normalizadas,
são 0,761826/0,159402 em SS1 e 0,532405/0,262825 em SS3. A condição de pressão
na saída tem projeção não nula no modo entrante. Isto verifica a contagem e a
independência local nestes estados; não prova estabilidade global nem ausência
de reflexão da implementação por extrapolação.

O código fixa `Gl` e `Gg` na entrada e usa a pressão interior. Na saída fixa
`Ps`, retém alfa e momento total interiores e recupera o estado sob o fechamento
selecionado. Ghosts são auxiliares; o fluxo físico de contorno é calculado à
parte. Ambos os contornos passam no replay SS1/SS3. A recusa terminal SS1 é
interna, na face53/54, portanto não há evidência de que uma condição extra no
separador deva ser retirada para resolver esse snapshot.

[Baliño, Burr e Nemoto (2010)](https://doi.org/10.1016/j.ijmultiphaseflow.2010.04.003)
prescrevem `Ql0`, `mg0` e `Ps` (p.647, seção2.3), usando gás ideal isotérmico,
líquido incompressível e momento sem inércia. Seu sistema NPW tem uma única
característica finita e duas infinitas (p.647, Eqs27–31); essa contagem não pode
ser transplantada para as três equações evolutivas inerciais deste código.

Confronto específico com o experimento e a referência:

- `Ps` é absoluto: SS1 usa208000 Pa, SS3 usa201000 Pa. Correspondem às médias
  de2,08 e2,01 bar a da Tabela3/p.653, casos2 e11. No experimento as dispersões
  são0,03 e0,07 bar; uma pressão exatamente fixa é a idealização do modelo,
  não a reprodução da série medida. A instalação tem separador, controle de
  nível líquido e válvula atuada; o modelo atual não representa sua dinâmica
  (descrição do equipamento, p.652).
- A conversão de gás STANDARD usa101300 Pa/293,15 K e conserva a massa
  imposta, em acordo com a estrutura da Eq78/p.649. Não se impõe `jg` constante
  em todas as pressões locais.
- Para a água, o código fixa massa e admite compressibilidade. Assim,
  `jl=Gl/rhol(P)` difere do `Ql0/A` constante da referência: na entrada da
  linearização, −0,0077699% em SS1 e −0,0058246% em SS3. Essa diferença foi
  medida e mantida; não foi corrigida alterando os decks.
- Os decks D mantêm EOS barotrópica com `Cg=343 m/s` e Bhagwat; Baliño2010
  usa EOS isotérmica e Bendiksen. SS1 também preserva `jl0=0,530 m/s` e
  `Tgas=20,1 C` da dissertação; a Tabela3/caso2 informa0,528±0,008 m/s e
  19,7±0,1 C. SS3 usa0,28/0,48 m/s e19,8 C, compatíveis com o caso11.
  A divergência de período SS3, dissertação26 s versus tabela76 s, permanece.

## P1: replay mínimo novo da recusa SS1

`replay_ss1_minimum.py` foi executado nesta etapa, ligando o probe existente aos
objetos congelados e relendo o snapshot novo de `v1/runs/D_SS1_fixed`. A pasta
`SS1_minimal_replay` contém log e proveniência próprios. A diferença máxima
para o vetor residual salvo é zero. Não se reutilizou um log antigo como uma
execução nova.

No tempo aceito78,571788476739314 s, viscosidade de Roe, todas as células,
faces e contornos são admissíveis. Na linearização salva o residual é
2,9764063309210322e-7. As 18 propostas, de lambda1 até2^-17, preservam as
células admissíveis, mas a primeira recusa ocorre na face53, ao recuperar o
traço direito de MATCHED: status12 → face1012 → residual3012. O término
real6000 significa falha de line search. Normas de propostas cuja montagem
falhou são indisponíveis; o sentinela `huge` não é uma medição de residual.

No ponto salvo, o fluxo gasoso é a diferença entre parcelas convectiva e
dissipativa de aproximadamente1,2070213e-4 kg/s, deixando+1,07438e-15 kg/s.
As células vizinhas têm alfa0,55998/0,01867 e velocidades gasosas positivas
0,05308/0,89901 m/s. A proposta lambda1 produz `jg*=-6,0549846e-9 m/s`, ou
−3,14346e-11 kg/s. Mesmo lambda2^-17 produz `jg*=-3,2402370e-14 m/s`.
O cancelamento é severo nos menores passos; a proposta integral também é
negativa, portanto maior precisão sozinha não resolve a travessia.

O traço riser salvo tem alfa2,31924e-13 e espectro real aproximadamente
`[-1497,4705;0,892312;1498,5295] m/s`: aproxima-se do limite acústico líquido.
Ainda não existe mudança de contagem das características enquanto alfa é
positivo. Não se calculou um espectro artificial em alfa0.

Como contraprova limitada, foi resolvido um sistema de acoplamento
**linearizado** com uma onda entrando no pipe, duas no riser e continuidade
de pressão/fluxos de fase. A matriz3x3 tem posto completo; seu resultado SS1
tem `Gg=-0,00241723 kg/m2/s` e massa gasosa direita negativa. Este candidato
não foi aplicado. Ele omite curvas de onda não lineares e a integração das
fontes entre centros e face; não prova inexistência de um Riemann admissível.
Mostra que substituir o fluxo mapeado por esse matching linear não constitui
uma correção admissível demonstrada. O resultado SS3 correspondente é apenas
diagnóstico da junção; sua recusa de célula é investigada separadamente.

## O que a fonte permite corrigir, e o que permanece aberto

A proposta negativa está fora do domínio cocorrente implementado de Bhagwat.
Isso é uma recusa de domínio do traço numérico; não demonstra que o escoamento
físico tenha gás reverso nem prova ausência global de solução de interface.
Trocar o fechamento por Chexal não autoriza contracorrente horizontal.

Baliño2010 trata esse evento com um estado físico adicional. Na condição
`x>0`, a frente líquida no pipe satisfaz `jgb=0`,
`dx/dt=(Ql0/A-jlb)/alpha_p`, uma equação de compressão de gás e
`Pb=Pg+rhol*g*x*sin(beta)` (Eqs1–4/p.645). As Eqs8–14/p.646 definem a troca
reversível entre penetração e acumulação, preservando pressão, derivada de
pressão e velocidade superficial total. As Eqs39–42/p.647 acoplam os
subsistemas; Eqs59–67/p.648 discretizam os dois estados. A passagem física é
de gás positivo para gás bloqueado enquanto a frente/inventário evolui;
não é a aceitação de uma raiz negativa em alfa.

Duas linhas de trabalho são justificáveis, mas nenhuma é uma correção local
já validada: (a) solver de interface não linear com curvas de onda e transição
de fase ausente matematicamente derivada para os fechamentos atuais; (b) modelo
alternativo explícito da frente/compressão de Baliño, com suas hipóteses e
inventários próprios. A segunda altera a física e deve ser selecionável,
com baseline preservada. Importar somente `jgb=0` de sua Eq1, sem a frente e
seus balanços, seria uma projeção incompleta. Uma extrapolação do espectro ou
da correlação para fabricar continuidade também não é sustentada.

Pesquisa primária complementar: Evje e Flåtten,
[On the Wave Structure of Two-Phase Flow Models](https://doi.org/10.1137/050633482),
relaciona a estrutura de ondas às derivadas da lei de slip. Foi consultado
somente o resumo oficial; ele não fornece aqui um solver pronto para a
junção estratificado/Bhagwat. A formulação inercial posterior de Baliño,
[DOI10.1016/j.jocs.2013.08.006](https://doi.org/10.1016/j.jocs.2013.08.006),
foi localizada, mas o texto integral não foi obtido nesta etapa. Nenhuma
equação não lida desses trabalhos foi implementada.

O teste justificável antes de mudar a junção deve verificar ondas/traços
admissíveis, limites de fase ausente, conservação por fase e pressão comum,
além de ida/volta através do evento e do balanço da frente, se introduzida.
Para esta base, a ação mínima sustentada é preservar e explicitar o status
de recusa, com captura da proposta e do estado aceito. Não há autorização
física demonstrada para forçar a passagem alterando apenas o fluxo gasoso.

Os PDFs permanecem fora do Git. O PDF local de Baliño tem SHA256
`c60944a348bf516f20ed7fd965cff7aa6b8417b07ef3b4b9433c64d28db8a3a5`.
Páginas PDF3/4/5/6/10/11 correspondem às impressas645/646/647/648/652/653.


A revisão suplementar e o alcance da exclusão SS3 posterior constam em
[ALTERNATIVES_REVIEW.md](ALTERNATIVES_REVIEW.md). A base dos espectros e do replay
acima permanece 9f2ed21, sem reatribuir resultados a fontes posteriores.

## Nova execução terminal em 69813ee

Conferência somente de leitura dos manifests `v1_2400/runs/D_SS1` e `D_SS3`.
Seal, plano, decks, 23 fontes e binário foram verificados por SHA. O binário
é `7442291fd6cdc311627e79500dc5a4804813b42c7bf4daecf88b4b8b2322288a`.
O plano preserva seu rótulo original `PREPARED_UNSEALED`; o `seal.json`
separado existe e seus hashes coincidem com os manifests e arquivos executados.

| Caso | Tempo aceito (s) | Aceitos/rejeitados | Proposta terminal | Iteração | lambda |
|---|---:|---:|---|---:|---:|
| SS1 | 78,571788476739314 | 662/24 | 3012, FACE, índice53 | 13 | 2^-17 |
| SS3 | 97,131852338700099 | 809/24 | 1017, CELL, índice69 | 8 | 2^-17 |

Os dois processos retornam3 e status6000 `NEWTON_GLOBALIZATION_FAILURE`.
Sidecar e última linha de `rejected_proposals.csv` coincidem: Jacobiano
completo, backtrack18, lambda real2^-17 e norma da proposta indisponível
(`NaN`, registrada como null no resumo). As últimas18 propostas de cada
caso têm o mesmo operador/índice/status terminal. Todos os2353/2288 registros
de recusas passaram nos checks de enums, índices e normas disponíveis.
Em FACE o índice designa a célula esquerda; em CELL, a célula recuperada.

Os snapshots completos são **byte-idênticos** aos de Burr15-v1, inclusive
estado aceito, linearização e residual: SS1 SHA
`1e9629e5293c8ee55ae324b6094b0a73d46aee47b6673dda2e5ee456917f27dc`;
SS3 SHA `e2cafe7245ac084fe8b8db764d4d19bbc6d4dc6e0f891ec4a58ad68f054fac3c`.
Não houve avanço físico adicional. A diferença SS3 é diagnóstica:
`1017=1000+10+NO_PHYSICAL_ROOT(7)` substitui a recusa genérica1014, conforme
o [certificado cocorrente](../../../docs/SS3_COCURRENT_BOUND_20261119.md).
Os números e o replay independente desse certificado permanecem na sua
evidência própria; os espectros anteriores não foram reatribuídos a69813ee.

Comandos e hashes desta conferência: `check_new_terminal_manifests.py` e
`terminal_69813ee_summary.json`. Nenhum solver novo foi executado pelo checker.

Dados brutos dos probes e log de compilação preservados fora do Git em `PRIMS-F-numerics-results-20261119/operator_evidence/p0p1`.
