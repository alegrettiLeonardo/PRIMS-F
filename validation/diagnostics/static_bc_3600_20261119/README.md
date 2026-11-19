# Seis sensibilidades estáticas de contorno, 3600 s

Operador congelado `69813ee70551fa4a3f075919aeafc45de9693f11`, binário
`7442291fd6cdc311627e79500dc5a4804813b42c7bf4daecf88b4b8b2322288a`.
As seis simulações terminaram 3600 s com zero rejeições. Foram executadas
sequencialmente, um worker/thread, conservando EOS, geometria, malha,
passo, fechamento e critérios de Newton. Ps/jg/jl mudam um eixo por vez;
as escalas da Tabela 3/caso 13 são desvios de sinais sem definição estatística
suficiente para tratá-las como erro, sigma ou intervalo de confiança.

[Método e fontes](../../../docs/STATIC_BC_AND_REGIME_CRITERIA_20261119.md),
[índice/hashes](campaign_index.json), [métricas completas](metrics_3600.json).
Os decks copiados são registros exatos dos casos executados; seus caminhos
de saída apontam para a coorte preservada. O gerador permite preparar uma
nova pasta sem sobrescrever essa evidência.

A pressão é `p_junction_Pa`. Os valores da tabela usam os mesmos dez
intervalos completos de cada execução; cada canal conserva cobertura,
janela e ponderação temporal no JSON. Pico a pico por ciclo não é a definição
experimental “Max 0,69 bar”. O detector usa histerese 20/80 e piso 100 Pa;
esta medição numérica não classifica golfada nem valida o experimento.

| Variante | Período médio (s) | Excursão média por ciclo (Pa) | Estado da comparação cíclica |
|---|---:|---:|---|
| Ps- | 23.13372656 | 22712.364 | Acomodado |
| Ps+ | 23.06398387 | 15379.414 | Provisório |
| jg- | 23.99813365 | 34502.529 | Acomodado |
| jg+ | n/a | n/a | Sem ciclos resolvidos |
| jl- | 23.18031925 | 20977.295 | Acomodado |
| jl+ | 23.03697038 | 18568.062 | Acomodado |

Referência nominal pela mesma análise: T = 23.10861760 s,
excursão média = 19818.034 Pa e pressão média = 266193.075 Pa.
Ps+ ainda apresenta crescimento de amplitude e não satisfaz acomodação;
baixa dispersão do período não substitui esse teste. Em jg+, o auxiliar de
picos encontra apenas oscilação residual decrescente abaixo do piso:
pico a pico 0,43674 Pa em 3355-3579 s, semiamplitude final 0,14128 Pa.
Não reportamos seu período como qualificado nem interpretamos ausência de
ciclos como prova de estabilidade física.

Todos seis passaram o contrato aritmético P2 e os seis balanços globais
(cada fase no sistema, oleoduto e riser), tanto normalizados por inventário
quanto por massa processada, no critério 1e-6. O pior erro/inventário é
2,6941e-7 em Ps+. Para os cinco casos com ciclos, os seis balanços de cada um
dos dez intervalos passaram; jg+ registra zero intervalos e métricas null,
sem aprovação vazia. As quatro amostras de contorno por frame têm cobertura
completa, status 0, contagens características 2/2/1/1 e nenhuma velocidade
de fase negativa. Esses diagnósticos se referem às amostras aceitas, não a
toda proposta Newton.

A versão do índice rejeita snapshots em andamento e identidades divergentes,
separa término de acomodação e marca comparações provisórias; sete testes de
corrupção/contrato passaram. A preparação da extensão foi também exercitada
contra casos pendentes e recusou criar a nova pasta. Só depois dos seis
auditores terminais foram preparados Ps+ e jg+ para 7200 s, em coorte distinta,
com mudança exclusiva de time_dim/out_dir. Esse adendo terá evidência própria.

A comparação experimental continua limitada pela instrumentação/ponderação
de alpha_p, pela diferença entre EOS e conversão STANDARD, pela dinâmica
dos contornos experimentais e pela geometria/configuração histórica.
Nenhum sinal temporal foi inventado nem parâmetro ajustado para produzir ciclos.
