# Comparação numérica e sequência de sinais físicos

Os scripts desta nota usam apenas a biblioteca padrão do Python e leem
`monitor.csv` por nomes de coluna. Não alteram o solver, os parâmetros físicos,
os resultados nem o analisador principal `analyze_validation.py`.

## Execução

Depois de atualizar o `metrics.json` do analisador principal:

```sh
python scripts/compare_validation.py --metrics validation/analysis_final/metrics.json --output validation/analysis_final
python scripts/analyze_phase_sequence.py --metrics validation/analysis_final/metrics.json --output validation/analysis_final
```

Para reunir as campanhas com prefixos `final_` e `refined_`, acrescente
`--pattern '[fr]*'` a ambos. `--results` escolhe outro diretório de resultados;
`--run-dir` pode ser repetido para acrescentar diretórios individuais. O arquivo
de métricas é opcional: sem ele, os monitores fornecem estatísticas e um detector
de período por cruzamentos com histerese, identificado separadamente no JSON.

A janela padrão começa em **900 s**. Um run que termina antes disso aparece
como janela insuficiente: não recebe automaticamente uma janela curta de piloto.
Para inspecionar explicitamente o transiente completo, use `--window-start 0`
e, de preferência, outro diretório de saída para preservar o relatório principal.
`--window-end` limita o fim da observação.

Saídas:

| Script | JSON | Markdown |
|---|---|---|
| `compare_validation.py` | `comparison.json` | `comparison.md` |
| `analyze_phase_sequence.py` | `phase_sequence.json` | `phase_sequence.md` |

## Estatísticas no tempo

Os tempos são os tempos físicos efetivamente observados. As médias integram a
interpolação linear por trechos, e o desvio padrão usa o segundo momento da
mesma função, com centralização para evitar cancelamento ao subtrair pressões
grandes. Assim, regiões com mais amostras devido ao passo adaptativo não recebem
peso temporal indevido.

Dados faltantes são preservados como faltantes e reduzem a cobertura. O programa
não os converte em zeros. O JSON identifica a coluna usada para pressão,
estatísticas de vazão de cada fase, inventários totais e por trecho, frações
volumétricas e os passos registrados no monitor.

Quando disponíveis e compatíveis com toda a janela solicitada, os picos do
analisador principal fornecem os períodos. Um arquivo de métricas desatualizado
não é usado para representar o trecho adicional de um run já concluído. Na
ausência de picos apropriados, o detector alternativo observa uma passagem
ascendente por 80% da faixa de pressão depois de ter observado um valor abaixo
de 20%. Ele registra somente candidatos de período, exige excursão mínima de
100 Pa e desvio de pelo menos 1 Pa, e não cria um período a partir de um sinal
plano ou de uma subida monotônica. Os métodos são identificados, e diferenças
de período de métodos diferentes não são apresentadas como comparação direta.

## Comparação de refinamentos

O programa monta pares de dois tipos:

1. Mudança de número de células mantendo `dt_max_s` e as demais escolhas iguais.
2. Mudança de `dt_max_s` mantendo as duas contagens de células e as demais
   escolhas iguais.

Geometria, fluidos, EOS, vazões, pressão, temperatura, fechamento, junção,
dissipação, fluxo numérico, perturbação, reconstrução e controladores numéricos
fazem parte da assinatura da família comparada. Alterar malha e teto temporal
simultaneamente não produz um par de refinamento independente. O hash do
documento de referências não faz parte dessa assinatura: corrigir uma ligação
documental não muda a física de uma execução.

Hashes diferentes das fontes do solver e dos executáveis são sinalizados.
As métricas de cada par são recalculadas sobre a interseção de suas janelas.
Execuções incompletas podem fornecer comparações descritivas de seu trecho
comum, mas permanecem marcadas como incompletas. Pressões medidas em locais
diferentes não recebem uma diferença direta.

`dt_max_s` é o **teto** do passo adaptativo. O relatório também mostra os passos
observados; a simples mudança do teto não garante que todos os passos tenham
sido reduzidos proporcionalmente. Os monitores podem não registrar cada passo
aceito.

Os avisos dinâmicos usam os seguintes limiares de triagem declarados:

| Comparação | Aviso |
|---|---|
| Oscilação resolvida presente somente em uma resolução | Mudança de presença de oscilação sob refinamento |
| Desvio ou amplitude de pressão diferem mais de 10%, acima dos respectivos pisos | Dinâmica não estável entre os controles comparados |
| Períodos diferem mais de 5% | Diferença material de período |
| Desvios nos três terços da janela têm dispersão relativa acima de 20% | Amplitude variável na janela; não assumir ciclo limite estacionário |

Esses limiares são escolhas do pós-processamento, não tolerâncias experimentais
nem ajuste de parâmetros. O programa não calcula ordem de convergência ou GCI
e nunca aprova convergência automaticamente. Um balanço de massa correto não
substitui concordância de amplitude, pressão, inventários, produção e período.

## Sequência física: observações e proxies

O detector preserva a distinção entre as séries calculadas pelo solver e sua
interpretação física. Para construir a sequência completa, exige pressão,
`Ml_riser_kg`, `Mg_pipe_kg`, produção de líquido e `mdot_g_base_kg_s`. Não
substitui os inventários por trecho pelo inventário total quando um canal falta.

| Indicador | Procedimento padrão | Limite da interpretação |
|---|---|---|
| Acumulação/compressão | Aumento simultâneo de líquido no riser, gás no pipeline e pressão, avaliado em uma janela de 2 s | Aumento de pressão é uma evidência indireta de compressão; não mede frente ou volume de bolha |
| Produção de líquido | Excursão acima da entrada efetiva, com histerese de amplitude e duração mínima | Mede produção; não identifica, sozinho, a passagem de um slug |
| Retomada de gás | Fluxo na base abaixo de 5% da entrada por pelo menos 0,5 s contínuo, seguido de subida acima de 25% | Proxy de bloqueio forte/liberação; fluxo baixo não é automaticamente ausência de gás |
| Despressurização/blowdown | Queda de pressão com transferência positiva de gás pela base | Proxy; inventário de gás e saída de gás acrescentam contexto |

A produção volumétrica é `A*jl_out` quando esse canal e o diâmetro estão
disponíveis. A entrada volumétrica de referência prioriza o fluxo mássico
efetivamente imposto dividido por `rho_l0`; o valor nominal `jl` não sobrepõe um
override de `mul`. Na ausência de volume calculável, a produção mássica é
identificada como proxy, com unidades explícitas.

Os limiares absolutos e relativos completos aparecem no JSON. O encadeamento
permite sobreposição de até 2 s entre etapas, por causa da resolução temporal do
detector, e limita a extensão total de um candidato a 300 s. Esses parâmetros
podem ser configurados e não são períodos experimentais prescritos.

### A sequência com bloqueio não é um critério universal

A cadeia que inclui retomada após fluxo baixo é útil para uma narrativa de
**bloqueio forte e liberação, do tipo SS1**. **OSC e SS3 podem manter passagem
gasosa:** a ausência dessa cadeia não exclui essas instabilidades nem todo
comportamento chamado slugging. Por isso, oscilações de pressão acompanhadas
de fluxo de gás positivo nos tempos observados recebem um resultado descritivo
próprio, sem exigir uma etapa de bloqueio.

Uma única sequência é um transiente candidato. Mesmo várias sequências não
constituem validação experimental de frente, regime ou amplitude. O código
retém explicitamente `automatic_slugging_validation=false` e
`automatic_flow_regime_classification=false`.

## Verificação dos scripts

```sh
python tests/test_validation_evidence.py
```

Os testes incluem estatística exata de uma função linear amostrada
irregularmente, canais faltantes, sinais planos, pressão monotônica, pares
independentes e parâmetros incompatíveis, runs incompletos, mudanças de código,
mudança de amplitude apesar de conservação, sequência sintética ordenada,
oscilações com gás continuamente positivo, prioridade da massa imposta e
rejeição de uma falsa duração obtida somando visitas separadas a fluxo baixo.
Os sinais sintéticos verificam o software; não são benchmarks físicos.
