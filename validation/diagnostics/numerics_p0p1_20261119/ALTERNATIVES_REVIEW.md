# Revisão suplementar: SS1, SS3 e alternativas de continuação

Consulta de 2026-11-19, sem alteração do solver. Complementa o README desta
pasta e a prova em `docs/SS3_COCURRENT_BOUND_20261119.md`, commit `f1f9ff0`.

## Alcance da evidência

O diagnóstico do README continua referido à base `9f2ed21`: SS1 falha no
traço de junção; SS3 falha na recuperação da célula 69. A prova posterior
acrescenta um resultado mais forte para SS3: exclui uma raiz em **todo o
intervalo cocorrente**, para cada uma das 18 propostas exatas, nos dois modos
de Jacobiano. A linearização ainda tem raiz, com líquido positivo minúsculo.
Os valores positivos do limite inferior foram confirmados a 80 dígitos.

Esse certificado não exclui todo estado conservativo possível, todas as
direções de Newton nem toda solução do próximo passo temporal. Também não
constitui uma solução contracorrente. A distinção entre busca mal resolvida,
raiz excluída, domínio constitutivo e falha de globalização deve permanecer.

## Alternativas e hipóteses necessárias

| Alternativa | O que a fonte/evidência oferece | Limite para estes casos |
|---|---|---|
| Frente de Baliño2010 | Eqs. 1–14, pp. 645–646, definem acumulação, compressão e troca reversível com gás bloqueado; Eqs. 59–67, p. 648, discretizam os estados. | É um modelo com frente e inventário próprios, pipeline simplificado e riser NPW. Importar apenas `jg=0` não conserva esse conjunto de hipóteses. |
| Retorno líquido no modelo de Baliño/Bendiksen | A Eq. 26/p. 647 permite líquido negativo conforme `alpha`, fluxo total e drift; Eq. 68/p. 648 acompanha o nível líquido, e p. 650 discute fallback. | Isso não valida toda contracorrente inclinada para Bhagwat, nem transforma o problema inercial atual no NPW. Trocar o fechamento altera a física e requer nova comparação com propriedades e entradas explícitas. |
| Chexal1992 contracorrente | A p. 262 admite gás ascendente/líquido descendente **somente na vertical**; p. 268 admite duas raízes, CCFL e seleção por informação adicional. | A célula 69 está a aproximadamente 51,64° desde a horizontal, não na vertical. A p. 264 exclui contracorrente horizontal; portanto também não resolve o traço SS1 na base. |
| Riemann não linear com fase ausente | Pode ser formulado para conservar massas e selecionar estados/ondas admissíveis, se o limite constitutivo e as condições da junção forem derivados. | Conservação isolada não prova positividade, admissibilidade de ondas, entropia ou consistência do acoplamento. O matching linear testado já produz massa gasosa direita negativa em SS1. |
| Contagem de características P0 | Nos snapshots, duas condições entram pela alimentação e uma pelo separador; as restrições projetadas são independentes. | Esse resultado local não amplia o domínio constitutivo, não cria uma frente e não garante estabilidade global da junção. |

As páginas de Chexal foram relidas diretamente em `chexal1992(1).pdf` nesta
revisão. Mantêm-se Reynolds assinados, ângulo desde a vertical e a versão de
1992. Não foi encontrada aqui uma fonte primária que autorize sua extensão
contracorrente à inclinação da célula 69.

Para pesquisa numérica adicional, o resumo oficial de
[Baudin, Coquel e Tran, *A Semi-implicit Relaxation Scheme for Modeling
Two-Phase Flow in a Pipeline*](https://doi.org/10.1137/030601624) descreve um
esquema explícito para os autovalores lentos e implicitamente linearizado para
os rápidos, em um modelo com leis algébricas de pressão e hidrodinâmica.
Isso é uma direção de pesquisa para a separação de escalas observada; o resumo
consultado não demonstra que tal método resolva esta junção com Bhagwat ou
seu endpoint de fase ausente. A transferência exige a leitura e verificação
das hipóteses do artigo completo. Nenhum coeficiente desse esquema foi usado.

## Busca por raio, trust region e matching linear

A busca ao longo de uma direção modifica apenas o comprimento do passo.
Testar mais lambdas pode encontrar uma proposta admissível ainda não testada;
isso não autoriza assumir que o residual completo diminuirá ou que a solução
atravessará uma fronteira física. Nos 18 lambdas efetivamente ensaiados há
recusas verificadas; a conclusão não se estende automaticamente a um contínuo
de lambdas não examinado.

Trust region pode mudar a direção e contornar uma dificuldade algorítmica
dentro do domínio existente. É uma alternativa numérica legítima a estudar,
mas sucesso exige resolver **as mesmas equações**, com residual e balanços
convergidos, sem aceitar uma raiz incompatível com o ramo. Uma minimização
restrita pode terminar no limite com residual não nulo; isso não é uma solução
transiente. Introduzir uma condição complementar de bloqueio sem derivá-la de
um modelo físico acrescentaria uma equação, não apenas globalização.

O sistema característico linear ensaiado tem posto completo, mas sua solução
SS1 viola a massa gasosa positiva. Curvas de onda não lineares podem produzir
outro resultado; essa possibilidade permanece aberta, sem garantia. No
endpoint de fase ausente, velocidades da fase inexistente e a dimensão do
problema precisam de tratamento derivado. Não se deve extrapolar um espectro
ou substituir raízes para manter artificialmente três ondas admissíveis.

O caminho técnico seguinte deve separar dois trabalhos: primeiro, testar uma
globalização alternativa dentro do domínio atual, conservando as equações e
documentando eventual estagnação; segundo, derivar e verificar um operador
de interface/limite de fase, ou implementar um modelo de frente selecionável
com as equações primárias completas. Cada alternativa precisa de benchmarks
de ondas/limites, conservação por fase, pressão de interface, reversibilidade
do evento e comparação independente. Oscilar pressão não substitui esses
testes. Nenhuma dessas alternativas foi implementada nesta revisão.

## Conferência preparada para os novos SS1/SS3

Na leitura inicial, o plano `v1_2400/plan.json` estava `PREPARED_UNSEALED` e
não havia manifests de execução. Quando existirem, a conferência será feita
apenas sobre os artefatos novos, preservando os anteriores:

1. Conferir seal, SHA do binário/fontes, deck e metadados efetivos; separar
   explicitamente mudança de tolerância de qualquer alteração física.
2. Conferir status de processo, tempo alcançado, passos aceitos/rejeitados e
   término real. Não deduzir causa terminal de um código de proposta anterior.
3. Relacionar sidecar, estado aceito, linearização e proposta realmente
   recusada por tempo, iteração, lambda e hash; norma não montada permanece
   indisponível, nunca igual ao sentinela `huge` como valor físico.
4. Localizar o primeiro operador recusado e seu status interno. Em SS1,
   verificar sinais/pressão/fluxos dos traços; em SS3, certificar ou não o
   intervalo cocorrente usando o estado novo, sem reutilizar a prova de outro U.
5. Conferir os balanços aceitos e defeitos de Newton globais e por trecho,
   incluindo cancelamento do fluxo comum de base. Uma mudança de código de
   erro ou captura melhor não significa avanço físico do transiente.

Referência física principal: [Baliño, Burr e Nemoto2010](https://doi.org/10.1016/j.ijmultiphaseflow.2010.04.003).
Os PDFs seguem fora do Git; nenhuma reprodução de página foi adicionada.

## Compatibilidade de métodos numéricos

- O capítulo 7 de LeVeque trata ondas entrantes e células fantasmas
  em problemas hiperbólicos. A contagem 2/1 foi calculada para o
  Jacobiano do PRIMS-F e conferida nos estados efetivamente impostos.
  Uma contagem correta é uma condição local de compatibilidade. A
  construção não linear do traço e as condições em reversão precisam
  de verificação própria. [LeVeque, 2002](https://doi.org/10.1017/CBO9780511791253).
- A região de confiança interior de Coleman e Li, pp. 418–445, trata
  minimização com limites nas variáveis. Aqui a restrição envolve uma
  raiz constitutiva e o fluxo reconstruído em uma face. Uma aplicação
  precisaria incluir essas funções de admissibilidade e exigir resíduo
  conservativo convergido. Os replays dos dois snapshots são testes
  iniciais para qualquer formulação desse tipo.
  [Coleman e Li, 1996](https://doi.org/10.1137/0806023).
- O Teorema 4.1 de Guermond e Popov assume conjunto convexo invariante,
  esquema explícito e condição CFL. Proteger a massa média da célula
  por um método conservativo desse tipo não certifica o traço entre
  duas leis nem as equações de fase ausente no FV implícito atual. Uma
  extensão exigiria teste do domínio dos dois traços, cancelamento dos
  fluxos por fase, limite monofásico e consistência do Jacobiano.
  [Guermond e Popov, 2016](https://doi.org/10.1137/16M1074291).
