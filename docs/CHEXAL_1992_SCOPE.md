# Chexal 1992 e Ozaki 2015: escopo para continuidade da auditoria

**Estado em 18/11/2026: fontes recebidas e escopo documental verificado;
Chexal não implementado.** O solver da campanha continua com o fechamento
Bhagwat circular da transcrição de Azevedo, Baliño e Burr (2015). Receber uma
correlação que trata algumas situações de contracorrente não demonstra que ela
resolva as falhas atuais de junção, discretização ou seleção de raízes.

## 1. Qual fonte é `CHEXAL_1992`

B. Chexal, G. Lellouche, J. Horowitz e J. Healzer, **A void fraction correlation
for generalized applications**, *Progress in Nuclear Energy*, volume 27,
número 4, páginas 255–295, 1992. O cabeçalho bibliográfico indica 1992; o
copyright impresso indica 1993. Esta identificação não deve ser misturada
silenciosamente com variantes posteriores de 1997 ou com um pacote de
implementação para outra geometria.

PDF fornecido na sessão: `project_sources/06-chexal1992.pdf`.
SHA-256: `b0b6b68f5e13fec48e931d89d87f9d65cb9c3bb134235c2307ded24e9c00700e`.
A página PDF é a página impressa menos 254. O PDF de terceiros não acompanha
o repositório nem o ZIP de evidências.

## 2. Domínio de orientação e sinais

| Situação | O que a fonte de 1992 estabelece | Localização |
|---|---|---|
| Ângulo | Graus medidos a partir da vertical, entre 0° e 90°; 0° é vertical e 90° é horizontal | p.262 / PDF8, Fig.3-1 e texto |
| Vertical cocorrente | Fluxos superficiais positivos para ambas as fases ascendentes e negativos para ambas descendentes | p.262 / PDF8 |
| Vertical contracorrente | Gás positivo, ascendente; líquido negativo, descendente | p.262 / PDF8 |
| Contracorrente inclinada | Não é qualificada por essa extensão: o texto restringe a contracorrente considerada ao caso vertical | p.262 / PDF8 |
| Horizontal | Apenas cocorrente; a prescrição usa fluxos volumétricos positivos e exclui fluxos negativos | pp.264–265 e 269 / PDF10–11 e 15 |
| Reynolds | Preservar as definições e sinais empregados nos ramos; aplicar módulo somente onde a fórmula o pede | p.264 / PDF10 |

O PRIMS-F usa ângulo em radianos a partir da horizontal e coordenada positiva
do inlet ao separador. No riser ascendente, a transformação geométrica seria
`theta_Chexal_deg = 90 - theta_PRIMSF*180/pi`; ela precisa ser acompanhada
pela convenção de sinais dos fluxos. Isso **não** autoriza usar o ramo de
contracorrente vertical em toda a catenária. A catenária D começa horizontal
e chega ao topo com cerca de 76,47° a partir da horizontal: nem mesmo o topo
é uma seção exatamente vertical. Uma extensão de contracorrente inclinada
exigiria outra fonte ou uma hipótese de modelo explicitamente identificada.

O fator de fluido `L` também depende da família física. Para ar e água,
conferir as equações 3-7, 3-8 e 3-19; não copiar automaticamente o fator de
vapor e água que depende da pressão crítica. A descrição do artigo e suas
expressões devem ser mantidas como conjunto coerente, incluindo C0 e Vgj.

## 3. Duas raízes e limite de contracorrente

A página 268 / PDF14, equações 3-39 a 3-41, distingue a linha de limitação de
contracorrente, **CCFL**, da região abaixo dela. O artigo afirma que nessa
região existem duas soluções para a fração de vazio. A escolha da solução
desejada precisa de informação anterior, como queda de pressão ou outra
informação física; a fórmula de C3 depende do ramo escolhido.

Consequências para uma implementação futura:

- Não escolher sempre a menor ou a maior raiz por conveniência numérica.
- Informar número de raízes, ramo, domínio e falhas de solução ao chamador.
- Se houver continuação por um estado anterior, declarar essa regra e
  verificar sua compatibilidade com a informação física exigida pela fonte.
- Tratar a determinação do fluxo na linha CCFL e a expressão correspondente
  de C3 como parte do modelo; não substituir isso por um simples `abs(jl)`.
- Testar limites de fase e de fluxo sem confundir fluxo nulo com ausência de
  inventário. Uma fase pode estar presente e estacionária.

A seção 4.0, p.270 / PDF16, registra que a base de **vapor e água** cobre
cocorrente ascendente, descendente e condições sobre a linha CCFL, mas que
não foram localizados dados fora da linha CCFL para verificar a modelagem
da região de contracorrente. Essa ressalva pertence à base descrita nessa
seção; não deve ser generalizada como afirmação de inexistência de qualquer
experimento de ar e água em contracorrente.

## 4. Cuidados de transcrição

A equação 3-1, p.261 / PDF7, apresenta uma aparente repetição de `jg` na soma
dos fluxos. A auditoria deve confrontá-la com a definição de fluxo total,
com a identidade drift-flux e com as expressões do apêndice A, pp.289–291 /
PDF35–37, equações A12–A16. Qualquer correção tipográfica deve ser anotada,
com fonte e justificativa, em vez de ser introduzida silenciosamente.

Para o novo módulo, congelar antes uma especificação rastreável de cada
constante, unidade, ângulo, ramo, critério de transição e limite. A
derivação ou transcrição de uma correlação não equivale à validação do
solver transiente que a utiliza.

## 5. O que Ozaki e Hibiki (2015) acrescentam

T. Ozaki e T. Hibiki, **Drift-flux model for rod bundle geometry**, *Progress
in Nuclear Energy*, volume 83, páginas 229–247, 2015.
DOI: [10.1016/j.pnucene.2015.03.015](https://doi.org/10.1016/j.pnucene.2015.03.015).

PDF fornecido: `project_sources/07-ozaki2015.pdf`.
SHA-256: `103b31e1ad5df5c3499f577d569bcefe6432620ddfb55acdd26282d4d7a96cf3`.
A página PDF é a página impressa menos 228.

Na seção 2.2.2, p.230 / PDF2, os autores descrevem a existência de vários
pacotes Chexal–Lellouche conforme condição térmica, orientação e direção de
escoamento, incluindo contracorrente. A Tabela1, p.231 / PDF3, resume
correlações para **feixes de varetas**. A avaliação experimental do artigo
usa as bases FRIGG e NUPEC; ela não é uma validação de uma catenária de ar e
água e não remove a restrição de orientação da fonte Chexal de 1992.

A seção 4.3, pp.241–243 / PDF13–15, chama atenção para a compensação de erros
entre C0 e a velocidade de drift: combinar C0 de um modelo com Vgj de outro,
ou validar cada um isoladamente apenas por coincidência de alfa, pode
produzir uma conclusão enganosa. Para a comparação futura com Bhagwat,
manter cada par C0/Vgj com sua própria definição, domínio e fonte.

## 6. Relação com as falhas observadas no PRIMS-F

Os diagnósticos atuais distinguem: fluxo gasoso negativo em uma **trace
numérica** de junção, aproximação de `ul=0` em uma célula, estado intermediário
de Roe fora do domínio e raiz positiva fora do intervalo pesquisado. Nenhum
desses sinais isolados comprova reversão física no experimento.

O candidato de cauda de alfa permanece separado, não aplicado à produção.
Ele corrige uma limitação de busca para ID3 cocorrente; não acrescenta
contracorrente a Bhagwat. A hipótese de um novo `CHEXAL_1992` deve ser
auditada separadamente e só comparada depois de verificados os estados
admissíveis, as raízes, as derivadas e a consistência da junção.

Para resultados físicos, manter estudos independentes de malha e passo,
pressão média, amplitude, período, produção e inventários. Um evento de
bloqueio forte seguido de retomada do gás é um detector específico útil
para SS1; sua ausência não exclui OSC ou SS3, que podem manter penetração
gasosa. Conservação de massa, conclusão de uma execução e validação
experimental são verificações diferentes.

Os alvos primários agora conferidos visualmente são documentados em
`BENCHMARK_REFERENCE.md` e `benchmark_sources.json`: SS3 reportado de 76s
na Tabela3, caso11, de Baliño et al. (2010); OSC com faixa de19,3–33,8s no
caso13. Os 26s de SS3 e os33,9s de OSC impressos na dissertação continuam
registrados separadamente, sem alteração dos parâmetros D executados.
