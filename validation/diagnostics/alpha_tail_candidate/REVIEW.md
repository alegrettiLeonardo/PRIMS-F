# Revisão independente da correção isolada da cauda de alfa

**Conclusão:** não encontrei defeito bloqueante na proposta mínima. A correção
resolve o falso `STEADY_NO_ROOT` demonstrado para raízes positivas abaixo da
malha original, sem estender a correlação ao escoamento contracorrente nem
alterar os sucessos anteriores. Esta é uma revisão local do helper, não uma
aprovação dos resultados transientes ainda não repetidos.

Fonte revisado: versão produzida por `alpha_tail.patch`, SHA-256
`76311d3a412b6699ccb25a93c74062b2f993181bbecd9aff0928530a02c01624`.
Produção preservada: SHA-256
`73fbc6f2359b450d2d582439e65fd7942084a4d6c68a3c1426c4a6b668d40450`.
Não alterei o fonte de produção nem a proposta.

## Verificação do algoritmo

A chamada adicional ocorre somente depois da varredura antiga e somente com
`ns==0`, ID3, `theta>=0`, `jl>=0` e `jg>0`. O domínio restante é verificado pelo
próprio fechamento. Declarações monofásicas, escolhas entre múltiplas raízes,
entradas inválidas e raízes já encontradas conservam suas operações anteriores.

Com fluxos superficiais fixos, Reynolds, fração de fluxo e qualidade são
constantes. Na faixa angular admitida, a forma implementada é
`C0(a)=A+w*B**((1-a)/5)` e `Ud(a)=K*sqrt(1-a)`, com
`1/2<=B<=1`, `0<=w<=1`, `K>=0`. Portanto
`0<=C0'(a)<=ln(2)/5`. A margem verificada no limite superior da cauda garante
`C0>0` em toda a cauda. A derivada do residual é estritamente positiva para
`a<2/3`, pois `d[a*sqrt(1-a)]/da>0` nessa faixa. Essa justificativa sustenta a
unicidade **dentro da cauda adicional**.

A busca exige residual positivo no limite antigo e desce geometricamente até
um residual negativo ou zero representável. A bisseção usa largura relativa
ao valor de alfa, ou encerra quando não existe ponto médio representável,
escolhendo o extremo válido com menor residual absoluto. O mínimo de
`xtol=8*epsilon` e o máximo de `digits+8` iterações são compatíveis com esse
intervalo inicial de razão dois. O teste de velocidade impede a divisão
adicional por uma fração que causaria overflow adimensional. Falhas do
fechamento ou falta de intervalo representável não produzem uma raiz
substituta. O número máximo de descidas cobre a faixa de expoentes binários
incluindo subnormais.

Como `record` é chamado somente quando `ns` era zero, a cauda acrescenta no
máximo um registro, mantém `nroots=1` e não envolve desempate por `alpha_hint`.
Os dois extremos admissíveis de hint, zero e um, retornaram a mesma raiz. O
contador continua significando raízes **detectadas**; a prova local não deve
ser apresentada como prova de inexistência de raízes tangentes ou na cauda
superior não pesquisada.

## Teste independente

O runner da revisão histórica compila a proposta e a cópia anterior separadamente com GNU
Fortran 14, `-std=f2018 -O2 -fcheck=all` e traps de invalid, zero e overflow.
As entradas combinam as sete referências independentes de 80 dígitos com
seis escalas de densidade/velocidade e tolerâncias `1e-6`, `1e-12` e `1e-14`.
Também verificam as exclusões de fluxo gasoso nulo e negativo.

| Contrato verificado | Resultado |
|---|---:|
| Raízes novas abaixo do limite antigo | 108/108 |
| Sucessos antigos preservados byte a byte | 18/18 |
| Exclusões preservadas byte a byte | 3/3 |
| Erro relativo máximo nas raízes novas | 3,26295e-13 |
| Residual relativo máximo nas raízes novas | 3,26588e-13 |
| Hash do fonte de produção preservado | Sim |

Os resultados completos, as opções de compilação e os hashes dos resumos
históricos estão em `prior_evidence.json`. Os logs brutos, as entradas e o
runner original pertencem à revisão histórica; não são arquivos disponíveis
neste pacote. A reprodução portátil usa `python3 run_tests.py --review`.
As raízes novas testadas
chegam a aproximadamente `1e-300`; os testes próprios da proposta também
cobrem `1e-320`, o menor positivo representável e a recusa de uma raiz abaixo
da representabilidade. Conferi os logs e hashes dessa última suíte: 55 guards
e 398 verificações estacionárias, O0/O3, com 101 registros antigos idênticos.

Reprodução portátil: consulte `README.md` e execute `python3 run_tests.py --review`
a partir deste diretório. Os resultados históricos completos estão em
`prior_evidence.json`; esta nota descreve a revisão já realizada, não uma nova execução.

## Limitação anterior identificada e preservada

Com `jg=6,445998075334944e-13 m/s`, raiz exata aproximada
`7,22357939817937e-13`, e **tolerância frouxa `1e-6`**, o Brent antigo aceita
o extremo `6,914400106935423e-13`, com erro relativo de 4,2801397%. Esse
comportamento ocorreu nas seis escalas e foi **byte a byte igual** na versão
anterior e na proposta: o critério absoluto antigo de intervalo é amplo
comparado a essa raiz minúscula. A tolerância de produção `1e-12` não
apresentou esse erro nos casos exercitados. A proposta deliberadamente não
muda o caminho antigo. Esse fato aparece separado, como limitação anterior,
no relatório JSON; não foi ocultado como aprovação de precisão.

## Correção de evidência durante a revisão

O runner inicial sobrescrevia o log do modo `emit` ao executar `check`, embora
comparasse corretamente os textos em memória. Comuniquei o problema ao autor,
que separou os diretórios por modo e repetiu a suíte. Na revisão histórica,
conferi os dez logs e seus hashes. O resumo original chamado `test_summary.json`
não é um arquivo deste pacote: sua síntese, os resultados e o hash do resumo
original estão preservados em `prior_evidence.json`.

Nenhuma execução transiente foi realizada por esta revisão. A aplicação
posterior e a comparação dos casos interrompidos devem manter versões e
resultados anteriores identificados, sem atribuir ao patch a resolução dos
limites distintos de `ul=0`, junção com fluxo gasoso negativo ou trajetória de
Roe fora do domínio.
