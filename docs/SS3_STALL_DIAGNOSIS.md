# Diagnóstico das interrupções do piloto D_SS3

## Resultado e escopo

Foram identificados **dois mecanismos distintos**, em estados congelados do
piloto D_SS3 com 24 volumes no pipeline e 24 no riser. O primeiro surge no fluxo
numérico de Roe ao atravessar uma mudança da derivada do atrito estratificado.
O segundo está na própria relação algébrica de velocidades do fechamento
Bendiksen por partes utilizado como alternativa explícita à correlação de
Bhagwat citada na dissertação.

| Estado | Tempo físico | Local relevante | Evidência | Consequência |
|---|---:|---|---|---|
| SS3 inicial | 50,002215052 s | Pipeline, índices estendidos 23/24 | Fluxo físico contínuo; fluxo de Roe instantâneo com salto finito ao cruzar `Re_g=2000` | A derivada local de fechamento correta, isoladamente, não torna contínua a dissipação instantânea |
| `pilot5_D_SS3_MATCHED` | 78,6073065235 s | Riser, índice estendido 49 | `Fr=3,49999999409`; lacuna de momento de 30,8996984 kg m⁻² s⁻¹ sem raiz algébrica | Nenhuma precisão adicional do solver encontra uma raiz que a relação por partes não admite |

Esses resultados diagnosticam o modelo e a discretização implementados.
**Não demonstram captura experimental de SS3, período, amplitude ou slugging.**
Ambos os transientes terminaram antes de uma janela útil de ciclos comparáveis.
A regularização de uma correlação, a seleção de outro fechamento e a escolha
de dissipação são mudanças distintas e devem ter nomes e validações próprias.

## Dados congelados e reprodução

Os índices incluem duas células fantasmas à esquerda: as células físicas são
`3:50`, o pipeline é `3:26`, e o riser é `27:50`. O índice 49 corresponde ao
penúltimo volume físico do riser. Os tempos nos arquivos `snapshot.dat` são
adimensionais; neste caso, `omega_t=Lr/0.48`.

Todos os caminhos abaixo são relativos à raiz do pacote.

| Arquivo | Função |
|---|---|
| `validation/diagnostics/ss3_stall_snapshot.dat` | Estado congelado a 50,002215052 s; não é o arquivo de rejeição vivo do piloto |
| `validation/diagnostics/ss3_stall_metadata.txt` | Parâmetros do primeiro diagnóstico |
| `validation/diagnostics/ss3_stall_regimes.csv` | Reynolds, Froude e primitivas por célula |
| `validation/diagnostics/probe_ss3_stall.f90` | Comparação do resíduo/Jacobiano com tolerâncias e passos de diferenças finitas distintos |
| `validation/diagnostics/probe_ss3_branch.f90` | Aproximação dos dois lados de `Re_g=2000`, preservando as massas locais e a célula vizinha |
| `validation/diagnostics/ss3_branch_baseline.csv` | Fechamento com derivada por diferenças finitas anterior |
| `validation/diagnostics/ss3_branch_analytic.csv` | Mesmo fechamento primal; derivada estratificada analítica |
| `validation/diagnostics/ss3_late_snapshot.dat` | Estado congelado do segundo diagnóstico, em 78,6073065235 s |
| `validation/diagnostics/ss3_late_metadata.txt` | Configuração do segundo piloto; dissipação `LAGGED`, fechamento do riser 2 |
| `validation/diagnostics/ss3_late_field_summary.json` | Último perfil aceito completo, com Froude calculado |
| `validation/diagnostics/ss3_late_bendiksen_gap.json` | Cálculo independente dos dois ramos e seus limites de admissibilidade |
| `validation/diagnostics/probe_ss3_bendiksen_gap.f90` | Verificação executável da recuperação conservativas → primitivas dentro e fora da lacuna |
| `validation/diagnostics/ss3_bendiksen_gap_probe.csv` | Estados amostrados e códigos de retorno dessa verificação |

Os scripts `build_probe.py`, `build_branch.py` e `build_gap.py` compilam os
probes e registram comandos e hashes SHA-256 das fontes efetivamente usadas nos
respectivos diretórios de compilação. O probe de ramo não depende do módulo de
integração temporal. A cópia `baseline_sources` preserva as fontes constitutivas
anteriores à derivada analítica; não deve ser descrita como uma imagem binária
integral do piloto, pois havia alterações paralelas de diagnóstico no módulo de
simulação.

Exemplos, a partir da raiz do pacote, com um compilador Fortran 2018:

```sh
FC=gfortran python validation/diagnostics/build_branch.py baseline
FC=gfortran python validation/diagnostics/build_branch.py analytic
FC=gfortran python validation/diagnostics/build_gap.py
```

O compilador utilizado nesta sessão foi GNU Fortran 14.2. Os scripts dos probes
antigos de Roe ainda contêm `-ftrampoline-impl=heap`; essa opção não é uma
necessidade física nem uma exigência do probe algébrico da lacuna. O executável
do pacote usa sua configuração de compilação própria.

## 1. A interrupção próxima de 50 s

### 1.1 Estado local

As células adjacentes relevantes tinham:

| Índice | Fração de gás | `Re_l` | `Re_g` | `Fr` da mistura |
|---|---:|---:|---:|---:|
| 23 | 0,7237030521 | 43036,2075 | **2000,0269066** | 0,78534356 |
| 24 | 0,7237031959 | 43036,0777 | **1999,9874026** | 0,78507345 |

O gás de uma célula estava imediatamente acima do início da transição de
atrito, e o da outra imediatamente abaixo. No perfil completo, os valores de
Froude eram aproximadamente 0,61–0,93; o limiar de Bendiksen `Fr=3.5` não
explica esse primeiro estado.

A regra de Fanning implementada é contínua no início de sua transição, mas sua
derivada muda de ramo. Portanto, o fechamento primal estratificado pode ser
contínuo sem ter a mesma derivada dos dois lados. Essa diferença propaga-se ao
Jacobiano físico recuperado e às velocidades características.

### 1.2 Tolerância da inversão não explica o piso do resíduo

O estado congelado registrava resíduo normalizado `1.61399e-8`. A recomposição
durante o probe deu aproximadamente `1.39491473e-8`; apertar simultaneamente as
tolerâncias de recuperação para `tol=1e-13` e `tola=1e-14` deu
`1.39488653e-8`. As diferenças absolutas nos componentes do resíduo foram da
ordem de `1e-14`.

Isso elimina, para esse estado, a hipótese de que o piso fosse dominado pela
precisão insuficiente da raiz de velocidade. Não é uma afirmação universal
sobre todas as tolerâncias de todos os estados do solver.

Diferenças finitas centradas com passos `1e-5`, `1e-6` e `1e-7` atravessavam
ramos constitutivos em algumas dessas células. Nessa situação, a diferença
finita centrada ampla mede uma secante entre ramos e não é um oráculo para a
derivada unilateral de um ramo. A comparação acusou diferenças da ordem de
6,5% em componentes de Jacobiano nas células 23/24.

### 1.3 Experimento que distingue fluxo físico e dissipação numérica

Fixaram-se as massas conservativas da célula 23 e o estado da célula 24. Apenas
o momento da célula 23 foi aproximado pelos dois lados do valor crítico

\[
U_{3,\mathrm{crit}}=0.58479072058637738,
\]

que produz `Re_g=2000`. O fluxo numérico examinado foi

\[
\widehat F=\frac{F_L+F_R}{2}
-\frac12\left|\frac{A_L+A_R}{2}\right|(U_R-U_L).
\]

Com a derivada estratificada **analítica**, mantendo o fechamento primal
inalterado, obteve-se:

| Perturbação `delta_U3` em cada lado | Diferença do fluxo físico de gás | Diferença do fluxo numérico de gás |
|---:|---:|---:|
| `1e-6` | `1.11689648e-6` | `1.30263402e-5` |
| `1e-7` | `1.11695003e-7` | `1.25229124e-5` |
| `1e-9` | `1.11695589e-9` | `1.24675327e-5` |
| `1e-11` | `1.11695152e-11` | `1.24669789e-5` |
| `1e-13` | `1.11654826e-13` | **`1.24669734e-5`** |

São fluxos nas unidades adimensionais do solver. O fluxo físico tende
continuamente ao mesmo limite. A dissipação instantânea de Roe mantém um
salto finito, com a célula vizinha fixa. As velocidades características
máximas dos dois lados tendem, respectivamente, a aproximadamente
`652.201010` e `629.388889`, nas unidades de velocidade adimensionais.

Na implementação anterior, a diferença finita interna de `Dlei` espalhava
esse salto por uma faixa estreita. Por exemplo, a diferença de fluxo numérico
era `1.26722307e-5` para `delta_U3=1e-7`, mas `3.38542835e-11` para
`delta_U3=1e-13`. Isso produz uma dependência muito íngreme da viscosidade
numérica, cuja escala pode ser menor que o passo da diferença finita externa
do Newton. Não constitui uma suavização física deliberada e qualificada.

### 1.4 Consequência para a discretização

Corrigir `Dlei` era necessário para representar a derivada do fechamento em
cada ramo. Essa correção, isoladamente, torna o salto da dissipação mais
explícito; não remove a diferença entre os limites de `|A|`.

O piloto subsequente usou a opção declarada `DISSIPATION_MODE='LAGGED'`: a
dissipação é mantida fixa durante a solução implícita de um passo, enquanto
os fluxos físicos permanecem avaliados no estado corrente. Assim, o Newton
não precisa atravessar uma viscosidade que muda de ramo a cada avaliação
do mesmo resíduo. O piloto ultrapassou 50 s. Esse é um resultado de
robustez local da opção numérica, ainda sujeito a refinamento de malha e tempo.
Conservar o mesmo fluxo em ambos os lados de cada face continua necessário
para o balanço global de massa.

## 2. A interrupção a 78,6073 s: lacuna primal de Bendiksen

### 2.1 Estado observado

No arquivo `pilot5_D_SS3_MATCHED`, o último estado aceito tinha, na célula 49:

| Quantidade | Valor |
|---|---:|
| `s` | 69,1621131186 m |
| Pressão | 207124,602298 Pa |
| Fração de gás | 0,367260166932 |
| Inclinação | 1,31961371869 rad |
| `j_l` | 1,28491342508 m/s |
| `j_g` | 1,22782673830 m/s |
| `Fr=abs(j_l+j_g)/sqrt(gD)` | **3,49999999409** |

A fração de gás variava de 0,355526339 a 0,723830717 no perfil. Portanto, a
falha não ocorre por desaparecimento de uma fase. O último passo rejeitado
tinha `dt=3.32150625e-8 s`, resíduo normalizado `1.02828042e-8` e código de
falha de busca linear 6000. A proximidade do resíduo da tolerância, sozinha,
não caracteriza ruído numérico: aqui há um limite constitutivo demonstrável.

### 2.2 Relação implementada e inversão exata por ramo

A variante de Bendiksen atribuída a Baliño, Burr e Nemoto (2010), equações
32–36, é a ID 2 do pacote. Com `v=sqrt(gD)` e `j=j_l+j_g`, ela usa

\[
u_g=C_0j+U_d.
\]

| Ramo | Condição | `C0` | `Ud` |
|---|---|---|---|
| Baixo Froude | `abs(j)/v < 3.5` | `1.05+0.15 sin(theta)` | `v(0.35 sin(theta)+0.54 cos(theta))` |
| Alto Froude | `abs(j)/v >= 3.5` | `1.2` | `0.35 v sin(theta)` |

O uso de `abs(j)` é a extensão explícita de fluxo com sinal presente no
código; o estado em análise tem ambas as velocidades superficiais positivas.
O diagnóstico, portanto, não depende de extrapolação para fluxo reverso.

Fixadas as duas massas conservativas, a EOS fixa `p`, `rho_l`, `rho_g` e
`alpha`. Para o momento total por volume

\[
m=(1-\alpha)\rho_lu_l+\alpha\rho_gu_g
 =\rho_lj-\alpha(\rho_l-\rho_g)u_g,
\]

o candidato de cada ramo é obtido diretamente:

\[
j_b(m)=\frac{m+\alpha(\rho_l-\rho_g)U_{d,b}}
{\rho_l-\alpha(\rho_l-\rho_g)C_{0,b}}.
\]

Nessa célula, os denominadores de ambos os ramos são positivos. O baixo
ramo é válido somente quando seu próprio candidato tem `Fr<3.5`; o alto
somente quando seu candidato tem `Fr>=3.5`. Não basta escolher o coeficiente
de um ramo antes da inversão e deixar de verificar a condição depois dela.

### 2.3 Limites numéricos e domínio ausente

As densidades no último estado aceito são `rho_l=997.097147689 kg/m3` e
`rho_g=1.76053007079 kg/m3`, com `Cg=343 m/s`. Os coeficientes de ramo são:

| Quantidade | Baixo Froude | Alto Froude |
|---|---:|---:|
| `C0` | 1,19529287435 | 1,2 |
| `Ud` | 0,339746526077 m/s | 0,243388827627 m/s |
| Denominador da inversão | 560,160834860 kg/m³ | 558,440156881 kg/m³ |
| Momento correspondente a `jcrit=2.51274016763 m/s` | 1283,34513944 kg m⁻² s⁻¹ | 1314,24483781 kg m⁻² s⁻¹ |

Logo, neste estado de massas,

\[
1283.34513944\le m<1314.24483781
\]

é um intervalo **sem solução algébrica admissível** do fechamento por partes.
A largura é 30,89969836 kg m⁻² s⁻¹, ou 2,40774655% do limite inferior.
Em variáveis adimensionais, os limites são `2.68154626867` e `2.74611110653`.

O estado aceito tem momento `1283.34513707 kg m-2 s-1`, somente
`2.38e-6 kg m-2 s-1` abaixo da borda inferior. No ponto médio da lacuna:

| Candidato | Froude calculado | Condição que deveria satisfazer |
|---|---:|---|
| Baixo ramo | 3,5384177535 | `< 3.5` — falha |
| Alto ramo | 3,4614638729 | `>= 3.5` — falha |

O probe executável chama a rotina real `compute_primitive_variables`, com as
duas massas positivas preservadas. Ela devolve sucesso no estado aceito,
imediatamente abaixo da borda inferior e acima da borda superior. Dentro da
lacuna, em três pontos independentes, devolve
`PRIMITIVE_NO_ROOT=4`. O estado de tentativa salvo ainda é recuperável: a
busca linear pode recuar para perto da borda, mas não consegue construir a
continuação suave desejada através do intervalo ausente.

### 2.4 O que a evidência permite concluir

A lacuna existe antes de calcular um fluxo de Roe, um Jacobiano numérico ou
um passo temporal. Portanto, nem apertar tolerâncias nem reduzir o passo
cria a raiz ausente. A evidência explica a barreira de admissibilidade
encontrada por esse piloto; não prova que todo problema transiente com essa
correlação deva falhar, nem exclui outras formulações de transição ou
soluções com mudanças de inventário que evitem a lacuna.

Não se deve preencher o intervalo retornando a melhor aproximação, aceitando
somente a largura pequena do intervalo de raiz, congelando permanentemente o
ramo errado ou recortando a velocidade. Isso substituiria a relação física
sem registrar a mudança.

A dissertação cita Bhagwat para o riser. A implementação e qualificação de
uma correlação completa rastreável desse tipo é uma mudança de fechamento
declarada; o caso com ID 2 deve ser preservado como alternativa e como
regressão desta limitação conhecida. Os resultados do novo fechamento não
devem ser misturados aos do piloto ID 2 sem identificar essa diferença.

## Referências e limites de atribuição

1. Dissertação fornecida, capítulo 3, equações de estado 3.10–3.12 e discussão
   dos fechamentos; capítulo 6, dados da instalação e casos. A especificação
   detalhada do benchmark está em `analysis/implementation_reference_spec.md`
   e no JSON correspondente entregues com a auditoria.
2. Baliño, J. L.; Burr, K. P.; Nemoto, R. H. (2010). *Modeling and simulation
   of severe slugging in air-water pipeline-riser systems*. International
   Journal of Multiphase Flow, 36, 643–660.
   DOI: `10.1016/j.ijmultiphaseflow.2010.04.003`, equações 32–36. A extração
   dessa referência foi feita dos textos indexados; a limitação anterior de
   acesso visual a sua tabela experimental de períodos permanece separada
   do diagnóstico constitutivo aqui reproduzido.
3. Código auditado em `scr/module_closure_law.f90`, `scr/module_find.f90` e
   `scr/module_roe_solver.f90`, com hashes registrados nos probes. As
   conclusões numéricas deste documento vêm dos arquivos congelados e dos
   experimentos algébricos citados, e não de períodos simulados tomados como
   dados experimentais.
