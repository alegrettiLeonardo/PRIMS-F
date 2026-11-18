# Limites constitutivos encontrados nos transientes

## Escopo e interpretação

Este documento registra diagnósticos locais das **execuções finais** e, separadamente, dos pilotos anteriores. O propósito é distinguir uma limitação do operador numérico de uma trajetória que atinge o limite do fechamento empregado. Uma execução interrompida não comprova a reprodução de um ciclo de severe slugging. Os tempos abaixo pertencem aos casos identificados; não são alvos experimentais nem resultados de convergência em malha ou passo temporal.

O identificador `closure_riser=3` significa **Bhagwat–Ghajar conforme a transcrição de Azevedo, Baliño e Burr (2015), seção circular, domínio cocorrente**. A fonte efetivamente conferida é a seção 6.3, equações 116–128, p.121 impressa dos [anais oficiais da IMAACA 2015](https://www.msc-les.org/proceedings/imaaca/2015/IMAACA2015.pdf). O [registro da fonte](BHAGWAT_SOURCE_STATUS.md) preserva a identificação da versão consultada e a localização das equações. A conferência integral do artigo original de 2014 continua pendente.

No domínio implementado, as velocidades de fase são não negativas e pelo menos uma é positiva. Uma fase pode estar estacionária, com fração volumétrica interior. Isso é diferente de uma fase ausente. Não há extensão verificada da fórmula para filme líquido descendente ou fluxo de gás negativo. As limitações de contracorrente discutidas na dissertação estão registradas em [BENCHMARK_REFERENCE.md](BENCHMARK_REFERENCE.md), seção 5.

## Execuções finais: três causas distintas

Os resultados desta seção pertencem às execuções `final_`, com perfil estacionário corrigido, ID3, interface `MATCHED`, fluxo Roe e dissipação congelada por passo. A EOS é identificada em cada linha. Cada diagnóstico lê a malha e as escalas do próprio `run_metadata.txt`; os hashes das fontes utilizadas são confrontados com o manifesto da execução. O [resumo estruturado consolidado](../validation/diagnostics/final_domain_summary.json) identifica os arquivos de evidência de cada linha. Nenhum dos casos abaixo atingiu os 2400 s solicitados.

| Caso / EOS | Células oleoduto + riser | Passo máximo (s) | Tempo final alcançado (s) | Status final | Causa localizada |
|---|---:|---:|---:|---:|---|
| SS1 / barotrópica | 51 + 51 | 0,5 | 77,8135719720 | 6000 | Traço gasoso `MATCHED`/Roe fora do domínio cocorrente |
| SS1 / barotrópica | 101 + 101 | 0,5 | 80,7176203174 | 3012 | Traço gasoso `MATCHED`/Roe fora do domínio cocorrente |
| SS1 / barotrópica | 201 + 201 | 0,5 | 83,8347381667 | 3012 | Traço gasoso `MATCHED`/Roe fora do domínio cocorrente |
| SS3 / barotrópica | 51 + 51 | 0,5 | 97,9664102050 | 6000 | Célula do riser no limite `ul=0`; Newton exige raiz fora do domínio cocorrente |
| SS3 / barotrópica | 101 + 101 | 0,5 | 104,8738277208 | 6000 | Célula do riser no limite `ul=0`; Newton exige raiz fora do domínio cocorrente |
| SS3 / barotrópica | 201 + 201 | 0,5 | 109,5081070983 | 3012 | Traço gasoso `MATCHED`/Roe fora do domínio; células continuam admissíveis |
| SS3 / barotrópica | 51 + 51 | 0,25 | 97,1927147578 | 9214 | Caminho de quadratura de Roe com inclinação média fora do domínio |
| SS3 / barotrópica | 51 + 51 | 0,125 | 97,1319590113 | 6000 | Célula do riser no limite `ul=0`; Newton exige raiz fora do domínio cocorrente |
| SS3 / isotérmica | 51 + 51 | 0,5 | 72,6323225104 | 6000 | Traço gasoso `MATCHED`/Roe fora do domínio; células continuam admissíveis |

Esses resultados não autorizam tratar toda interrupção SS3 como retorno físico do líquido. A causa imediata muda com a malha ou o passo máximo. Também não permitem inferir período, amplitude convergida ou ordem de convergência de ciclos de slug a partir dos tempos de parada.

### SS3 final: limite celular em N51, N101 e passo máximo 0,125 s

Os [resumos dos probes SS3](../validation/diagnostics/final_ss3_cases) registram os valores completos. Nos três casos abaixo, todas as células e faces base são admissíveis, a preparação da dissipação congelada retorna sucesso e todos os espectros locais calculados são reais. O impedimento aparece na recuperação celular após a correção de Newton.

| Malha / passo máximo | Célula estendida / índice no riser | Posição desde entrada (m) | Fração de gás | Pressão (Pa) | `ul` (m/s) | `ug` (m/s) |
|---|---:|---:|---:|---:|---:|---:|
| N51 / 0,5 s | 75 / 22 | 62,6891071 | 0,5723646521 | 228 657,9321 | 3,2853 × 10⁻¹⁴ | 0,555070443 |
| N101 / 0,5 s | 174 / 71 | 66,1575337 | 0,5060111507 | 218 038,2659 | 3,9824 × 10⁻¹⁵ | 0,476385676 |
| N51 / 0,125 s | 69 / 16 | 61,2130772 | 0,6382755099 | 228 440,8099 | 3,8695 × 10⁻¹⁴ | 0,624562795 |

Na execução N51/0,5 s, o momento adimensional supera o limite `ul=0` em apenas `2,933113821268485e-14`, enquanto o incremento de Newton no momento é `-7,50705884722369e-9`. Recalculando também o limite para as novas massas propostas, a margem candidata é `-7,541608915270434e-9`. O residual do sistema linear é `3,33770204442631e-24`. As dezoito tentativas de busca linear previstas na produção falham na célula 75 com `PRIMITIVE_NO_ROOT=4`. O probe realiza duas tentativas adicionais apenas para diagnóstico; passos menores podem permanecer do lado admissível, mas não estabelecem uma continuação contracorrente da lei.

Nas outras duas linhas, as margens candidatas também são negativas: `-6,61190801545461e-10` para N101/0,5 s e `-9,617470303436779e-9` para N51/0,125 s. As frações volumétricas interiores afastam a hipótese de desaparecimento de fase. A explicação constitutiva por meio de `x` e `(1-x)^1.5`, apresentada na seção do piloto SS3 abaixo, aplica-se a esses estados finais porque foi novamente verificada nos respectivos snapshots, não porque os tempos do piloto tenham sido reutilizados.

### SS3 final N201: limite do traço da junção

O caso N201/0,5 s para em 109,50810709826744 s com outra causa. A menor velocidade líquida nas células do riser é aproximadamente **0,0650731 m/s**, e todas as primitivas continuam admissíveis. A [análise independente da junção](../validation/diagnostics/final_junction_cases/final_D_SS3_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N201_201_dt0.5/summary.json) identifica a interface 203 → 204.

No estado trial, a pressão comum é 267 922,9433533 Pa; `jl*=0,279476944544 m/s`; e `jg*=5,555801532e-13 m/s`. A fração de gás do traço direito é `8,905572977e-13`, enquanto a menor fração de gás das células é `0,00275400782806`. O fluxo mássico gasoso resulta de `2,8609098419672e-5 - 2,8609098416929e-5 kg/s`, restando `2,743079579e-15 kg/s`.

A correção de Newton mantém todas as primitivas celulares admissíveis, mas produz `jg*=-4,669618865e-10 m/s` e fluxo mássico gasoso de `-2,305542428e-12 kg/s` no traço aproximado. A lei ID3 rejeita o traço com status de domínio 12, propagado como face 1012 e montagem 3012. Como no mecanismo da junção SS1, o sinal da estimativa de Roe não prova retorno físico do gás nas células.

O caso **isotérmico N51/0,5 s**, diagnosticado separadamente em [seu próprio resumo](../validation/diagnostics/final_ss3_cases/final_D_SS3_C3_ISOTHERMAL_MATCHED_LAGGED_ROE_N51_51_dt0.5/summary.json), também atinge o limite do traço de junção, em 72,63232251044695 s. Nesse estado, todas as primitivas e faces base são admissíveis, os espectros locais são reais e a menor velocidade líquida é `0,1058006784 m/s`. Dezesseis das vinte propostas locais examinadas falham na junção com status 1012, preservando a admissibilidade celular. Portanto, a falha do caso isotérmico não deve receber a explicação de velocidade líquida nula observada no caso barotrópico N51.

### SS3 final N51/0,25 s: caminho de Roe com inclinação média

Este caso exige distinguir dois estados. O último snapshot **rejeitado** é de 97,1926957734239 s. O integrador aceitou um passo depois e terminou em **97,1927147577989 s**, com status 9214 durante a preparação da dissipação. Portanto, classificar a parada final apenas pelo último arquivo rejeitado daria uma causa incompleta e um tempo incorreto.

O diagnóstico final usa as primitivas impressas com precisão completa no último bloco de `fields.csv`, reconstrói as conservativas e registra explicitamente essa origem no [resumo separado](../validation/diagnostics/final_ss3_cases/final_D_SS3_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N51_51_dt0.25__final_fields/summary.json). Essa reconstrução reproduz a falha de preparação. Não foi inventado um passo temporal a partir dos campos: os valores de `dt` e norma nesse arquivo de entrada diagnóstico são marcadores zero, identificados como tais.

Todas as 102 células são admissíveis. A menor velocidade líquida é `1,255538159e-7 m/s`, ainda positiva, na célula 73. A face 72 → 73 falha porque o caminho conservativo é avaliado com **uma inclinação fixa igual à média das inclinações das células**, `theta=0,9972091540023094 rad`. A célula direita possui inclinação própria `1,0092051213023341 rad`.

No nó de Gauss de parâmetro `t=0,8872983346207417`, o caminho empregado tem:

\[
U_3=0{,}0015194042751592322,
\qquad U_{3,\min}(\theta_{\mathrm{média}})=0{,}0015206985074773922,
\]

com margem `-1,2942323181599266e-6`. A recuperação não encontra raiz cocorrente. Na extremidade direita reavaliada com a inclinação média, a margem é ainda `-5,204425824320005e-6`, embora o estado seja admissível em sua inclinação própria. O status se propaga por recuperação 4 → quadratura 14 → face 1214 → preparação global 9214.

A [tabela dos nós do caminho](../validation/diagnostics/final_ss3_cases/final_D_SS3_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N51_51_dt0.25__final_fields/roe_path_nodes.csv) inclui uma comparação diagnóstica com inclinação interpolada ao longo do segmento. Nesse teste local, os três nós e as duas extremidades permanecem admissíveis; no nó crítico, a margem passa a `2,79213345973214e-6`. Isso identifica o papel da escolha do caminho geométrico. **Não foi alterado o operador de produção.** Uma mudança exigiria verificar a consistência do fluxo com a dependência geométrica do fechamento e seus próprios testes; esse contraste local não valida automaticamente outro esquema.

### SS1 final: o impedimento de junção persiste nas três malhas

Os [probes SS1 por caso](../validation/diagnostics/final_ss1_cases) confirmam o mecanismo nos estados finais, separadamente do piloto anterior. Em N51, a interface 53 → 54 tem pressão `287927,89367980015 Pa`, menor fração de gás celular `0,0121348655655842` e `jg*=7,722951612e-13 m/s`. A tentativa de Newton mantém as primitivas celulares admissíveis, mas leva o traço a `jg*=-1,86316051647e-9 m/s`.

Em N101 e N201, as menores frações de gás celulares são, respectivamente, `0,00344218738154` e `0,00108947879030466`. Os traços gasosos propostos por Newton tornam-se negativos e são rejeitados pelo mesmo contrato ID3. Os tempos finais 80,7176203174 s e 83,8347381667 s pertencem aos seus próprios snapshots. A persistência desse mecanismo em três malhas demonstra uma restrição recorrente do acoplamento atual; não demonstra convergência de um ciclo SS1 nem desaparecimento do inventário gasoso das células.

### Reprodução das análises finais

O [probe dinâmico](../validation/diagnostics/probe_bg_domain.f90) e seu [script de compilação](../validation/diagnostics/build_final_domain_probe.py) leem as dimensões da malha, as propriedades e as escalas de cada caso. Por exemplo:

```sh
FC=gfortran LIBS='-llapack -lblas' python3 validation/diagnostics/build_final_domain_probe.py \
  validation/results/final_D_SS3_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N51_51_dt0.5
```

Para a parada de preparação no caso N51/0,25 s, use a opção explícita:

```sh
FC=gfortran LIBS='-llapack -lblas' python3 validation/diagnostics/build_final_domain_probe.py \
  --snapshot-source final-fields \
  validation/results/final_D_SS3_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N51_51_dt0.25
```

Cada diretório de saída preserva seu `summary.json`, o estado de entrada, metadados, manifesto, hashes e `build_provenance.json`. O resumo registra tanto o tempo do snapshot quanto o término da execução, evitando confundir um estado rejeitado anterior com o estado final. Os logs dos probes são locais e compactos; nenhuma reprodução inicia um transiente. Os [scripts independentes de junção](../validation/diagnostics/build_bg_ss1_case.py) complementam essa análise para SS1 e SS3 N201.

## Piloto SS3: início de retorno do líquido dentro do riser

### Caso e proveniência

O caso diagnosticado é `bgpilot_D_SS3_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N24_24_dt0.5`: 24 células no oleoduto e 24 no riser, EOS gasosa barotrópica com `Cg=343 m/s`, fechamento estratificado no oleoduto, fechamento ID3 no riser, interface `MATCHED`, dissipação de Roe congelada por passo e limite superior de passo de 0,5 s. Os [metadados congelados](../validation/diagnostics/ss3_bg_stall_metadata.txt) contêm todos os parâmetros.

O snapshot foi copiado do último estado rejeitado do piloto e preservado como [ss3_bg_stall_snapshot.dat](../validation/diagnostics/ss3_bg_stall_snapshot.dat). Seu SHA-256 é `6b6734b7936f205d72c7f1deb355731f7f36bb689e1a497ca7fc04e4aeca7432`. A cópia evita que a escrita posterior do piloto altere o objeto da análise; o [arquivo de hashes](../validation/diagnostics/ss3_bg_stall_snapshot_hashes.json) também identifica os metadados.

Este piloto antecede a correção posterior da contribuição da inclinação no perfil estacionário inicial. Portanto, seu tempo de chegada ao limite não deve ser transferido para um caso reinicializado com a correção. O diagnóstico abaixo é local: reavalia as equações e a direção de Newton no estado congelado, com a implementação disponível no registro de compilação do probe.

### Evidência numérica

O [probe](../validation/diagnostics/probe_ss3_bg_stall.f90) reavalia as primitivas de todas as células, os espectros dos Jacobianos, todas as faces, a montagem com dissipação congelada, o sistema linear de Newton e a busca linear. O [log completo](../validation/diagnostics/ss3_bg_stall_probe.log) e a [tabela por célula](../validation/diagnostics/ss3_bg_stall_cells.csv) sustentam os valores seguintes.

| Quantidade | Valor no snapshot |
|---|---:|
| Tempo físico | 94,48170685397338 s |
| Passo rejeitado | 3,073513716802781 × 10⁻⁸ s |
| Norma do residual antes de Newton | 1,185701433370051 × 10⁻⁸ |
| Código de rejeição do piloto | 6000: busca linear sem passo admissível aceito |
| Célula que atinge o limite | 41, dentro do riser |
| Inclinação da célula | 1,1918637472270275 rad |
| Fração de gás | 0,5257935835525902 |
| Pressão | 222 827,91222934905 Pa |
| Velocidade do líquido | 1,286907777437712 × 10⁻¹⁴ m/s |
| Velocidade do gás | 0,4982050942771629 m/s |
| `C0` | 0,9749007377943918 |
| `Ud` | 0,2428268804963167 m/s |
| Reynolds da mistura na fonte | 15 103,64703211046 |
| Reynolds usado no coeficiente ID3 | 15 419,27768328448 |

Ambas as fases permanecem presentes. A velocidade líquida chega a zero enquanto o gás continua subindo. No lado admissível, todas as recuperações de primitivas tiveram sucesso. Todos os espectros locais calculados pelo LAPACK foram reais: a maior parte imaginária retornada foi zero. Todas as faces base e a preparação da dissipação congelada tiveram status zero. A junção oleoduto–riser também teve sucesso, com fluxo comum de gás positivo, aproximadamente `4,3148514e-4` nas unidades adimensionais do fluxo. Esses resultados excluem, nesse snapshot, uma falha base da junção ou uma perda detectada de hiperbolicidade como causa imediata da rejeição.

O fator `C4` da transcrição ID3 tem uma mudança de sinal associada a inclinação negativa. As células do riser deste caso têm inclinação positiva; o oleoduto usa ID1. Portanto, esse salto de `C4` não está ativo no ponto diagnosticado. Os Reynolds da célula crítica também estão afastados dos limiares 2000 e 2300 da lei de atrito de parede.

### Ausência de raiz cocorrente no passo de Newton

Para as densidades e a fração de gás da célula 41, o probe resolve a mesma lei ID3 com `ul=0` e encontra o limite do momento total adimensional:

\[
U_{3,\min}=0{,}0010366845397163551.
\]

O estado salvo possui

\[
U_3=0{,}0010366845397290843,
\qquad U_3-U_{3,\min}=1{,}2729184026283313\times10^{-14}.
\]

O incremento de Newton calculado para essa célula é

\[
\delta U=\left(
3{,}1352300711581445\times10^{-12},
-1{,}7734242497486271\times10^{-12},
-1{,}0961662827351605\times10^{-8}
\right).
\]

O sistema linear foi resolvido com residual máximo `2,1252102113271311e-24`. A direção tenta reduzir o momento além do limite cocorrente. As vinte avaliações da busca linear, com multiplicadores de 1 até `1,9073486328125e-6`, falharam na recuperação da célula 41: `PRIMITIVE_NO_ROOT=4`, propagado como status 14 do ponto e 1014 da montagem. As massas de ambas as fases continuam positivas; o impedimento é constitutivo.

Essa conclusão não depende apenas da mensagem final 6000. Ela resulta de identificar o ponto que falha, resolver separadamente o estado limite com `ul=0`, conferir o sistema linear e reproduzir as tentativas de busca linear.

### Por que a fórmula não pode ser continuada por recorte ou módulo

A transcrição ID3 usa a fração mássica

\[
x=\frac{\rho_gj_g}{\rho_lj_l+\rho_gj_g}
\]

e o fator \((1-x)^{1{,}5}\) no coeficiente de distribuição. Se o líquido passa a descer, com \(j_l<0\), enquanto \(j_g>0\) e o fluxo mássico total permanece positivo, então \(x>1\). Essa potência deixa de ser real. Substituir o fator por um módulo, limitar `x` a 1 ou prender a velocidade líquida em zero definiria outra lei. A fonte de 2015 consultada não fornece essa extensão de filme contracorrente; sua implementação exige uma formulação física adicional, com sua própria verificação.

A redução do passo não remove esse limite. No integrador atual, um passo suficientemente pequeno pode fazer o residual de `Utrial=Uold` ficar abaixo da tolerância absoluta, permitindo avanço muito pequeno do relógio sem resolver uma continuação contracorrente. O padrão de milhares de passos próximos de `1e-8 s`, frequentemente aceitos na primeira avaliação, não constitui evidência de que a solução atravessou o limite. Não se deve relaxar a tolerância para declarar o transiente concluído.

### Regressão e reprodução

[test_primitives_numerics.f90](../tests/test_primitives_numerics.f90) inclui uma regressão específica desse snapshot: o estado base e seu Jacobiano são admissíveis; a fração de gás permanece interior; as tentativas de Newton com multiplicadores 1 e `1e-4` mantêm massas e pressão positivas, mas retornam ausência de raiz cocorrente, sem fabricar velocidades. A execução local desse conjunto passou **283 assertions**, incluindo dez verificações novas para o limite SS3. Essa contagem se refere ao teste de primitivas/números; não substitui o resultado de uma execução global posterior.

Para reproduzir a análise local, a partir da raiz do repositório:

```sh
FC=gfortran LIBS='-llapack -lblas' python3 validation/diagnostics/build_bg_stall_probe.py
```

O [script de compilação](../validation/diagnostics/build_bg_stall_probe.py) usa apenas a biblioteca padrão do Python, compila em Fortran 2018 com verificações e traps de ponto flutuante, registra o comando e os hashes das fontes em `validation/diagnostics/build_bg_stall/provenance.json` e executa somente o probe do snapshot. Não inicia um transiente. Em um ambiente que fornece somente as bibliotecas versionadas, `LIBS='-l:liblapack.so.3 -l:libblas.so.3'` permite selecioná-las explicitamente.

## Piloto SS1: traço de gás da interface `MATCHED` fora do domínio

### Caso e evidência independente

O segundo diagnóstico usa o snapshot do piloto `bgpilot_D_SS1_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N24_24_dt0.5`, preservado no diretório [bg_ss1_stall](../validation/diagnostics/bg_ss1_stall). O [estado rejeitado](../validation/diagnostics/bg_ss1_stall/last_rejected_state.dat), os [metadados](../validation/diagnostics/bg_ss1_stall/run_metadata.txt) e a [proveniência da cópia](../validation/diagnostics/bg_ss1_stall/snapshot_provenance.json) identificam a execução. O [probe independente](../validation/diagnostics/bg_ss1_stall/probe_bg_ss1_stall.f90), seu [log](../validation/diagnostics/bg_ss1_stall/probe.log) e a [tabela de primitivas](../validation/diagnostics/bg_ss1_stall/point_states.csv) foram produzidos pela auditoria de fechamento. Como no SS3 acima, trata-se de um piloto anterior à correção do perfil estacionário; o tempo não deve ser atribuído a uma execução nova.

O [resumo estruturado](../validation/diagnostics/bg_ss1_stall/summary.json) reúne a identificação e os valores. As [instruções](../validation/diagnostics/bg_ss1_stall/README.md) e o [script de compilação](../validation/diagnostics/bg_ss1_stall/build_probe.py) reproduzem o probe sem executar um transiente.

| Quantidade | Valor no snapshot / perfil trial |
|---|---:|
| Tempo físico | 41,95011050178417 s |
| Passo rejeitado | 7,300434350813846 × 10⁻⁹ s |
| Código de rejeição salvo | 3012 |
| Interface crítica | Células 26 → 27, oleoduto → riser |
| Menor fração de gás entre as células físicas | 0,005079894502560056 |
| Pressão do traço comum | 298 980,84692980215 Pa |
| Velocidade superficial líquida do traço | 0,5296144723724787 m/s |
| Velocidade superficial gasosa do traço | 1,278751303050040 × 10⁻¹² m/s |
| Fração de gás do traço direito | 1,433164309330149 × 10⁻¹² |

As 48 células físicas e os dois contornos permanecem admissíveis, assim como a avaliação base das faces. A fração de gás pequena é a do **traço da junção**, não um inventário de célula que tenha desaparecido. Os estados adjacentes têm frações de gás aproximadamente 0,559798 no oleoduto e 0,00507989 no riser. Suas velocidades gasosas são positivas, aproximadamente 0,112057 e 0,909662 m/s, respectivamente.

O sistema de Newton retorna incremento máximo de `8,6225527690114464e-11`. Após aplicá-lo, todas as primitivas celulares continuam admissíveis, mas a reconstrução do traço `MATCHED` produz

\[
j_g^*=-6{,}800580025319073\times10^{-11}\ \mathrm{m/s},
\]

equivalente a fluxo mássico de gás de `-3,7468913638705913e-13 kg/s`. A inversão da lei ID3 para o traço rejeita esse fluxo negativo com status de domínio 12. O operador de face propaga 1012, e a montagem propaga 3012. Multiplicadores de busca linear 0,5; 0,25; 0,125; e 0,0625 reproduzem o sinal negativo e a mesma falha, sem invalidar as primitivas das células.

### A diferença entre fluxo numérico e retorno gasoso físico

O gás estimado na junção resulta da diferença entre o fluxo central e a contribuição dissipativa de Roe. No estado trial base, ambas são aproximadamente `1,854557944e-4 kg/s`, restando apenas `7,0454905e-15 kg/s`. Na tentativa completa de Newton, a parcela central é `1,8545579414442178e-4 kg/s` e a parcela subtraída é `1,8545579451911091e-4 kg/s`.

Assim, o sinal negativo identificado pertence à **estimativa do fluxo numérico da junção**. Ele não prova que o gás das células passou a escoar fisicamente para trás. O diagnóstico correto é uma incompatibilidade entre o traço produzido pelo acoplamento `MATCHED`/Roe atual e o domínio cocorrente exigido para inverter a lei ID3 no riser. A admissibilidade dos estados celulares, sozinha, não garante a admissibilidade do traço aproximado.

Uma investigação posterior pode formular um acoplamento de interface que produza traços admissíveis, respeitando pressão, conservação de fase e as condições de onda ou de frente do modelo escolhido. Simplesmente trocar o fluxo negativo por zero preservaria a igualdade algébrica dos fluxos de fase se aplicado nos dois lados, mas introduziria uma condição unilateral adicional: isso seria uma nova opção de modelo a justificar e verificar. A alternativa experimental `BLOCKING` existente para outros fechamentos não foi silenciosamente aplicada ao ID3.

### Comparação dos dois impedimentos

| Aspecto | SS3 congelado | SS1 congelado |
|---|---|---|
| Onde a admissibilidade é perdida | Recuperação de uma célula interior do riser | Reconstrução do traço da junção |
| Estado celular antes do passo | Líquido praticamente estacionário; gás ascendente | Fases presentes e velocidades gasosas positivas |
| Tentativa de Newton | Exige momento abaixo do limite cocorrente da célula | Mantém células admissíveis; fluxo Roe do traço torna-se negativo |
| Código localizado | `PRIMITIVE_NO_ROOT=4`, ponto 14, montagem 1014 | Domínio do traço 12, face 1012, montagem 3012 |
| Conclusão permitida | É necessária uma lei que represente retorno do líquido para continuar essa trajetória | É necessário qualificar o acoplamento/traço na região de fluxo gasoso quase nulo |
| Conclusão não sustentada | Que uma tolerância maior resolveria a limitação física | Que o fluxo negativo de Roe comprova retorno gasoso físico |

## Consequência para a qualificação dos casos

O código deve manter a falha explícita fora do domínio verificado. Uma solução física para continuar a trajetória SS3 acima precisa representar a descida do líquido com uma lei apropriada — por exemplo, um fechamento de filme contracorrente ou um modelo com dinâmica adicional devidamente documentado. Para SS1, o traço da interface precisa ser examinado juntamente com o modelo de bloqueio ou de frente escolhido; a falha encontrada não demonstra impossibilidade de todos os acoplamentos compatíveis com ID3. Escolher essas formulações e compará-las com dados compatíveis é trabalho de modelagem; não é uma tolerância numérica a ajustar.

Novos transientes com o perfil estacionário corrigido e refinamentos devem registrar separadamente: tempo efetivamente alcançado, situação de domínio das fases, balanços de massa, passo mínimo e motivo de interrupção. O sucesso dos testes de operadores e a conclusão de outros regimes não demonstram que os ciclos SS1/SS3 foram reproduzidos.
