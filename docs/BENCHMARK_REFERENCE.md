# Especificação congelada para implementação e testes: catenária da dissertação

Data: 18/11/2026. Este documento reúne a especificação de física, entradas e testes da catenária. Valores verificáveis estão em `benchmark_sources.json`. Não foram alterados os arquivos originais do usuário.

## 1. Decisão de referência e alcance

A referência principal é a **dissertação de Leonardo**, sobretudo cap.3, tabela6.1 (p108), tabela6.2 (p109) e tabela6.7 (p138). A página física do PDF é a página impressa +2. Os casos principais recebem identificadores **D_Stable, D_OSC, D_SS3 e D_SS1**. Os valores originais e suas discrepâncias permanecem registrados; nenhuma correção de parâmetro ou período é silenciosa.

A fonte que efetivamente corresponde à geometria e aos casos citados é **Baliño, Burr e Nemoto (2010), International Journal of Multiphase Flow36,643–660, DOI10.1016/j.ijmultiphaseflow.2010.04.003**. Ela fornece a definição padrão das vazões ausente nas tabelas da dissertação. A fonte experimental original citada é Wordsworth et al.(1998)/CALtec; não foram obtidos os sinais brutos desse relatório.

O PDF novo `04-j.ijmultiphaseflow.2012.06.004.pdf` é **Malekzadeh, Henkes e Mudde(2012)**, outra instalação, com pipeline longo, riser vertical e buffer. Ele não resolve a incompatibilidade geométrica do antigo deck “Montgomery CaseA” e não será combinado com estes casos.

### Acesso e qualidade da evidência

Baliño2010 foi inicialmente lido em duas transcrições indexadas. Em 18/11/2026 o usuário forneceu o PDF completo, permitindo a conferência visual direta. Os links abaixo permanecem como registro de localização bibliográfica:

- [PDF do artigo](https://electronicsandbooks.com/edt/manual/Magazine/I/International%20Journal%20of%20Multiphase%20Flow/2010%20Volume%2036/8/643-660.pdf).
- [Segunda transcrição do artigo](https://www.scribd.com/document/897648646/10-1016-j-ijmultiphaseflow-2010-04-003).
- [Registro da editora](https://www.sciencedirect.com/science/article/abs/pii/S0301932210000807).

O PDF fornecido está identificado no JSON por `project_sources/05-Sci-Hub.-Modeling-and-simulation-of-severe-slugging-in-air-water-pipeline-riser-systems-_-International-Journal-of-Multiphase-Flow-2010.pdf`, SHA-256 `c60944a348bf516f20ed7fd965cff7aa6b8417b07ef3b4b9433c64d28db8a3a5`. Foram conferidos visualmente Eq.78–79,p.649/PDF7, Tabelas3–4,p.653/PDF11 e Tabela5,p.654/PDF12. A restrição anterior de acesso foi superada; o sinal bruto original de Wordsworth não foi obtido. O PDF e imagens de terceiros não são redistribuídos no repositório.

## 2. Geometria e unidades

| Grandeza | Perfil principal D | Perfil comparativo Baliño2010 | Fonte |
|---|---:|---:|---|
| Lp |57.4m|57.4m|Dtab6.1; Baliño§7,p654|
| Inclinação do pipeline no sentido do fluxo |−2°|−2°|mesmas fontes; beta positivo de entrada=2°|
| X do topo do riser |6.435m|6.435m|idem|
| Z do topo |9.886m|9.886m|idem|
| D interno |0.05254m|0.0525018m|valores distintos, não substituir silenciosamente|
| Rugosidade absoluta |4.6×10⁻⁵m|4.6×10⁻⁵m|idem|
| Comprimento extra de buffer |não especificado; não introduzir|Le=0m|Baliño§7|

Catenária da dissertação, eq3.26–3.30, pp66–67:

\[
z=a[\cosh(x/a)-1],\quad Z=a[\cosh(X/a)-1],\quad s_r=a\sinh(x/a).
\]

\[
\theta(s_r)=\arctan(s_r/a),\qquad \theta'(s_r)=a/(a^2+s_r^2).
\]

Resultam aproximadamente a=3.018181937174462m, Lr=12.546253993197072m e θtop=76.47369465598748°. Na malha global, s_r=s−Lp no riser. A junção é uma face com pressão e fluxos de massa contínuos; os fechamentos de cada trecho podem produzir dois valores distintos de α. Pressões são absolutas. Velocidades e gravidade usam o sentido positivo do fluxo do inlet ao separador.

## 3. EOS: duas configurações que não devem ser confundidas

### 3.1 Perfil D declarado

A dissertação, eq3.10–3.12,p62 e tab6.1,p108, especifica:

\[
\rho_l=997.05+(p-101325)/1498^2,\qquad \rho_g=p/343^2.
\]

μl=0.89×10⁻³Pa·s, μg=1.85×10⁻⁵Pa·s, g=9.81m/s². A tabela fixa **Cg=343m/s**. O texto também menciona Cg=√(γRT) para propagação isentrópica. Esses fatos não justificam recalcular Cg em função da temperatura e simultaneamente alegar reprodução literal da tabela. O perfil D usa Cg343 fixo; T serve para identificar o ensaio e para comparações com outra EOS.

A relação ρg=p/Cg² com Cg²=γRT não é a equação de densidade do gás ideal isotérmico p/(RT). A razão γ é pertinente à derivada acústica isentrópica, não uma correção multiplicativa arbitrária na densidade ideal. Uma EOS barotrópica linearizada com intercepto diferente requer outra declaração de modelo.

### 3.2 Perfil comparativo isotérmico

Baliño2010 considera líquido incompressível e gás ideal isotérmico (§2,p645; §7,p654): ρl=1000kg/m³, Rg=287J/(kg·K), μl=0.001Pa·s, μg=1.8×10⁻⁵Pa·s. Para comparação em um solver conservativo com ondas, Cg²=RgT implica a mesma lei de densidade gasosa. **O solver completo com inércia e líquido compressível continua sendo diferente do modelo NPW do artigo**, mesmo que compartilhe propriedades e fechamento.

A tensão superficial não consta na tab6.1. O valor0.0728N/m do repositório deve ser rotulado como hipótese adicional se o fechamento exigir σ; não deve ser apresentado como dado da tabela.

## 4. Vazões: jg padrão, não velocidade local a Pb ou Ps

Baliño2010 eq78–79,p649 define Pstd=1.013bar absoluto e Tstd=20°C:

\[
G_g=\dot m_g/A=\frac{101300}{287\,293.15}j_{g0},\qquad Q_l=A j_{l0}.
\]

A massa de gás de entrada é um dado físico e permanece igual ao mudar a EOS usada para calcular densidades internas. Não usar \(\rho_g(P_b)j_{g0}\) ou \(\rho_g(P_s)j_{g0}\). A velocidade local resulta de \(j_g=G_g/\rho_g(p)\).

Para os casos D, fixa-se Gl=ρl0·jl, isto é, a vazão volumétrica indicada é referida à densidade líquida de referência. Essa convenção completa uma ambiguidade da dissertação e está identificada no JSON; a pequena compressibilidade altera jl local. As massas abaixo empregam a área D=0.05254m. O JSON fornece separadamente as massas calculadas com o diâmetro e a densidade de Baliño2010.

| Caso D | Ps(Pa abs) | T(°C) | T(K) | jg padrão(m/s) | jl(m/s) | Gl(kg/m²s) | Gg(kg/m²s) |
|---|---:|---:|---:|---:|---:|---:|---:|
|D_Stable|420000|35.3|308.45|0.690|0.594|592.24770|0.8307813543|
|D_OSC|206000|20.0|293.15|1.160|0.547|545.38635|1.3966758999|
|D_SS3|201000|19.8|292.95|0.480|0.280|279.17400|0.5779348552|
|D_SS1|208000|20.1|293.25|0.110|0.530|528.43650|0.1324434043|

Dtab6.2 arredonda jlStable a0.59 e jlOSC a0.55; os valores0.594 e0.547 são os apresentados na comparação da própria Dtab6.7. O deck GitHub Stable usava55.4°C/328.55K, além de Ps e D diferentes: ele não representa a linhaStable da dissertação.

## 5. Fechamentos e alcance contracorrente

### 5.1 Intenção da dissertação

Dpp62,65 escolhe **Bhagwat(2014)** no riser e fechamento estratificado no pipeline. A referência bibliográfica correta tem como autores **Swanand M. Bhagwat e Afshin J. Ghajar**, DOI10.1016/j.ijmultiphaseflow.2013.11.001. O nome de função `Swananda` no legado aparenta referir-se ao primeiro nome do autor; não comprova uma terceira correlação independente.

A dissertação não transcreve o conjunto completo dos coeficientes, e o texto integral do artigo original de 2014 não foi obtido. Entretanto, foi recuperado e conferido visualmente um artigo primário de aplicação: **G. R. Azevedo, J. L. Baliño e K. P. Burr, “Linear stability analysis for severe slugging: sensitivity to void fraction correlations”, IMAACA 2015, pp.115–124, seção 6.3, equações 116–128, p.121 impressa / p.135 do PDF**. Ele fornece a transcrição completa usada pela implementação circular `BHAGWAT=3`. A identificação correta é **“Bhagwat–Ghajar conforme Azevedo, Baliño e Burr (2015), duto circular”**; não se alega auditoria integral do original de 2014 nem reprodução exata de todos os detalhes da implementação da dissertação.

O PDF foi obtido de `https://www.msc-les.org/proceedings/imaaca/2015/IMAACA2015.pdf` e salvo em `analysis/IMAACA2015.pdf`. SHA-256: `51891ab67e550804e732603faf37bfac53504ee2d1e36b4f589a4921aefbd171`. A página de equações foi conferida visualmente; consulte a seção 6.3 no PDF da fonte. O PDF e o recorte de terceiros não são redistribuídos no repositório. Esse artigo fornece o fechamento; sua aplicação de estabilidade linear a riser vertical não substitui a geometria ou os dados experimentais da catenária D.

Dpp139–140 discute limitações para SS1/SS2 e para filme líquido descendente/contracorrente. Cobertura angular de−90°a90° não equivale a validação com velocidades de fases de sinais opostos. Inversões jl·jg<0 exigem domínio físico explícito e seleção de ramo; não se deve tornar frações, radicais ou potências admissíveis aplicando `abs` ou recorte sem fundamento.

### 5.1.1 Equações implementadas de Bhagwat conforme a fonte de 2015

Definam-se \(r=\rho_g/\rho_l\), \(j=j_l+j_g\), \(Re_{tp}=\rho_ljD/\mu_l\), \(x=\rho_gj_g/(\rho_gj_g+\rho_lj_l)\) e \(\beta=j_g/j\). A letra beta nessa correlação é uma fração de fluxo volumétrico; não é o ângulo positivo de entrada do pipeline.

\[
C_0=\frac{2-r^2}{1+(Re_{tp}/1000)^2}
+\frac{\left[(1+r^2\cos\theta)/(1+\cos\theta)\right]^{(1-\alpha)/5}+C_{0,1}}
{1+(1000/Re_{tp})^2},
\]

\[
C_{0,1}=0.2(1-\sqrt r)\left[(2.6-\beta)^{0.15}-\sqrt{f_{tp}}\right](1-x)^{1.5},
\qquad
\frac1{\sqrt{f_{tp}}}=-4\log_{10}\left(\frac{\epsilon}{3.7D}+
\frac{1.256}{Re_{tp}\sqrt{f_{tp}}}\right).
\]

Estas são as eqs.116–121. O coeficiente \(f_{tp}\) é resolvido pela equação 118 tal como transcrita. Ele não é substituído pela função de atrito de parede que contém ramos laminar/transição; tal substituição mudaria a correlação de deriva.

As eqs.122–128 dão

\[
U_d=(0.35\sin\theta+0.45\cos\theta)
\sqrt{gD(1-r)}\sqrt{1-\alpha}\,C_2C_3C_4,
\qquad j_g=\alpha(C_0j+U_d).
\]

| Correção ou grupo | Expressão e condição |
|---|---|
| \(\mu_l^*\) | \(\mu_l/(0.001\,\mathrm{Pa\,s})\) |
| \(C_2\) | \([0.434/\log_{10}\mu_l^*]^{0.15}\) para \(\mu_l^*>10\); 1 caso contrário |
| \(La\) | \(\sqrt{\sigma/[g(\rho_l-\rho_g)]}/D\) |
| \(C_3\) | \((La/0.025)^{0.9}\) para \(La<0.025\); 1 caso contrário |
| \(Fr_{sg}\) | \(\sqrt{\rho_g/(\rho_l-\rho_g)}\,j_g/\sqrt{gD\cos\theta}\), avaliado no ramo descendente |
| \(C_4\) | −1 se \(-50^\circ\le\theta<0\) e \(Fr_{sg}\le0.1\); 1 caso contrário |

O domínio inicialmente qualificado no código é: \(\rho_l>\rho_g>0\), \(g,D,\mu_l,\sigma>0\), \(0<\alpha<1\), \(j_l,j_g\ge0\), \(j>0\), \(-50^\circ\le\theta\le90^\circ\). A catenária usa apenas inclinações não negativas no riser. O salto de \(C_4\) no ramo descendente não é alisado; a derivada no ponto de mudança não é inventada. A correlação continua sem qualificação para velocidades de fases de sinais opostos.

A recuperação de \(\alpha\) a fluxos superficiais fixos usa o resíduo completo, recalculando \(C_0\) e \(U_d\) a cada candidato. Como ambos dependem de \(\alpha\), não se usa a fórmula explícita de Bendiksen com coeficientes congelados. Os limites monofásicos exigem inventário ausente declarado pelo chamador; são limites do modelo conservativo e não avaliação da correlação bifásica fora de seu domínio.

### 5.2 Alternativa rastreável: Bendiksen conforme Baliño2010

Eq32–36,p647:

\[
Fr_j=|j_l+j_g|/\sqrt{gD},\quad j_g=\alpha(C_0j+U_d).
\]

\[
Fr_j<3.5:\quad C_0=1.05+0.15\sin\theta,\quad U_d=\sqrt{gD}(0.35\sin\theta+0.54\cos\theta).
\]

\[
Fr_j>3.5:\quad C_0=1.2,\quad U_d=0.35\sqrt{gD}\sin\theta.
\]

A igualdadeFrj=3.5 não está definida: adotar ramo superior é decisão numérica que deve ser declarada. Esta variante difere da fórmula reproduzida por Montgomeryp151, que usa sin²θ, Froude baseado em jl e outra regra de Ud. Não misturar variantes pelo mesmo nome “Bendiksen”.

Perfis executáveis com esta alternativa devem identificar **“parâmetros/EOS D + Bendiksen-Baliño”**. Seus resultados não são reprodução integral da dissertação. Os oito decks anteriores com ID 2 foram preservados em `input/benchmarks/bendiksen/`; os oito perfis principais agora selecionam ID 3, com atribuição explícita à transcrição de 2015. A lacuna algébrica da ID 2 observada perto de \(Fr=3.5\) e a dificuldade distinta de dissipação de Roe estão documentadas em `docs/SS3_STALL_DIAGNOSIS.md`.

## 6. Alvos experimentais e divergências preservadas

| Caso | Alvo declarado em D | Informação do artigo2010 correspondente | Tratamento |
|---|---|---|---|
|Stable|sem Texp|tab5caso1,p654: experimentalSF; modeloStable; sem Texp e sem amplitude numérica|Estabilidade global não significa ausência de slugs hidrodinâmicos ou ruído|
|OSC|Texp33.9s; Tsim3.3s|tab3caso13,p653: Texp19.3–33.8s, sinal modulado; ΔPexp máximo0.69bar|Não usar3.3s como experimento; manter33.9 declarado e contraponto primário|
|SS3|Texp26.0s; Tsim6.9s|tab3caso11,p653/PDF11: Texp76s confirmado visualmente; ΔPexp0.86bar|Leitura primária resolvida:76s. Conflito documental com D26s preservado; não usar6.9s como experimento|
|SS1|não aparece em tab6.7|tab3caso2,p653: Texp85.6s; ΔPexp0.51bar;T19.7°C,jl0.528|Difere de D20.1°C,jl0.530; comparação externa aproximada, identificada|

A inspeção visual confirma que a fonte primária imprime Texp=76s, Tsim=72.2s e erro=−1%. Esses períodos implicam erro de−5%, portanto há uma inconsistência interna na coluna de erro. A dissertação imprime26.0s na Tabela6.7,p.138/PDF140, também conferida visualmente. A origem do erro de transcrição não foi demonstrada;76s é o período reportado pelo artigo primário, enquanto26s permanece preservado como declaração da dissertação.

ΔP é a amplitude da pressão na base do riser, não a queda permanente Pb−Ps. “Max” significa máximo observado em sinal sem amplitude de ciclo clara; não é uma média experimental. As variações ± das condições constam no JSON e são dispersões do ensaio, não uma tolerância automaticamente aplicável ao erro do solver.

Não há neste material uma tabela experimental de picos de produção de gás/líquido no topo para os quatro casos. As curvas da fig20Baliño2010 são **simuladas**. Não transformá-las em dados experimentais. As saídas de produção do novo solver podem ser verificadas por conservação e convergência, mas não qualificadas experimentalmente sem medidas adicionais.

## 7. Protocolo de2400s

Propõe-se total2400s, descarte inicial de900s e análise dos1500s restantes, duração compatível com os registros experimentais citados no artigo (§6,p652,10Hz). Esses tempos de descarte são escolha deste protocolo, não prescrição da fonte. Verificar que as estatísticas já não derivam antes de fixar a janela; se não estiverem estacionárias em sentido estatístico, a execução precisa continuar ou a conclusão permanece indeterminada.

Para SS1 com85.6s, a janela contém aproximadamente17.5ciclos. Para OSC, analisar os dois picos/modulação e não exigir uma senoide única. Para SS3, o período primário reportado é76s, agora conferido visualmente, e a divergência com D26s deve permanecer explícita. ParaStable, testar ausência de uma instabilidade global crescente; pequenas oscilações hidrodinâmicas não constituem, por si, falha.

Registrar em intervalos físicos regulares: Pb,Ps, α dos dois lados da junção, fluxo de cada fase na base e no topo, inventários de massa, balanços cumulativos e contagem de passos rejeitados. Extrair período de picos e espectro da mesma janela; conservar também variabilidade ciclo a ciclo. Confirmar classificação por bloqueio/penetração de gás e produção, além da pressão.

Executar refinamento espacial e temporal separadamente. Oscilação dependente de CFL, passo de saída, tolerância de Newton ou malha, ou acompanhada de erros de massa/estado inadmissível, não demonstra captura física de severe slugging. Um ciclo visível sem convergência não basta. Não calibrar coeficientes ou vazões para reproduzir um período e depois chamar a mesma comparação de validação independente.

A divergência de modelo em relação ao solver NPW e a ausência de auditoria integral do artigo original de2014 permanecem limitações. A página primária de SS3 agora foi obtida e conferida:76s. O conflito com D26s e a inconsistência da coluna de erro do artigo são registrados. Os novos PDFs Chexal1992/Ozaki2015 são fontes para uma auditoria futura; seu recebimento não altera o fechamento implementado. Consulte `docs/CHEXAL_1992_SCOPE.md`.

Os decks principais usam o seguinte protocolo numérico explícito: `tol=1e-12`, `tola=1e-13`, `dt_start_s=0.001`, `dt_max_s=0.5`, `dissipation_mode='LAGGED'` e `junction_model='MATCHED'`; `save_every=0`, `field_interval_s=15` e `output_interval_s=0.1` controlam as saídas. São configurações de solução e observação, não ajustes dos parâmetros físicos ou dos alvos experimentais. O passo de tempo continua adaptativo; `dt_max_s` é um teto. As massas impostas, Ps, temperaturas, diâmetro, viscosidades, \(C_g=343\), pares isotérmicos e hipótese adicional de \(\sigma=0.0728\) permanecem documentados e inalterados. Os decks Bendiksen arquivados preservam o conteúdo anterior, incluindo seu protocolo de saída; ao executá-los, deve-se escolher uma saída própria para distinguir as alternativas.

## 8. Implementação estacionária entregue nesta etapa

Foi substituído apenas o módulo estacionário do checkout de trabalho `PRIMS-F-implementation/scr/module_steady_state.f90`, além de seu teste `tests/test_steady.f90`. A assinatura antiga foi preservada e os parâmetros opcionais acrescentados. `face_positions` e `face_pressures` retornam os mesmos espaços adimensionais usados por `vns` e `vnp`.

A integração emprega fluxos de massa constantes e a identidade

\[
H=G_lu_l+G_gu_g+\Pi_Pp,\qquad
\frac{dp}{ds}=\frac{S_3-H_\theta\,d\theta/ds}{H_p}.
\]

Os termos de derivada são calculados sobre o próprio fechamento, com a pressão explícita diferenciada analiticamente. A derivada angular usa stencil unilateral de segunda ordem na base do riser (e junto ao limite vertical), mantendo todas as avaliações no domínio geométrico do riser. Isso evita atravessar o salto real de C4 da correlação de Bhagwat em θ=0 para baixo Frsg. O teste do traço direito da base falhou no stencil central anterior e passou com duas resoluções independentes após a correção; a lei de C4 não foi suavizada. A fonte de atrito usa o fator de Fanning e o multiplicador\(-2L_r/D\). O estado inicial é integrado da **face** de saída, com Ps imposto nessa face; pipeline e riser são integrados separadamente para não cruzar a descontinuidade do fechamento em um estágio de Runge–Kutta. A razão geométrica e o número independente de células do riser são respeitados.

`solve_alpha_from_superficials` é público. Velocidades superficiais, densidades e velocidades de referência seguem o contrato adimensional existente. Raiz ausente, domínio inválido e correlação não suportada produzem status e NaN; nenhum ponto de menor resíduo é devolvido como se fosse solução. A ausência de fase conhecida pelo inventário é declarada com `phase_hint` (1:líquido puro;2:gás puro), produzindoα=0 ou1 exatamente. Fluxo de fase zero sozinho não indica ausência: o fechamento estratificado preserva a raiz de gás estagnado, inclusive na passagem porjg=0.

O ramo estratificado cocorrente descendente com Reynolds líquido mínimo acima de2200 usa bracket nas extremidades e Brent protegido. A varredura completa em logit permanece disponível: é executada para contracorrente, estados com possível transição interfacial, para Bhagwat dependente de \(\alpha\) e sempre que `nroots` é solicitado. Múltiplas raízes retornam status22 por padrão. O argumento opcional `alpha_hint` permite continuação explícita: seleciona a raiz física mais próxima em distância logit; hints0 e1 selecionam os respectivos limites. Empate permanece ambíguo. Uma varredura finita identifica raízes com mudança de sinal; não prova ausência de uma raiz tangente degenerada. Essa escolha de ramo não qualifica experimentalmente uma correlação contracorrente.

A hipótese de acoplamento P contínua/α com dois traços vem da referência. Ela não implica que as parcelas cinéticas de H sejam iguais nos dois lados; o tratamento discreto da junção deve incorporar a hipótese de acoplamento consistentemente. Os testes do módulo verificam cada trecho e seus traços, sem alegar continuidade automática do fluxo de momento total por impor somente continuidade de P.

**Verificação executada do módulo estacionário:**398 verificações passaram com compilação Fortran2018, verificações de acesso e traps de ponto flutuante ativadas. Incluem solução hidrostática compressível analítica, pressão na face versus centro, duas fases com invariantes de massa, independência da pressão em relação à malha de saída, contribuição explícita da catenária em dH/ds, concordância de raízes rápidas e varredura em oito pressões para cada um dos quatro casos, duas raízes contracorrentes estratificadas e convergência de diferenças finitas em cada ramo selecionado. As 19 verificações mais recentes conferem a identidade dH/ds=S no traço direito da base, continuidade local, roundoff da face e convergência das pressões de face entre duas malhas. A extensão ID 3 acrescentou inicialização dos quatro perfis D, vinte raízes em inclinações distintas, invariância às escalas adimensionais, preservação das massas e continuidade da raiz e de sua derivada na vizinhança do antigo limiar de Froude da ID 2. Esses são testes da implementação; não substituem os ensaios transientes experimentais definidos acima.

## Atualização documental e rastreabilidade

A atualização com os novos PDFs altera somente metadados de referência e a leitura dos alvos reportados; não altera os parâmetros físicos D, EOS, namelists, código ou resultados da campanha. Os snapshots de referência utilizados nas execuções permanecem em `validation/evidence/benchmark_sources_at_run.json` e `validation/evidence/benchmark_sources_at_refined_run.json`. O último conserva o hash `8c81309e2d584304e6876528ee766007c0342a40ed841c048affe89f285f3d95` utilizado nos refinamentos.
