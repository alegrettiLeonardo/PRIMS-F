# Auditoria avançada SS1, SS3 e OSC — P0–P2 (20/11/2026)

Branch local `implementation/p0p1-continuation-20261119` (worktree
`~/Documentos/PRIMS-F-p0p1c-20261119`), base histórica `f3717dd`. `main`, PDFs e
resultados anteriores preservados; **sem push/merge**. Relatórios anteriores:
`QUALIFICATION_P0P1_CONTINUATION_20261120.md`, `QUALIFICATION_P0P1_CONTINUATION2_20261120.md`.
Dados brutos em `~/Documentos/PRIMS-F-p0p1c-results-20261119/{e1,p03,p0e1,baro_combined}`.

## Conclusão científica

O bloqueio gasoso de SS1 e do OSC isotérmico tem três componentes distintos.

1. **Instabilidade física do domínio unificado.** O estado estacionário discreto é
   linearmente instável. A partir dele, com ΔP = 0, a solução permanece **bit a bit
   estacionária** por 300 s. Qualquer perturbação cresce exponencialmente (σ ≈ 0,026 s⁻¹,
   período 26,3 s) até o vale do fluxo gasoso na junção chegar a zero.
   - A trajetória até a vizinhança do bloqueio é a **mesma** com Roe e com o fluxo
     experimental E1: Δp_in RMS ≤ 0,6% e A\*(0,25) idêntico, 11,59 kPa.
   - A aproximação ao bloqueio **não** é, portanto, artefato do fluxo de junção.
2. **Operador numérico na junção.** No estado aceito terminal, o Roe anula o fluxo
   gasoso da face com gás presente e subindo nas duas células: classe **C**, com
   u_g a montante igual a 26% de |U_d|. O E1 com traços estacionários leva a dinâmica
   até u_g ≈ 1,3% de |U_d| no fim do oleoduto, perto da estagnação física (**A**).
   O Roe corta cedo, mas apenas 0,2–2,3 s antes.
3. **Representação/fechamento.** A execução termina sempre na recusa do traço do
   riser para j_g = 0 (3021) ou j_g < 0 (3012). O fechamento cocorrente do riser não
   representa a face com gás parado e líquido subindo, que é o início do selo. Essa é
   a extensão E2 (passagem misto→puro), ainda não implementada.

A falta de convergência do instante do evento vinha da **inicialização**: a semente
O(h^0,8) era criada pela amostragem da EDO contínua e pela gaussiana legada. Com IC
estacionária discreta mais perturbação física fixa (ΔP = 250 Pa), a falha cai no
**mesmo ciclo** nas cinco discretizações refinadas (até N401), entre 173,3 e 174,2 s,
convergindo em h para ≈ 173,75 s. Com a IC
legada, a faixa ia de 198 a 276 s. O resíduo é o erro de primeira ordem de σ,
quantizado por ciclo.

SS3 continua BLOCKED. A ausência de raiz é do ramo cocorrente declarado, e existe um
ramo contracorrente contínuo e hiperbólico, mas fora de domínio qualificado.
BAROTROPIC: período e média próximos da convergência; amplitude **não** converge
(erro estimado de 19–22% no nível mais fino). Validação experimental: **NOT_VERIFIED**.

## 1. Diagnóstico j_g = 0 × α_g = 0 (P0.1)

Ferramentas:
- `junction_state_compare.f90` avalia, nos **mesmos** estados esquerdo/direito
  salvos (aceito e proposta recusada), Roe/E1 × CENTER/STEADY. Registra, em SI:
  α_g, ρ_g, u_g, j_g, U₂, p, U₃, G_g de cada célula, a velocidade de deriva U_d do
  fechamento ativo, e na face p\*, G_l\*, G_g\*, α\* de cada lado.
- `junction_classification.py` classifica.

Critério, sem limiar fixo, com razões a escalas físicas da própria execução:
r_α = α/α(t=0), r_u = u_g/|U_d|, r_G = G_g,face / min(G_g das células).
- **B** se r_α < ε;
- **A** se u_g a montante tem |r_u| < ε;
- **C** se o gás está presente e subindo nas duas células e r_G < ε.

A classe é listada para ε ∈ [10⁻¹⁵, 1].

| Estado (aceito, salvo indicação) | Fim do oleoduto α / u_g (m/s) | Base do riser α / u_g | G_g face (kg/m²s) | Classe por faixa de ε |
|---|---|---|---|---|
| SS1 Roe (78,577 s) | 0,560 / 0,0531 | 0,0188 / 0,899 | 2,6·10⁻⁸ (células: 0,0711; 0,0404) | **C** em [6,6·10⁻⁷; 0,25]; A só se ε ≥ 0,26 |
| SS1 Roe, proposta recusada | idem | idem | 6,8·10⁻¹⁵ | C em [1,7·10⁻¹³; 0,25] |
| OSC‑ISO Roe (223,645 s) | 0,5505 / 0,0567 | 0,00565 / 0,915 | 1,1·10⁻⁸ | **C** em [6,5·10⁻⁷; 0,018]; B se ε ≥ 0,018 |
| SS1 E1/CENTER (80,861 s) | 0,559 / 0,0053 | 0,0152 / 0,899 | 8,7·10⁻⁹ | C em [1,2·10⁻⁶; 0,025] |
| OSC E1/CENTER (223,883 s) | 0,549 / 0,0061 | 0,0060 / 0,915 | 4,3·10⁻⁹ | C em [4·10⁻⁷; 0,019] |
| SS1 E1/STEADY (79,886 s) | 0,559 / **0,0026** | 0,0164 / 0,899 | 1,6·10⁻⁹ | C em [4,5·10⁻⁷; 0,013]; **A** em [0,013; 0,3] |

Em nenhum estado aceito o gás desaparece (B exigiria chamar de "ausente" um α de
2–34% do valor inicial). Com o Roe também não está parado (A exigiria tratar
26–27% de |U_d| como zero). A situação é **C**: fluxo numérico nulo. Com E1/STEADY,
o estado terminal fica a 1,3% de |U_d| da estagnação, isto é, aproxima‑se de **A**.
As propostas recusadas aparecem separadas e não entram na classificação física.

## 2. E1 × Roe (P0.2)

### 2.1 Formulação (`007359d`, opção `junction_flux='PHASE_UPWIND'`; default `ROE` inalterado)

O vetor de junção MATCHED efetivo é (G_l\*, G_g\*, p\*). Os fluxos de momento de cada
lado vêm da força pontual existente, `momentum_trace(p*, G*, α*)`, com α\* do
fechamento de cada trecho. O E1 substitui **o vetor inteiro**, não apenas G_g.

- Problema de Riemann acústico linearizado da mistura entre o traço do oleoduto e o
  traço direito mapeado ao mesmo fechamento. Impedância de Wood
  Z = ρ_m c_m, com 1/(ρ_m c_m²) = α/(ρ_g c_g²) + (1−α)/(ρ_l c_l²), média dos lados:
  - P\* = (P_a + P_b)/2 − Z(j_b − j_a)/2;
  - j\* = (j_a + j_b)/2 − (P_b − P_a)/(2Z).
- Deslocamento acústico comum das velocidades de fase, δj = j\* − j̄, que preserva o
  deslizamento médio.
- Separação de fluxo material por fase (limite supersônico das massas tipo AUSM):
  G_k = U_k,a·max(u_k,a + δj, 0) + U_k,b·min(u_k,b + δj, 0).
- Sem clipping de massa, velocidade ou sinal. Fluxo único por fase, compartilhado.
  Para traços iguais, o vetor é exatamente (U₁u_l, U₂u_g, p).

Uma primeira versão (massa a montante × velocidade média) foi descartada pelo próprio
teste: ampliava G_g para 4,3× o fluxo do traço a montante, violando os limites do
mesmo modo que o Roe.

### 2.2 Testes (`test_junction_phase_upwind`, 50 asserções, PASS)

| Teste | Resultado |
|---|---|
| Traços iguais (consistência) | erro relativo 0 em G_l, G_g, p\* |
| Reflexão (lados e sinais trocados) | antissimetria exata de G, p\* simétrico |
| u_g > 0 dos dois lados (salto de 73 Pa) | G_g = fluxo a montante + U₂δj (1,03× o traço); sem amplificação |
| u_g < 0 dos dois lados | G_g < 0, massa a montante (direita) |
| u_g → 0 por ± | G_g → 0 contínuo, sinal acompanha u_g |
| Junção transparente (riser = imagem do traço do oleoduto pelo fechamento do riser) | erro ≤ 8·10⁻¹⁵ (Bendiksen e Bhagwat) |
| Perfil estacionário independente | CENTER: erro de G_g\* 2,0% (Roe 29%), mas G_l\* 1,2% (Roe 0,09%); STEADY: ≤ 4·10⁻⁹ |
| Conservação | telescopagem de massa ≤ 2·10⁻¹⁹ por fase |
| Jacobiano montado × DF central | ≤ 4,8·10⁻⁷ |

O termo acústico do E1 também lê a diferença hidrostática entre centros. A
configuração coerente é, portanto, **E1 + traços STEADY**.

### 2.3 Mesmos estados esquerdo/direito (estados aceitos terminais do Roe)

| Caso | Roe/CENTER | Roe/STEADY | E1/CENTER | E1/STEADY | G_g célula do oleoduto |
|---|---:|---:|---:|---:|---:|
| SS1 G_g\* (kg/m²s) | 2,6·10⁻⁸ | 0,0198 | 0,0686 | 0,0701 | 0,0711 |
| OSC G_g\* | 1,1·10⁻⁸ | 0,0223 | 0,0960 | 0,0980 | 0,1003 |

(Correção registrada: uma avaliação anterior de Roe/STEADY que retornou 1012 usava
um fixture offline com inclinação do oleoduto nula. Com β = 2° o resultado é o da
tabela.)

### 2.4 Replays dinâmicos (binário congelado `007359d`, mesmos decks)

| Referência → E1 | Fim ref. (s) | Fim E1 (s) | aceitos/rejeitados | Δt mín. E1 | CPU ref./E1 (s) | Δp_in máx/RMS | máx \|ΔM_l\|, \|ΔM_g\| (kg) | A\*(0,25) ref./E1 |
|---|---:|---:|---|---:|---|---|---|---|
| SS1 Roe/CENTER → E1/CENTER | 78,577 | 80,861 | 706/29 → 678/23 | 2,5·10⁻⁶ | 213/177 | 2,0%/0,63% | 1,53; 0,0016 | — |
| SS1 Roe/CENTER → E1/STEADY | 78,577 | 79,886 | 706/29 → 669/23 | 1,0·10⁻⁶ | 213/754 | 1,2%/0,39% | 0,95; 0,0010 | — |
| SS1 Roe/STEADY → E1/STEADY | 78,838 | 79,886 | 661/23 → 669/23 | — | 629/754 | 0,96%/0,30% | 0,73; 0,0008 | — |
| OSC‑ISO Roe/CENTER → E1/CENTER | 223,645 | 223,883 | 1829/25 → 1825/24 | 2,2·10⁻⁷ | 311/329 | 1,0%/0,21% | 0,73; 0,0015 | 11 594 / 11 594 |
| OSC‑ISO Roe/CENTER → E1/STEADY | 223,645 | 223,842 | 1829/25 → 1823/24 | 3,7·10⁻⁷ | 311/1581 | 0,54%/0,11% | 0,39; 0,0008 | 11 594 / 11 584 |
| OSC‑ISO N101 Roe → E1/CENTER | 249,400 | 226,307 (*) | 2037/24 → 2423/— | 4,6·10⁻⁵ | 811/— | 0,37%/0,08% | 0,26; 0,0005 | 11 330 / 11 300 |
| OSC‑ISO N101 Roe → E1/STEADY | 249,400 | 226,884 (*) | 2037/24 → 2160/— | 3,2·10⁻⁵ | 811/— | 0,19%/0,04% | 0,14; 0,0003 | 11 330 / 11 320 |

(*) Interrompidas por estagnação de custo, não por falha (`OPERATOR_TERMINATED.txt`).
No vale do fluxo gasoso de ~226 s, que o Roe N101 atravessou (vale de 0,043·ṁ_in,
recuperado), as duas variantes E1 passaram a aceitar passos de 3–5·10⁻⁵ s com
~10 iterações de Newton, avançando ~0,03–0,1 s físico em 20 min. Nesse momento,
ṁ_g na junção era 5–8·10⁻⁶ kg/s e α na base do riser 0,0012–0,0024. Nenhum
TIMEOUT (teto de 4 h) nem FAIL; resultado não terminal.

Todas as execuções E1 concluídas terminam na face da junção com **3021** (fechamento do riser
recusa j_g = 0), não em 3012.

**Critério do pedido:** o E1 remove um defeito demonstrado (o corte Roe com gás a
26% de |U_d|) sem impor transporte gasoso. Ele não desloca o bloqueio de forma
material e não altera a dinâmica nem o mapa amplitude→vale. Não foi promovido.

**Fluxo mais consistente.** Nenhum dos dois foi qualificado contra uma referência
independente do problema de Riemann drift‑flux com junção de fechamentos distintos.
Ela não existe no repositório (NOT_VERIFIED). O E1/STEADY é o único que recupera o
vetor estacionário (4·10⁻⁹) **e** mantém o fluxo gasoso de face dentro dos fluxos
dos traços.

## 3. OSC ISOTHERMAL com perturbação física fixa (P0.3)

Construção:
- **Estado estacionário discreto** de cada malha, pelo mesmo operador e EOS: Euler
  implícito com Δt grande, ponto fixo R(U) = 0 (`SS_N*`). Critério de convergência:
  p_in e ṁ_g constantes entre amostras sucessivas.
- **Perturbação** p(s,0) = p_ss,h(s) + ΔP·f(s), com f = ½[1 − tanh((s − L_p)/w)],
  w = 2 m (sobrepressão no oleoduto), independente da malha.
- **Projeção** (`make_perturbed_ic.py`): com α e u_k fixos, as EOS afins tornam δU
  linear em δp. As médias de célula de f são analíticas, e a integral de δp se
  conserva a 5·10⁻¹⁶ (N51) e 2·10⁻¹⁵ (N101).
- **Erro L2 da projeção constante por célula:** 1,27% (N51), 0,64% (N101), ou seja,
  O(h). O inventário de gás aumenta 0,33% com ΔP = 1000 Pa e 0,082% com 250 Pa.
- Opção `initial_u_file` (opt‑in, `007359d`): EOS, alimentação e geometria continuam
  montadas pelo `main`.

**Controle decisivo.** N51 a partir do estacionário discreto com ΔP = 0 permanece
**exatamente** estacionário por 300 s (semi‑amplitude 0,000 Pa; ṁ_g na junção igual
à alimentação). A oscilação e o bloqueio da campanha legada vinham, portanto,
**inteiramente da IC**: EDO contínua amostrada nos centros mais gaussiana de 1% no
momento amostrada pontualmente. Não vinham de arredondamento nem da propagação.

Campanha (`evidence/p03/P03.{json,md}`). Binário `007359d`, deck B_ISO_MASS e EOS
isotérmica (c_g = 290,0587 m/s) em todas. Definições pré‑registradas em
`p03_analysis.py`:
- regime linear = ciclos com vale ≥ 0,5;
- lei exponencial aceita se o resíduo log for ≤ 0,05.

**ΔP = 250 Pa (série principal):**

| Discretização | A₁ realizada (Pa) | A₀ do ajuste (Pa) | σ (1/s), resíduo log | t(A = 11 kPa) envoltória (s) | E_h(0,9) | E_g(0,25) | E_g(0,01) | Falha de Newton | Vale do evento (rel. aos refinados) |
|---|---:|---:|---|---:|---:|---:|---:|---:|---|
| N51 Δt 0,125 | 689 | 316 | 0,02568; 0,006 | 138,3 | 171,7 | 171,2 | 199,21 | 199,36 | +1 |
| N51 Δt 0,0625 | 725 | 316 | 0,02735; 0,007 | 129,8 | 148,0 | 147,3 | 173,58 | 173,77 | 0 |
| N51 Δt 0,03125 | 744 | 312 | 0,02845; 0,005 | 125,2 | 147,0 | 146,4 | 173,11 | 173,26 | 0 |
| N101 Δt 0,125 | 722 | 314 | 0,02729; 0,007 | 130,3 | 148,4 | 147,8 | 174,01 | 174,20 | 0 |
| N201 Δt 0,125 | 739 | 310 | 0,02840; 0,006 | 125,7 | 147,4 | 146,9 | 173,71 | 173,86 | 0 |
| N401 Δt 0,125 | 749 | 309 | 0,02886; 0,006 | 123,8 | 147,0 | 146,7 | 173,64 | 173,78 | 0 |

**ΔP = 1000 Pa (linearidade e confirmação):** A₁ = 2736 / 2879 / 2954 Pa (N51,
Δt 0,125 / 0,0625 / 0,03125), 2865 (N101), 2934 (N201). A₁/ΔP ≈ 2,9 contra 2,8–3,0
com 250 Pa, o que confirma a resposta linear.

| | N51 Δt 0,125 | N51 Δt 0,0625 | N51 Δt 0,03125 | N101 | N201 |
|---|---:|---:|---:|---:|---:|
| E_g(0,25) (s) | 117,9 | 94,4 | 93,7 | 95,0 | 94,4 |
| Falha de Newton (s) | 123,3 | 121,0 | 120,5 | 121,3 | 121,0 |

Leitura:

1. **Semente.** Com a IC estacionária discreta, A₀ = 310–316 Pa (±1%) para
   ΔP = 250 Pa em todas as discretizações. Com a IC legada, A₁ variava de 43 a
   250 Pa. A hipótese "A₁ ∝ h^0,8 desloca o evento" está **confirmada como efeito
   de inicialização**: o controle com ΔP = 0 não cresce.
2. **Lei exponencial.** A = A₀·e^{σt} vale nos ciclos lineares (resíduo log
   ≤ 0,007). A relação t_ev ≈ t₀ + ln(A_crit/A₀)/σ descreve a **envoltória**.
   O evento observado (vale de ṁ_g) só ocorre no vale seguinte ao cruzamento da
   envoltória, ou seja, é **quantizado por ciclo** (≈ 26,3 s).
3. **Convergência.**
   - σ em h (N51 → N401): 0,02568 → 0,02729 → 0,02840 → 0,02886. Ordem 0,54 no
     primeiro trio e **1,27 no último** (região assintótica de primeira ordem);
     extrapolado 0,02918 s⁻¹; erro de N401 1,1%.
   - σ em Δt (N51): 0,02568 → 0,02735 → 0,02845, ordem 0,6; três níveis, região
     assintótica não demonstrada.
   - Envoltória t(11 kPa) em h: 138,3 → 130,3 → 125,7 → 123,8 s, ordem 1,24 no
     último trio, limite 122,3 s, erro de N401 1,2%. Em Δt: 138,3 → 129,8 → 125,2.
   - Falha de Newton em h: 174,20 → 173,86 → 173,78 s (N101 → N401), diferenças
     0,34 e 0,09 s, limite ≈ 173,75 s.
   - E_g(0,25): 147,8 → 146,9 → 146,7 s.
   - Exceto no nível mais grosso, que fica um ciclo atrás por ter σ 10% menor, as
     cinco discretizações refinadas bloqueiam **no mesmo ciclo**. Falha em
     173,3–174,2 s (dispersão de 0,5%), E_g(0,25) em 146,4–147,8 s (1%). Com
     ΔP = 1000 Pa, falha em 120,5–121,3 s.
4. **Evento ≠ falha.** A falha de Newton vem 0,15–0,2 s depois de E_g(0,01), em
   todas as execuções. E_h(0,9) e E_g(0,25) ocorrem ~26 s antes, no vale anterior.
5. **Diagnóstico.** A não convergência do instante decorria (i) da inicialização,
   agora removida, e (ii) do erro de primeira ordem de σ, amplificado pela
   quantização por ciclo. Não decorria da representação do bloqueio, que só fixa o
   término ~0,2 s após o evento físico. Mesmo com perturbação fixa, o instante pode
   saltar um ciclo enquanto o erro de σ×t não for pequeno frente à distância ao
   vale seguinte. Isso foi observado entre N51/Δt 0,125 e o resto; a convergência
   no mesmo ciclo vale para os quatro níveis refinados.

## 4. Mapa constitutivo de SS3 (P1.1)

`ss3_closure_map.f90`, célula 69 do estado aceito: α = 0,638, p = 228 440 Pa,
θ = 51,64°, G_mix = 0,774 kg/(m²·s), u_l = 0 em u_g = 0,6246 m/s. Varredura
u_g ∈ [−1, 5] m/s, 8001 pontos, cobrindo os dois lados de u_l = 0.

| Fechamento | Domínio admissível | Raízes no estado aceito | Na raiz | Ramo com U₃ → κU₃, κ ∈ [0,5; 1,5] |
|---|---|---|---|---|
| Bhagwat (id 3) | só 0 < u_l < 2,14 mm/s (0,0005 < u_g < 0,6253); recusa 12 dos **dois** lados | 0 | — | raiz cocorrente para κ ≥ 1,025; u_l → 0⁺ quando κ → 1. **O ramo termina na fronteira do domínio**, não numa dobra (sem tangência detectada) |
| Bendiksen (id 2) | todo o intervalo | 1, simples, contracorrente: u_l = −3,74 mm/s, u_g = 1,712 m/s | dΦ/du_g = −0,164; condicionamento 0,45; autovalores −37,4; 1,71; 37,4 m/s (**hiperbólico**) | contínuo para todo κ: u_l de −4,8 a −2,7 mm/s |
| Estratificado de dois fluidos (id 1) | todo o intervalo | 0 | Φ ≈ −0,437 quase constante | sem raiz: a 51,6° o cisalhamento interfacial não sustenta a camada |
| Chexal 1992 | não aplicado | — | — | PDF local, p. 8 (p. 262 do periódico): "Countercurrent flow is only considered for vertical flows" |

Literatura primária. Foram conferidos só registros bibliográficos e resumos; textos
integrais indisponíveis e não lidos:
- Ullmann, Zamir, Ludmer & Brauner (2003), IJMF 29:1583–1604: contracorrente
  estratificado em tubos inclinados, com **duas configurações de holdup estáveis**
  para as mesmas vazões.
- Ullmann & Brauner (2006), IJMF 32:82–105: fechamentos teóricos de parede e
  interface, estratificado liso e ondulado, horizontal e inclinado, com múltiplas
  soluções em inclinados.
- Barnea, Ben Yoseph & Taitel (1986), Can. J. Chem. Eng. 64:177–184: dados de
  início de *flooding* em toda a faixa de inclinações. É modelo de **limite**, não
  de holdup.
- Ohnuki, Adachi & Murao (1988), NED 107:283–294: CCFL em tubo horizontal ligado a
  riser inclinado. Faixas experimentais **NOT_VERIFIED** (página OECD‑NEA recusou o
  acesso, HTTP 403).

Equações, faixas e domínios: **NOT_VERIFIED**.

**Conclusão.** A ausência de raiz é do **ramo declarado** (Bhagwat cocorrente), não
uma inexistência física global. Existe um ramo contracorrente contínuo e bem‑posto
(Bendiksen), mas por extrapolação de uma correlação cocorrente de golfadas, portanto
não qualificável.

**Proposta para continuar SS3:**
1. Obter os textos integrais e os dados de Barnea 1986 e Ohnuki 1988.
2. Implementar no solucionador local um fechamento drift‑flux assinado com ramo
   contracorrente explícito, C₀ e U_d(θ) com domínio declarado, verificando o estado
   de SS3 contra o limite de *flooding* desses dados.
3. Avaliar o modelo de dois fluidos de Ullmann & Brauner, com as duas soluções de
   holdup e seleção documentada.
4. Integrar ao solver global só depois de reproduzir pontos experimentais e com
   aprovação.

## 5. BAROTROPIC: refinamento combinado (P1.2)

A rodada N101/Δt 0,0625 foi verificada no estado real:
- `run_status`: `COMPLETED`, 900 s, 14 418 passos, 0 rejeições, CPU 12 761 s;
- `manifest.json` presente;
- mesmo binário (`restart_baro_813da75`, byte a byte equivalente ao HEAD) e mesma IC
  (`common_initial_U.txt`).

| Execução | Acomoda v2 (14 canais) | Período (s) | Semi‑amplitude (Pa) | Média (Pa) | H₂/H₁ | CPU (s) |
|---|---|---:|---:|---:|---:|---:|
| N51 Δt 0,125 | sim | 23,108 | 9 906 | 266 175 | 0,101 | 2 268 |
| N51 Δt 0,0625 | sim | 23,265 | 13 764 | 266 574 | 0,138 | 6 482 |
| N51 Δt 0,03125 | sim | 23,344 | 15 354 | 266 776 | 0,152 | 10 243 |
| N101 Δt 0,125 | sim | 23,266 | 12 763 | 266 541 | 0,128 | 6 398 |
| N101 Δt 0,0625 | sim | 23,423 | 16 002 | 266 957 | 0,158 | 12 761 |

| Par | Período | Amplitude | Níveis ≤ 1% | Faixas | Formas | Defasagem (fração de T) | Estado |
|---|---:|---:|---|---|---|---|---|
| h (Δt 0,125) | 0,68% | 22,4% | 13/14 | 21–32% | 8–11% | −0,8 a −1,2% | FAIL |
| Δt 0,125 → 0,0625 | 0,67% | 28,0% | 11/14 | 27–41% | 10–15% | −0,8 a −1,2% | FAIL |
| Δt 0,0625 → 0,03125 | 0,34% | 10,4% | 14/14 | 10–18% | 3,8–6,4% | −0,4% | FAIL |
| combinado (N51/0,125 → N101/0,0625) | 1,34% | 38,1% | 10/14 | 36–54% | 14–20% | −1,6 a −2,0% | FAIL |

**Interação espaço–tempo** (cinco pontos; A = A∞ − a·h^p − b·Δt^1,28 + c·h·Δt):
- O efeito de refinar h cai de 2857 para 2238 Pa quando Δt é refinado.
- O efeito de refinar Δt cai de 3859 para 3240 Pa quando h é refinado.
- Interação: −619 Pa, ou seja, 22–27% dos efeitos isolados. A soma não é aditiva.
- Limite estimado da amplitude: 19,7–20,6 kPa (p = 0,87–1,28). Erro do nível mais fino
  (N101/0,0625): **19–22%**.
- Período: limite 23,63–23,71 s, erro de ~1%. Média: 0,2–0,3%.

O ajuste deixa um grau de liberdade, e o eixo h tem dois níveis. A região assintótica
**não** está demonstrada, e não se afirma que algum Δt específico baste. O alvo de 1%
para a amplitude é **FAIL**. A distorção da forma de onda cresce com o refino
(H₂/H₁ de 0,10 para 0,16), coerente com menor dissipação numérica.

## 6. Comparação experimental (P2)

O PDF original de Baliño, Burr e Nemoto (2010) **não está nesta máquina**. Usam‑se os
valores conferidos visualmente no original em 18/11/2026 e registrados em
`docs/benchmark_sources.json` (Tabela 3, p. 653), sem reconferência agora. A medida é
a pressão na base do riser.

| Benchmark | Experimental | Trajetória do modelo | Comparação | Estado |
|---|---|---|---|---|
| SS1, caso 2 | T = 85,6 s; Δp_pp = 0,51 bar; α_oleoduto = 0,55 | bloqueia em 78–81 s, antes do 1º ciclo completo | impossível | **NOT_VERIFIED** |
| SS3, caso 11 | T = 76 s (a dissertação traz 26 s, discrepância mantida); 0,86 bar; α = 0,76 | bloqueia em 97,13 s | impossível | **NOT_VERIFIED** |
| OSC, caso 13 | T = 19,3–33,8 s (modulado); 0,69 bar; α = 0,75 | ISO: crescimento até bloqueio, sem ciclo limite | ISO impossível | **NOT_VERIFIED** |
| OSC, caso 13 (BAROTROPIC, c_g = 343, só indicativo) | idem | ciclo limite T = 23,1–23,4 s, extrapolado 23,6–23,7 s; Δp_pp = 2A = 0,20–0,32 bar, extrapolado 0,39–0,41 bar | período dentro da faixa; amplitude ~0,6× a experimental; EOS diferente do experimento (isotérmico) e amplitude não convergida | **não é validação** |

Não foram usadas previsões numéricas publicadas como alvo experimental.

## 7. Balanços e regressões

O ledger independente (`analyze_accepted_conservation.audit_ledger`) usa fluxos aceitos
com pesos de Euler implícito, por fase, para sistema, oleoduto e riser. Índice completo,
com hashes de entradas e saídas, em `evidence/RUN_INDEX.json` (158 execuções).

| Conjunto | máx \|balanço\|/M₀ | Aritmética do ledger | Descasamento na junção |
|---|---:|---|---|
| E1 (SS1, OSC‑ISO concluídas) | ≤ 4,8·10⁻¹¹ | PASS | 0 exato |
| P0.3 trajetórias dinâmicas (11) | ≤ 4,8·10⁻¹⁰ | PASS | 0 exato |
| P0.3 cálculos estacionários auxiliares (Δt de 1 a 100 s) | ≤ 2,5·10⁻⁸ (N401) | PASS | 0 exato |
| BARO combinado N101/Δt 0,0625 | 1,6·10⁻¹¹ | PASS | 0 exato |
| BARO separados (4) | ≤ 9,4·10⁻⁹ | PASS | 0 exato |

Não há clipping, fonte artificial, pressão ajustada nem tolerância relaxada. `tol_res`
fica em 1e‑10 em todas as campanhas.

| Item (SHA congelado `007359d`) | Resultado |
|---|---|
| Builds release e `check_all` (Fortran 2018) | PASS |
| Suíte Fortran debug (`-fcheck=all -ffpe-trap`) | **25/25 PASS** (+ `test_junction_phase_upwind`) |
| Python `unittest` | **123 OK** |
| CLI | 6/6 PASS |
| Default byte a byte vs `1d71f00` (ISO 30 s; BARO N51/Δt 0,125, 30 s) | monitor, conservation, fields, boundary, checkpoint, rejected e results **idênticos**; difere só o eco das novas chaves na metadata |
| Jacobianos | E1: direcional montado ≤ 4,8·10⁻⁷; junção STEADY e testes anteriores mantidos |
| `tol_res` | `1e-10` em todas as campanhas |

## 8. Commits, binários e ambiente

| Commit | Conteúdo |
|---|---|
| `007359d` | E1 (`junction_flux`), `initial_u_file`, `test_junction_phase_upwind`, `junction_state_compare.f90` |
| `37b60c6` | classificação A/B/C, correção de β no fixture, mapa SS3, construtor da IC P0.3 |
| `d036d43` | evidências E1, P0.3, BARO combinado, RUN_INDEX |
| (este commit) | N401, relatório e handoff |

Ambiente:
- GNU Fortran 13.3.0;
- flags release `-std=f2018 -Wall -Wextra -fimplicit-none -ffree-line-length-none -fopenmp -O2`;
- BLAS/LAPACK: `-llapack -lblas`, resolvido para Intel MKL 2020.4.304 (`libmkl_rt`);
- `OMP/MKL_NUM_THREADS=1`;
- Python 3.12.3, NumPy 1.26.4.

Binários congelados em `frozen_007359d/` (SHA‑256):
- `prims_f` `d1184263…`;
- `restart_baro` `6bf92909…`.

Os binários de campanha anteriores ficam inalterados e documentados.

Comandos:
- estacionário: deck B_ISO_MASS com `perturbation_amplitude=0`,
  `dt_start_s=dt_max_s=1e4`;
- IC: `make_perturbed_ic.py SS_N<N> ΔP 2.0 290.0587009555135 U.txt U.json`;
- execução: o mesmo deck mais `initial_u_file`;
- E1: o mesmo deck mais `junction_flux='PHASE_UPWIND'`, `junction_trace=…`.

## 9. Matriz de estado

| Item | Estado |
|---|---|
| P0.1 classificação A/B/C nos estados aceitos | **PASS**: C (Roe); C→A (E1/STEADY) |
| P0.2 E1 implementado, opt‑in, testado (50 asserções), default inalterado | **PASS** |
| P0.2 E1 remove o corte prematuro sem impor transporte | **PASS** (diagnóstico); dinâmica inalterada |
| P0.2 desbloqueio | **BLOCKED** (3021: traço do riser com j_g = 0; requer E2) |
| P0.2 referência independente do problema de Riemann | **NOT_VERIFIED** |
| P0.3 origem da semente | **PASS**: IC (controle ΔP = 0 bit a bit estacionário) |
| P0.3 instante com perturbação fixa | **PASS**: mesmo vale para os 5 níveis refinados (dispersão 0,5%); em h, falha → 173,75 s, σ assintótico de ordem 1,27 (erro de N401 1,1%); eixo Δt sem região assintótica demonstrada; nível mais grosso um ciclo atrás |
| P1.1 SS3 mapa constitutivo | **PASS** (mapa); SS3 **BLOCKED**; experimento **NOT_VERIFIED** |
| P1.2 BARO combinado concluído e analisado | **PASS** |
| P1.2 BARO amplitude a 1% | **FAIL** (19–22% estimado) |
| P1.2 BARO período a 1% | ≈ limiar (~1% estimado; não demonstrado) |
| P2 validação SS1/SS3/OSC | **NOT_VERIFIED** |
