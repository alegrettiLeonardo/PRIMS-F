# Integração Chexal 1992 — contrato do solver

## Seleção e escopo

`closure_riser=5` seleciona `CHEXAL_1992` ar–água. IDs `6` e `7` selecionam,
respectivamente, o ramo de maior e de menor fração de vazio da mesma correlação.
São escolhas explícitas do usuário/modelo; não são um critério experimental para
decidir qual ramo ocorre. ID5 recusa selecionar implicitamente um ramo ambíguo.
Bhagwat permanece ID3 e mantém seu domínio cocorrente anterior.

O núcleo recebe SI e converte o ângulo do PRIMS-F (radianos desde a horizontal)
para graus desde a vertical. A inversão de orientação conserva os sinais dos
fluxos; Reynolds não é tornado positivo indiscriminadamente. A contracorrente
permitida é vertical, gás ascendente e líquido descendente. A catenária D não
contém seção vertical; o traço inferior MATCHED é horizontal. IDs6/7 não ampliam
esses domínios.

## Operadores e unidades

O residual é `R=alpha*(C0*(jl+jg)+Vgj)-jg`, em m/s, inclusive em `jl+jg=0`.
O adaptador divide R pela escala de velocidade. Para as derivadas em relação a
`[alpha,rhol,rhog,ul,ug]`, as escalas são, respectivamente,
`[1/wu,wr/wu,wr/wu,1,1]`. Densidades e velocidades físicas são reconstruídas
antes da avaliação de coeficientes; propriedades, gravidade e diâmetro são SI.
Área inconsistente com o diâmetro circular é rejeitada.

A recuperação conservativa conserva exatamente a identidade do momento e usa
o mesmo ramo selecionado no fechamento. Uma dica de velocidade só é enviada ao
núcleo quando o chamador fornece o estado anterior; a velocidade média de uma
mistura não é evidência para escolher um ramo. Os Jacobianos de fluxo usam a
mesma equação implícita de velocidades e a mesma EOS da recuperação.

As derivadas Chexal são diferenças de quarta ordem por extrapolação Richardson,
com redução de incremento e verificação da assinatura das peças constitutivas.
Para velocidades não nulas muito pequenas, um stencil unilateral de quarta
ordem aponta para o interior do domínio. Uma guarda bilateral adaptativa exige
a mesma peça em torno do estado; o stencil também preserva essa assinatura.
Isso remove o piso numérico artificial de aproximadamente 2,2586e-8 m/s do
stencil centrado anterior, mantendo a recusa em zero e nas quinas exercitadas.
Se não é possível obter um stencil na mesma peça ou ele não converge, o operador
retorna falha de derivada; esse status não é prova de domínio físico excedido.
Consulte `CHEXAL_VELOCITY_DERIVATIVE_REFERENCE.md` para a referência independente
e a regressão da guarda perto de Rel=Reg. Não se arredondam autovalores complexos para reais. O caminho
Roe com geometria variável é documentado e testado separadamente; sua matriz
é uma dissipação aproximada, sem alegação de propriedade Roe exata ao omitir
o termo geométrico explícito da identidade de fluxo.

## Fases ausentes e erros

Massa conservativa exatamente nula estabelece inventário ausente e permite
recuperação monofásica exata. Fluxo superficial nulo isolado não estabelece isso.
Na inversão estacionária, um endpoint algébrico exige `phase_hint` explícito;
caso contrário retorna `STEADY_PHASE_REQUIRED=25` e alfa inválido. O Jacobiano
bifásico completo no endpoint de fase ausente permanece indefinido e é recusado;
não se cria um terceiro modo artificial para passar pelo limite.

Status do núcleo: `60` entrada inválida, `61` domínio não coberto,
`62` informação de ramo necessária, `63` raiz não encontrada,
`64` limite CCFL excedido, `65` singularidade e `66` falha numérica.
Os status propagam pelos adaptadores; códigos compostos dos operadores FV
preservam o motivo da falha. `NO_ROOT` não é prova analítica de inexistência.

## Qualificação

Os testes de integração verificam ida/volta U–primitivas, dimensões e escalas,
derivadas de fluxo contra perturbações independentes de U, hiperbolicidade dos
estados exercitados, contornos, fluxo comum e ausência de fontes de massa.
Referências independentes do núcleo, CCFL e ramos são mantidas separadamente.
Esses resultados são verificação numérica; validação experimental requer a
campanha e a comparação das grandezas físicas, com a mesma EOS e geometria.

A seleção Chexal não resolve automaticamente o traço gasoso negativo de
MATCHED/Roe. Gás negativo na base horizontal permanece fora do domínio e não
é substituído por zero. Bloqueio gasoso completo não é requisito universal de
OSC/SS3; produção, penetração, inventários e balanços devem ser avaliados junto
com a pressão. A divergência SS3 26s (dissertação) versus76s (Baliño2010,
p.653, Tabela3, caso11) permanece identificada, sem calibração de parâmetros.
