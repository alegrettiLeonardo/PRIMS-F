# P0: contrato dos contornos — protocolo anterior aos testes

Base local: 2c3e3af8c99d490ffc1dcbf9357da6c0d59d633d; operador de produção
igual a a179e51 (histórico 69813ee identificado separadamente). Branch
audit/boundary-contract-20261119. Nenhuma alteração autorizada em scr/,
defaults, física, main ou remoto. Instrumentação/protótipos são cópias geradas
fora do Git. Todos os hashes das cópias e comandos ficam no manifesto.

## Hipóteses e discriminação

Auditar o fluxo da face, e não inferir sua lei a partir das células fantasma.
Separar: (a) imposição física; (b) balanço com fluxo efetivo; (c) Jacobiano;
(d) domínio constitutivo; (e) globalização. Correlação temporal não demonstra
causalidade. Nenhuma normalização cumulativa substitui a por inventário.

Contrato: jl=Ql0/A; rhog*jg=mdotg0/A; pB=Pt. Definir pB pela EOS do estado
retornado por boundary_flux; nenhuma média com a primeira/última célula.
ml_prescrita=rhol(pB)*Ql0; no limite incompressível =rho_l0*Ql0.
Registrar também o alvo mássico legado rho_l0*Ql0, para distinguir contratos.
Fluxo aceito e peso backward Euler dt são os usados no residual convergido.

## Limiares fixados

* Imposição algébrica: erro relativo de fluxo <=1e-10; pressão
  <=max(1e-7 Pa,1e-12*Pt); vazão zero: erro absoluto e relativo nulo/indefinido,
  respectivamente, sem denominador artificial. Limiares numéricos, não
  incerteza experimental. Reportar mesmo quando a EOS causa diferença pequena.
* Seis balanços (duas fases x global/oleoduto/riser): |erro|/M0<=1e-6.
  Reportar também kg, /massa processada, defeito por passo e cancelamento interno.
  M0=0: normalização indefinida; reportar kg sem piso artificial.
* Jacobianos: diferenças independentes centradas/4a ordem com h relativo
  de 1e-3 a 1e-7; reportar curva de erro, mudança de ramo e arredondamento.
  Comparação principal em estados suaves: erro escalado <=1e-5. Ausência de
  fase exige derivada tangente/status explícito; não exigir derivada transversal
  inexistente nem permitir divisão por h=0.
* Partição artificial: malha, U, theta, largura e fechamento idênticos;
  fluxos/resíduos relativos <=1e-12. Caso de modelos diferentes não é esse teste.
* Acústica: primeiro operador tangente com fontes, sinais e posto das condições;
  pressão fixa pode refletir. Não se aceita absorção como alvo de validação.
* Replays SS1/SS3: snapshots binário64 anteriores preservados; registrar
  status/face/célula, residual e efeito localizado de cada variante. Um replay
  não qualifica dez ciclos e não prova causa global.
* Janela OSC: 120 s a partir de um frame acomodado existente, reconstrução
  comum de U explicitamente identificada, dt_max=0.125 s e demais parâmetros
  idênticos. Limite por execução 180 s de parede. Não extrair período/regime
  qualificado dessa janela nem confundir reconstrução com restart exato.
* Testes de aproximação de fase ausente, drenagem/enchimento e reversão são
  controles locais/curtos, não experimentos de severe slugging.

Nenhuma alteração de tolerância/física para resgatar teste; falhas são evidência.
Custos limitados: replays com até 32 Newton/18 amortecimentos; processos com
timeout registrados como truncamento de infraestrutura, nunca falha física.
