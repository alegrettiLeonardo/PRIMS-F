# Extensões de Ps+ e jg+ para 7200 s

**Ps+ acomodou; jg+ reprovou o critério de balanço por inventário.** Ambos
terminaram 7200 s, com 57622 passos aceitos e zero rejeições. Os dados de
[3600 s](../static_bc_3600_20261119/README.md) permanecem preservados.

O mesmo operador `69813ee70551fa4a3f075919aeafc45de9693f11` e binário
`7442291fd6cdc311627e79500dc5a4804813b42c7bf4daecf88b4b8b2322288a`
foram reutilizados. Cada execução partiu da mesma prescrição inicial;
somente time_dim/out_dir mudaram. EOS, alimentação, geometria, fechamento,
malha, passo e tolerâncias permaneceram iguais. O novo executor permite até
dois workers externos, um thread por solver e teto global de quatro; isso
não altera max_workers=1 dos decks nem qualquer controle numérico.

[Índice e hashes](campaign_index.json), [métricas](metrics_7200.json) e
[seleção dos seis eixos físicos](selected_sensitivity_index.json).
A seleção usa a maior duração concluída nos controles originais
tol_res=tol_inc=1e-8. Conserva jg+ reprovado e não o substitui por um resultado
com outra tolerância. Cinco comparações periódicas atendem aos critérios
numéricos; isso não constitui validação experimental.
Aqui, “qualificada” significa término, acomodação, balanços e diagnósticos
de contorno atendidos na malha e no passo fixados. Não significa convergência
de discretização das sensibilidades; todas usam 51+51 células e dtmax=0,125 s.

## Ps+: dez ciclos acomodados

A pressão é o traço `p_junction_Pa`. Nos dez últimos intervalos completos:
período médio 23,07714863 s, CV 7,0030e-6, pressão média 271463,987 Pa e
excursão média pico a pico 15874,718 Pa. Acomodação dos canais, contrato
aritmético, seis balanços globais e seis balanços em cada ciclo passaram.
O pior erro global normalizado por inventário é 3,8436e-7. Os quatro lados
de contorno têm 26901 amostras cada, status 0 e contagens 2/2/1/1.

O valor anterior de 3600 s continua rotulado provisório; sua menor dispersão
do período não compensava o crescimento de amplitude ainda observado.

## jg+: limite da aceitação de Newton na duração maior

O detector não encontra dez ciclos acima do piso de 100 Pa; período e
métricas por ciclo permanecem null. O contrato aritmético P2 passa, mas
cinco dos seis balanços normalizados por inventário excedem 1e-6:

| Balanço | Máximo absoluto (kg) | Relativo ao inventário | Critério 1e-6 |
|---|---:|---:|---|
| Líquido, sistema | 1,429914e-3 | 1,998595e-5 | FAIL |
| Líquido, oleoduto | 4,949523e-6 | 8,885966e-8 | PASS |
| Líquido, riser | 1,424964e-3 | 8,992871e-5 | FAIL |
| Gás, sistema | 1,890585e-6 | 1,067604e-5 | FAIL |
| Gás, oleoduto | 4,219947e-6 | 2,737378e-5 | FAIL |
| Gás, riser | 2,329362e-6 | 1,016016e-4 | FAIL |

As seis normalizações por massa processada passam. Isso não apaga as
reprovações por inventário. O desvio acompanha o residual aceito acumulado;
a identidade aritmética não o converte em erro físico aceitável. Os quatro
contornos também passam, com 27074 amostras por lado e status 0: essas
verificações não substituem o balanço integral.

O [diagnóstico da cauda](accepted_tail_observation.json) usa diferenças dos
contadores cumulativos de soluções lineares e conserva hashes dos CSVs.
Desde 4818,93185 s, todos os passos foram aceitos sem nova solução linear;
foram 19049 passos desse tipo no total. Nas 2401 observações entre
6900,05685 e 7200 s, a pressão está constante à precisão gravada. A norma
residual fica entre 5,45e-9 e 9,9941e-9, abaixo da tolerância 1e-8. Esse
patamar numérico não demonstra estabilidade física e não justifica atribuir
um período aos picos residuais abaixo do piso.

O [estudo independente de tolerância](../../../scripts/prepare_static_bc_tol_study.py)
foi autorizado com apenas tol_res=1e-10, mantendo tol_inc=1e-8 e todos os
demais parâmetros. Ele tem nova coorte/selagem, não substitui o FAIL tol8
e não integra uma comparação física de controles numéricos idênticos.
Seu [resultado e limites](../static_bc_tolerance_20261119/README.md) ficam em
pacote próprio, com ambos os manifestos e a reprovação anterior preservada.

As ressalvas de EOS STANDARD/BAROTROPIC, instrumentação alpha_p/P4/P5 e
compatibilidade experimental continuam na [nota de fontes](../../../docs/STATIC_BC_AND_REGIME_CRITERIA_20261119.md).
