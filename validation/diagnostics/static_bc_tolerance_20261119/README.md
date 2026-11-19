# jg+7200: estudo isolado de tol_res

**Com tol_res=1e-10, os seis balanços passaram em 7200 s.** O
[FAIL anterior com tol_res=1e-8](../static_bc_7200_20261119/README.md)
permanece preservado e continua no índice dos seis eixos físicos.

O [gerador](../../../scripts/prepare_static_bc_tol_study.py) e a revisão
independente confirmaram que somente `tol_res` e `out_dir` mudaram no deck
e nos metadados efetivos. `tol_inc=1e-8`, duração 7200 s, passo, malha,
controlador, EOS, alimentação, geometria e inicialização permaneceram iguais.
Operador `69813ee70551fa4a3f075919aeafc45de9693f11`, binário
`7442291fd6cdc311627e79500dc5a4804813b42c7bf4daecf88b4b8b2322288a`.
Essa é uma nova integração desde a mesma prescrição inicial, não um restart.

[Índice e hashes](campaign_index.json), [comparação completa](tolerance_comparison.json)
e [diagnóstico da cauda](accepted_tail_observation.json). Os arquivos ligam
os dois manifestos, seus CSVs e as mesmas versões dos auditores. A revisão
conferiu todos os hashes de entrada do P2 anterior; nenhum FAIL foi apagado.
Os [contadores completos](solver_counters.json) preservam passos, rejeições,
iterações Newton, avaliações de Jacobiano/residual, soluções lineares,
tentativas de backtracking, preparação e tempos de CPU/parede de ambos os casos.

| Balanço relativo ao inventário | tol_res=1e-8 | tol_res=1e-10 |
|---|---:|---:|
| Líquido, sistema | 1,998595e-5 | 6,664490e-9 |
| Líquido, oleoduto | 8,885966e-8 | 2,060307e-10 |
| Líquido, riser | 8,992871e-5 | 2,936743e-8 |
| Gás, sistema | 1,067604e-5 | 2,009933e-9 |
| Gás, oleoduto | 2,737378e-5 | 7,339562e-9 |
| Gás, riser | 1,016016e-4 | 3,382714e-8 |

Critério: cada uma das seis combinações deve atender 1e-6 tanto por
inventário quanto por massa processada. O caso tol8 reprova cinco critérios
por inventário; o estudo tol10 passa todos. Ambos passam o contrato
aritmético, que é uma verificação separada. A nova execução completou
57622 passos, zero rejeições, CPU 990,341 s e parede do executor 993,151 s
(parede interna do solver: 992,969 s).
Os quatro lados de contorno têm 27074 amostras cada, status 0, contagens
características 2/2/1/1 e nenhuma velocidade de fase negativa.

A mudança também deslocou o início da cauda sem nova solução linear:
de 4818,93185 s e 19049 passos no caso tol8 para 7136,80685 s e 506 passos
no estudo tol10. Nos últimos 300 s observados, a excursão da pressão de
junção é 0,00036744 Pa, muito abaixo do piso de 100 Pa do detector. Não há
dez ciclos resolvidos nem período qualificado. Os balanços por ciclo
permanecem sem intervalos, sem aprovação vazia e sem inferir estabilidade
física pela ausência de oscilação observável.

O estudo demonstra sensibilidade da acumulação do residual à tolerância
de aceitação nesta condição e duração. Não modifica a regra de aceitação,
não prova conservação exata ou controle do erro em durações arbitrárias
e não garante o mesmo resultado para outras condições. Não foi usada
mudança de parâmetro físico para produzir ou suprimir ciclos.

O resultado tol10 não substitui um dos casos da comparação de contornos com
controles numéricos idênticos. Convergência de malha/passo e validação
experimental continuam sendo questões separadas, com as limitações de
EOS e observação registradas na [nota de fontes](../../../docs/STATIC_BC_AND_REGIME_CRITERIA_20261119.md).
