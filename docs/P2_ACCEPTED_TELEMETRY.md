# Telemetria conservativa por passo aceito — P2

Base: `9f2ed21`. A captura de fluxo/contadores reutiliza o desenho de
`4a35c1f` (coorte `14191fc`), ampliado para os dois lados da junção e os
residuais de massa. Não muda residual, Jacobiano, física, domínio, tolerâncias
ou critérios de aceitação; não altera os campos nem o diagnóstico Burr.

`steps.csv` tem uma linha por passo aceito: índice, tempo e dt em segundos,
iterações da tentativa aceita, rejeições acumuladas, norma aceita, totais de
iterações Newton, avaliações de Jacobiano/residual, soluções lineares,
tentativas de backtracking e preparação da dissipação, CPU e tempo de parede.
Avaliações de residual incluem aquelas com Jacobiano. O custo até uma aceitação
inclui tentativas rejeitadas desde a aceitação anterior; `run_status.txt` registra
também o custo e os contadores terminais, inclusive falhas posteriores.

`conservation.csv` contém uma linha inicial e uma por aceitação. Inventários
L/G são somas de U com as larguras FV exatas, convertidas em kg, para sistema,
oleoduto e riser. As integrais usam `dt` aceito e os fluxos de entrada/saída e
de ambos os lados da junção capturados da **mesma montagem que convergiu**,
inclusive quando a aceitação ocorreu no line-search. Nenhum fluxo de tentativa
rejeitada contribui. O fluxo inicial é diagnóstico, não contribui às integrais;
se indisponível, conserva NaN e status explícito.

Orientação positiva: entrada → oleoduto → riser → saída. Para cada fase:

- Sistema: `E = M - M0 - I_in + I_out`.
- Oleoduto: `E_pipe = M_pipe - M_pipe0 - I_in + I_base_left`.
- Riser: `E_riser = M_riser - M_riser0 - I_base_right + I_out`.
- Resíduo FV ponderado acumulado do sistema:
  `sum(R_mass) = E + I_base_left - I_base_right`, até arredondamento.
- Por trecho, o resíduo ponderado é seu próprio erro de balanço.

Integrais absolutas somam `dt*abs(fluxo)` somente para medir massa processada;
os fluxos físicos e as integrais de conservação permanecem assinados. O resíduo
Newton medido é um observável: não é subtraído para fazer o balanço passar.

`monitor.csv` mantém as colunas existentes e acrescenta acumuladores compatíveis
com a instrumentação anterior. Os fluxos de fase reportados após falha são os
últimos aceitos. As colunas de pressão/holdup do traço e `blocking_active`
continuam sendo diagnósticos recalculados; não constituem um histórico de todos
os traços aceitos.

`scripts/analyze_accepted_conservation.py` faz leitura streaming e verifica
índices, tempos/dt, sinais, inventários, integrais e identidade com o residual.
Associa `steps.csv` e `conservation.csv` por hash dos pares índice/tempo/dt/norma,
verifica o status final e registra hashes dos arquivos e do analisador. Recusa
sobrescrita. Na auditoria v1, o limite fixo de 256 eps foi ultrapassado após
2718 adições de fluxo gasoso constante: diferença de 3,45e-14 kg entre a soma
simples Fortran e a soma compensada independente. Os JSONs v1 permanecem intactos.

A auditoria v2 usa o limite de recorrência `gamma_n = n*u/(1-n*u)`, com
`u = epsilon/2`, multiplicado pela soma das magnitudes das contribuições do
acumulador correspondente. Isso decorre do modelo de arredondamento relativo
por operação; pressupõe aritmética finita normal e `n*u < 1`. O termo local
256 eps continua separado e é acumulado nas identidades de defeito por passo;
não constitui uma prova geral de erro para todo operador físico. Os limites e
erros máximos ficam no JSON junto do resultado do antigo limite fixo. O critério
físico de balanço 1e-6 é separado, em inventário inicial e massa absoluta processada;
escala nula é reportada indisponível, nunca substituída por valor artificial.

Testes específicos: 15 testes Python de sinais/defeitos/contratos e cinco
contratos CLI curtos (saída esparsa/frequente, término por limite de passos,
rejeição real e estado final invariável com amostragem). A fixture Fortran FV
inclui o telescópio por trecho em malha não uniforme, os dois fluxos de junção
e identidade dos operadores com/sem saída opcional. Logs locais em
`tmp/p2_cli_v3`, `tmp/p2_cli_freeze`, `tmp/p2_unit_tests` e
`/tmp/primsf_p2_python_final.log` e `/tmp/primsf_p2_gamma_tests.log`. O parser preserva também o campo textual
`termination_reason` da telemetria P1; 5/5 smokes foram repetidos com o binário
congelado `69813ee`, sem alterar esse solver.

```sh
python3 -m unittest discover -s tests -p test_accepted_conservation.py -v
python3 tests/test_accepted_telemetry_cli.py --binary bin/prims_f --output tmp/p2_cli_new
python3 scripts/analyze_accepted_conservation.py --run-dir <run> --output <novo.json>
```

Os testes verificam implementação e conservação discreta; não são validação
experimental ou classificação de regime. Revisão independente de sinais,
unidades e montagem aceita não encontrou bloqueante.
