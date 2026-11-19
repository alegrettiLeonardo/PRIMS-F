# Diagnóstico congelado: SS1 com fechamento 3

Este probe lê o estado congelado de uma tentativa rejeitada de `bgpilot_D_SS1_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N24_24_dt0.5`. Não modifica produção, estado nem modelo. O arquivo foi congelado em `snapshot_provenance.json`; o processo original podia continuar escrevendo em seu diretório de resultados.

```sh
FC=gfortran LIBS='-llapack -lblas' \
  python3 validation/diagnostics/bg_ss1_stall/build_probe.py
```

`build/provenance.json` identifica compilador, argumentos e hashes das fontes. `probe.log` e `point_states.csv` contêm os resultados. A configuração dimensional reproduz os valores congelados em `run_metadata.txt`.

## Resultado

O estado está em **41,9501105017842 s**, com tentativa de passo de **7,30043435081385×10⁻⁹ s**. Todas as 48 células, os contornos e as faces do estado congelado são admissíveis. A menor fração celular de gás é aproximadamente **0,00507989450**.

Na junção entre as células estendidas 26 e 27, o fluxo Roe resulta da diferença entre uma parcela central de gás e uma parcela dissipativa quase iguais, ambas de aproximadamente **1,854557944×10⁻⁴ kg/s**. Seu resultado positivo é somente **7,04549×10⁻¹⁵ kg/s**. O traço ascendente correspondente tem fração de gás **1,43316×10⁻¹²**, embora o inventário celular ainda seja finito.

A próxima correção de Newton tem incremento conservativo máximo **8,62255×10⁻¹¹** e mantém todas as primitivas celulares admissíveis. Entretanto, ela faz o fluxo gasoso proposto na junção ficar negativo: **−3,74689×10⁻¹³ kg/s**. A resolução do traço pelo fechamento 3 retorna domínio inválido. As reduções de correção por 1/2, 1/4, 1/8 e 1/16 também continuam do lado negativo nesse estado.

O código da face é `1012 = 1000 + CLOSURE_DOMAIN`; na montagem global resulta em `3012 = 2000 + 1012`. Portanto, o registro `3012` desta ocorrência vem da face interna, e não da pressão de saída. O estado escrito em `last_rejected_state.dat` é a última tentativa armazenada, que ainda pode ser admissível; a proposta de busca linear seguinte explica o diagnóstico de rejeição.

Esta evidência localiza o limite no acoplamento **MATCHED/Roe com traço ascendente cocorrente**. O sinal do fluxo numérico proposto não prova, isoladamente, retorno físico do gás. Aceitar um traço contracorrente, projetar a vazão em zero ou introduzir uma frente de bloqueio são alterações de modelo que exigem formulação e validação próprias.
