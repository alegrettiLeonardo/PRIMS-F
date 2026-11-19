# Auditoria independente das entradas da comparação final

Registro em **18/11/2026, 20:34:31 UTC**. Conferência somente de leitura; nenhum
solver ou teste foi reexecutado. O registro não atesta término, acomodação ou
validação experimental dos casos.

Worktree conferido:
`/home/leonardoab/Documentos/PRIMS-F-comparison-final-20261118`.
HEAD congelado: `490936a676c2a86f1348f0c867bda65685c27823`.
O estado Git não tinha alterações em arquivos rastreados; os decks e planos da
campanha eram arquivos novos. Este README fica no worktree de implementação.

## Arquivos e procedimento conferidos

- `validation/comparisons/comparison_final_20261118/plan.json`, iniciado em
  `2026-11-18T20:25:13.764078+00:00`: D_SS1, D_SS3, D_OSC e D_Stable, C3/C5,
  51+51 células, teto de passo 0,125 s, duração solicitada 2400 s.
- `validation/comparisons/comparison_final_20261118_h/plan.json`, iniciado em
  `2026-11-18T20:26:48.770794+00:00`: D_OSC, C3/C5, 101+101 células,
  teto de passo 0,125 s, duração solicitada 2400 s.
- Os dez `validation/results/comparison_final_20261118*/manifest.json`, os
  dez `run_metadata.txt` correspondentes e os dez decks em `validation/cases/`.
  Os hashes de namelist/referências declarados nos manifestos coincidem com
  os arquivos. As entradas efetivamente registradas pelo Fortran também foram
  comparadas, incluindo os valores derivados de geometria, vazões e escalas.
- As 23 fontes `scr/*.f90` coincidem com os bytes de `git show 490936a:<caminho>`.
  Também coincidem Makefile, `scripts/run_validation.py`,
  `scripts/run_closure_comparison.py`, `docs/benchmark_sources.json` e os quatro
  decks originais `input/benchmarks/D_{SS1,SS3,OSC,Stable}.nml`.
- `tmp/corrected_build.log` registra o build release Fortran 2018/O2. Os dois
  planos e os dez manifestos apontam para os mesmos hashes de fontes e binário;
  o hash do executável presente coincide com eles. Não foi feito novo build
  para tentar demonstrar reprodutibilidade binária.

## Resultado da comparação

| Pares C3/C5 | Malha | dt máximo | Diferenças nos decks e metadados efetivos |
|---|---:|---:|---|
| D_SS1, D_SS3, D_OSC, D_Stable | 51+51 | 0,125 s | Somente `closure_riser` e `out_dir` |
| D_OSC | 101+101 | 0,125 s | Somente `closure_riser` e `out_dir` |

EOS BAROTROPIC, propriedades ar–água, vazões líquida/gasosa, pressão,
geometria, perturbação de 0,01, condições padrão, escalas e controles numéricos
restantes coincidem dentro de cada par. A igualdade é das entradas e da
amplitude de perturbação; não exige que fechamentos distintos produzam o mesmo
estado estacionário inicial ou a mesma trajetória.

Para cada fechamento, o refino espacial D_OSC 51→101 muda somente `nodes`,
`n_riser` e `out_dir`; o teto de passo permanece 0,125 s. O passo adaptativo
efetivamente aceito ainda pode depender da solução e das restrições do solver.
Nenhum desvio ou mistura de coortes foi encontrado nos arquivos disponíveis.

## Identificadores SHA-256

| Artefato | SHA-256 |
|---|---|
| `bin/prims_f` | `290280e37841d3f768299035afc668590474871f066c3200abec341f198b1bb7` |
| `scr/module_chexal_1992.f90` | `63a60d888a288eee582ab61d418cfe8838d1aa0e79b145e8a163d25c026d0443` |
| `scr/module_closure_law.f90` | `2c7e9b7d81651faa6580be044fcf04e8430fc54670afd8b10dd614d2e955d233` |
| Plano base | `fdf6fa922254f86671dc4fc12716c4ae3d29824cf71cbce135962f4414d05002` |
| Plano de malha `_h` | `ad798f4fc539c116d395965e59e27b4c2f1afe10dc65657163a214c462f23003` |
| `tmp/corrected_build.log` | `0763f734a7a7cea1247629eaee1aeabf903949e2199fdb4a35339f84b6182ba4` |
| Helper independente de leitura | `5a7f2247da20bced8b2b9df2ced48ca45420abd062f74b59acea59e4a3850039` |

Os planos conservam o inventário completo dos 23 hashes de fontes. O helper
local usado foi `audit_final_comparison_readonly.py`, em
`C:/Users/leonardoab/Documents/rotordin_frontend/artifacts/primsf-continuation-20261118`;
seu resultado foi dez manifestos, cinco pares e dois refinamentos espaciais
conferidos, sem erros. Ele não executa nem modifica a campanha.

## Limite temporal nesta conferência

Às 20:34:31 UTC, ainda não havia diretório/plano/manifestos
`comparison_final_20261118_dt`. O orquestrador local
`launch_corrected_comparison.py` agenda D_OSC C3/C5, 51+51 células,
dt máximo 0,0625 s, depois de concluir a coorte base. Isso separa os controles
planejados de h e dt, mas não comprova que a etapa temporal foi executada.

Quando seus arquivos existirem, um adendo conferirá somente a etapa temporal
contra o freeze e o par base; não repetirá a auditoria dos dez casos acima.

## Adendo — etapa temporal e extensão disponíveis

Conferência concluída em **18/11/2026, 20:39:22 UTC**, limitada aos três novos
casos e aos seus pares base já auditados; nenhuma métrica dinâmica foi
reanalisada. A pendência temporal registrada acima fica resolvida quanto às
entradas e à proveniência, sem implicar conclusão das execuções.

Foram lidos os planos, manifestos, decks e `run_metadata.txt` de:

- `comparison_final_20261118_dt_D_OSC_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N51_51_dt0.0625`;
- `comparison_final_20261118_dt_D_OSC_C5_BAROTROPIC_MATCHED_LAGGED_ROE_N51_51_dt0.0625`;
- `comparison_final_extended_20261118_D_OSC_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N51_51_dt0.125`.

Todos conservam o freeze `490936a`, as 23 fontes e o binário `290280e3…b1bb7`
identificados integralmente acima. Os hashes de namelist conferem; a referência
física, EOS, vazões, geometria, propriedades, escalas e perturbação de 1%
coincidem com os respectivos pares base. Os dois casos temporais diferem entre
si somente em fechamento e diretório de saída.

| Comparação com o par base | Únicas diferenças nos decks e metadados efetivos |
|---|---|
| C3 e C5: N51, dt máximo 0,125→0,0625 s | `dt_max_s`, `output_interval_s`, `out_dir` |
| C3: N51, dt máximo 0,125 s, duração 2400→3600 s | `time_dim`, `out_dir` |

A cadência de observação temporal passa de 0,1 para 0,0625 s conforme o runner;
isso deve permanecer visível na análise das séries. A malha não muda no refino
de passo. A extensão altera somente a duração solicitada e o destino dos dados.

| Plano | Início UTC | SHA-256 |
|---|---|---|
| `validation/comparisons/comparison_final_20261118_dt/plan.json` | 20:37:06.626284 | `4506fe1446b1489853dc0b402a3373d0092f1f2017b53cdaa299b970be056b0b` |
| `validation/comparisons/comparison_final_extended_20261118/plan.json` | 20:36:48.314686 | `db76b14620d4a7b03bed18c71996f58d9e7fed2dbb98bd12e07ac05009d5f59f` |

Na leitura, temporal C3 e extensão C3 constavam `RUNNING`; temporal C5 constava
`FAILED`. Esses estados são um retrato do horário da conferência, não uma
classificação física do regime. O helper local `audit_final_dt_readonly.py`, no
mesmo diretório de artefatos do helper anterior, encerrou sem erros. Nenhum
solver/teste foi reexecutado e nenhuma fonte congelada foi alterada.

### Verificação do efeito da cadência de saída

No código congelado `490936a`, `output_interval_s` é copiado de `main.f90`
para `cfg` e usado somente para reagendar `next_output` em
[`module_simulate_pipe_riser.f90`](../../../scr/module_simulate_pipe_riser.f90),
linha 729. A condição da linha 655 escreve o monitor **depois** de aceitar `U`,
avançar o tempo e acumular os fluxos integrados. A seleção/redução/crescimento
do passo nas linhas 585, 643 e 663–666 não usa esse intervalo nem `next_output`.
As consultas `point_state` e `face_pair` do monitor recebem o simulador e os
estados com `intent(in)`; o inventário recalculado ali alimenta apenas a saída.
Assim, essa opção muda a gravação das amostras, sem mudar o passo numérico ou
o estado conservativo. O solver também não encurta passos para coincidir com
os instantes de saída: grava no primeiro estado aceito após atingir o limiar.
A diferença de amostragem permanece registrada para a análise das séries.
