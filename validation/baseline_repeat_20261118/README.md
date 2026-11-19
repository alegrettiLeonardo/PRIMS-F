# Fechamento da campanha congelada — 18/11/2026

As duas repetições antes interrompidas completaram **2400 s**, sem rejeição de
passos. A baseline permanece em `05ad799af3b82aa540cb6c87b58debe914e52681`,
**sem alpha-tail e sem Chexal**. Ambas reiniciaram em t=0; nenhuma série foi
concatenada à tentativa interrompida. N201 aceitou 19222 passos; N51, 76815.
Tempos de relógio: 2101,696 s e 1942,617 s, respectivamente, executando em paralelo.

## Dez intervalos finais e refinamentos independentes

Onze picos consecutivos encerrados entre mínimos delimitam dez intervalos.
O período e a semi-amplitude são os observáveis do protocolo de envelopes;
não constituem classificação experimental de slugging.

| Malha | dt máximo (s) | Período médio (s) | Semi-amplitude média (Pa) | Variação entre primeiro/último pico | Proveniência |
|---|---:|---:|---:|---:|---|
| 51 + 51 | 0.125 | 23.103011 | 9699.735 | +2.398601% | arquivo histórico |
| 51 + 51 | 0.0625 | 23.265630 | 13763.806 | +0.000523% | arquivo histórico |
| 51 + 51 | 0.03125 | 23.344185 | 15347.394 | -0.000031% | nova repetição |
| 101 + 101 | 0.125 | 23.266396 | 12769.678 | +0.021283% | arquivo histórico |
| 201 + 201 | 0.125 | 23.345470 | 14011.547 | +0.002630% | nova repetição |

Nas novas repetições, pressão, frações líquidas de base/saída, vazões mássicas
de ambas as fases na base/saída e inventários de ambas as fases no conjunto e
nos trechos foram avaliados em cada intervalo: **13 canais**. A maior diferença
relativa das excursões médias entre os primeiros cinco e os últimos cinco
intervalos é 0,002773% em N201 e 0,000640% em N51. Isso evidencia uma janela
numérica acomodada nessas grandezas, sem afirmar ciclo limite assintótico.

O refinamento espacial 101→201, a dt máximo 0,125 s, aumenta a semi-amplitude
final em aproximadamente 9,73%. O refinamento temporal 0,0625→0,03125 s,
mantendo N51, aumenta-a em aproximadamente 11,51%. Portanto a amplitude ainda
não está demonstrada convergente. As duas comparações mais finas mudam o
período em aproximadamente 0,34%; isso não prova convergência do campo ou do
regime físico. `last_ten_refinements.csv` guarda as diferenças exatas. A janela
comum 900–2400 s ainda inclui acomodação em algumas malhas; seus alertas são
preservados separadamente, sem confundir essa janela com os dez ciclos finais.

## Conservação e interpretação das fases

| Repetição | Máx. erro líquido (kg) | Máx. erro gás (kg) | Máx. erro relativo líquido/gás |
|---|---:|---:|---:|
| N201, dt 0,125 | 9,044925e-10 | 9,227064e-11 | 3,1873e-13 / 1,2697e-11 |
| N51, dt 0,03125 | 2,725528e-9 | 2,370512e-9 | 9,6043e-13 / 3,2619e-10 |

Os erros são os balanços integrais gravados pelo solver em todo o transiente,
normalizados pelo maior valor entre inventário inicial e massa acumulada de
entrada. Todos passam o critério numérico 1e-6 relativo; isso não valida a física.

Em 900–2400 s, N201 apresenta 65 episódios de acumulação/compressão, 62
excursões de produção líquida e 65 candidatos a blowdown; N51 apresenta
65/64/65. São proxies dos sinais nomeados e dos limiares declarados no JSON,
não eventos experimentais observados. A transferência de gás na base continua
positiva (mínimos 0,00136519 e 0,00124482 kg/s): não há reinício de gás contado
após bloqueio, e **bloqueio completo não é requisito universal de OSC/SS3**.
Penetração e compressão permanecem interpretações indiretas; posição de frente
e temporização experimental não foram medidas. A sequência e as amplitudes
não estabelecem validação experimental.

## Proveniência, testes e arquivos

Os cinco casos da coorte têm **todos os hashes Fortran idênticos**, conferidos
também contra o worktree congelado. EOS barotrópica, propriedades, geometria,
vazões, perturbação, fechamento ID3 e opções MATCHED/LAGGED/ROE foram mantidos.
As séries históricas foram recuperadas seletivamente do ZIP fornecido, com
hashes de arquivo e de membros em `historical_extraction.json`; não houve
sobrescrita. O compilador histórico era GNU Fortran 14.2.0; estas repetições
usam 13.3.0. Os binários diferem e as comparações preservam esse alerta: não
representam uma reprodução bit a bit no mesmo compilador.

Executados na baseline: `make test FC=gfortran LIBS='-llapack -lblas'`
(2185 verificações, seis programas) e os 16 testes Python da CI, todos PASS.
`test_evidence.json` mantém versões, hashes e contagens. O pós-processamento
usou NumPy 1.26.4 com integração trapezoidal equivalente à API NumPy 2; a
adaptação foi verificada numa série histórica completa antes destas análises.

`summary.json` contém métricas compactas, fontes, hashes, versões, fases,
limiares e as médias/excursões dos 13 canais em cada um dos dez intervalos.
`manifests/` e `cases/` preservam as entradas e os términos das novas execuções.
Séries, campos, logs e análises completas permanecem fora do Git em:

```
/home/leonardoab/Documentos/PRIMS-F-baseline-20261118/validation/results/
/home/leonardoab/Documentos/PRIMS-F-baseline-20261118/validation/historical_validation_20261118/
/home/leonardoab/Documentos/PRIMS-F-baseline-20261118/validation/baseline_repeat_20261118/final_analysis/
```

O arquivo local `PRIMSF_Baseline_Repeat_20261118.zip` reúne as cinco séries da
comparação, os dois novos campos e as análises finais; `artifact_manifest.json`
identifica caminho e hash. Não contém PDFs bibliográficos nem executáveis.

## Reprodução

Num worktree separado em 05ad799, compile com `make -j2 FC=gfortran
LIBS='-llapack -lblas'`. Execute os comandos abaixo a partir dessa raiz, usando
um rótulo novo para preservar a tentativa atual:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 scripts/run_validation.py --case D_OSC --cells 201 --dt 0.125 --time 2400 --closure 3 --eos BAROTROPIC --junction MATCHED --dissipation LAGGED --flux ROE --timeout 86400 --label NOVO_ROTULO
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 scripts/run_validation.py --case D_OSC --cells 51 --dt 0.03125 --time 2400 --closure 3 --eos BAROTROPIC --junction MATCHED --dissipation LAGGED --flux ROE --timeout 86400 --label NOVO_ROTULO
```

`analysis_commands.json` registra os comandos realmente executados. A seleção
do analisador de envelopes deve ser explícita por `--run-dir` ou `--pattern`,
agrupando apenas fontes e física compatíveis; as opções de seleção estão na
branch de continuação, não no solver congelado. `summarize.py --repo CAMINHO_DA_BASELINE
--analysis PASTA_DAS_ANALISES` reconstrói o resumo completo e recusa deriva de
fontes. O limite de relógio não altera os 2400 s físicos solicitados.
