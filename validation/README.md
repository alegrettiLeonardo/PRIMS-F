# Evidências e reprodução

`docs/VALIDATION_RESULTS.md` resume a campanha. O primeiro ponto de entrada para
continuar o trabalho é `docs/CODEX_HANDOFF_STATE.md`; os testes de operadores
ficam em `tests/` e os comandos de compilação no README da raiz.

O Git distribui os fontes, decks, testes, resumos e diagnósticos compactos.
As séries completas, campos, gráficos, logs e demais dados da campanha estão no
arquivo `PRIMSF_Evidencias_Implementacao.zip` entregue com o relatório. Para
reanalisar os dados antigos, extraia a árvore `validation/` desse arquivo na
raiz de uma cópia da branch de entrega, preservando resultados locais com os
mesmos nomes. O manifesto do ZIP lista os arquivos distribuídos separadamente
pelo Git e seus hashes. Os PDFs de terceiros e os executáveis não são incluídos.

Os manifests registram os hashes das fontes, do executável e do deck utilizados.
Caminhos absolutos históricos são proveniência do ambiente original; não são
pré-requisitos para compilar no WSL. A referência bibliográfica vigente está em
`docs/benchmark_sources.json`; os snapshots `evidence/benchmark_sources_at_run.json`
e `evidence/benchmark_sources_at_refined_run.json` preservam as versões registradas
antes das atualizações bibliográficas posteriores.

`runtime_limited/` preserva duas tentativas encerradas pelo limite de relógio.
As duas repetições correspondentes perderam a sessão de execução e ficaram
`INTERRUPTED`, com evidência em `runtime_interruption.json`. Nenhuma delas tem um
`run_status.txt` de término escrito pelo solver. A ausência desse arquivo não
foi preenchida com um status de sucesso ou falha física inventado.

Resultados de verificações locais, convergência de malha/passo, conservação de
massa e concordância experimental são evidências diferentes. Os nomes `D_SS1`
e `D_SS3` identificam condições de referência, não classificações impostas à
solução calculada.
