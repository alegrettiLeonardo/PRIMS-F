# Evidencia primaria P4/P5 - 2026-11-19

Fontes e limites: [nota](../../../docs/P4_P5_PRIMARY_EVIDENCE_20261119.md).

Executado localmente, sem solver:

    python3 scripts/assess_burr15_primary_bounds.py --output validation/diagnostics/primary_evidence_p4p5_20261119/burr15_primary_bounds.json
    python3 -m unittest discover -s tests -p test_burr15_primary_bounds.py -v

Resultado: 4 testes PASS (0.003 s). O JSON fixa os retangulos lidos da figura e os limites derivados. Nenhum PDF ou imagem da fonte esta no Git.

P5: scripts/digitize_balino2010_case13.py foi executado com o PDF local cujo SHA256 esta fixado no script. Foram obtidas 544 colunas, zero ausencias. CSV/JSON/figuras proprias estao fora do Git em /home/leonardoab/Documentos/PRIMS-F-numerics-results-20261119/primary_evidence/balino2010_case13/. A figura experimental foi inspecionada visualmente antes e depois da digitalizacao. Nao e recuperacao de dados brutos nem comparacao sincronizada de fase.
