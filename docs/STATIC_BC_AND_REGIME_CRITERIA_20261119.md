# Sensibilidades estáticas OSC e critérios de regime

## P1: seis contornos constantes, uma variável por vez

[Baliño2010](https://doi.org/10.1016/j.ijmultiphaseflow.2010.04.003), p.652,
informa que médias e desvios foram calculados dos sinais experimentais.
Não define nesse trecho o estimador do desvio: **não** atribuímos aos sinais
± da Tabela3/caso13 um desvio-padrão, erro instrumental, nível de confiança
ou distribuição probabilística. São escalas para sensibilidades estáticas.
Nenhuma frequência, fase, correlação entre contornos ou excitação temporal
foi inferida. Fonte e SHA256 estão na [auditoria P4/P5](P4_P5_PRIMARY_EVIDENCE_20261119.md).

Referência: OSC 3600 s, 51+51 células, dtmax=0,125 s, freeze
`69813ee70551fa4a3f075919aeafc45de9693f11`, binário
`7442291fd6cdc311627e79500dc5a4804813b42c7bf4daecf88b4b8b2322288a`.
Manifesto nominal: `PRIMS-F-numerics-results-20261119/v1_osc3600/runs/D_OSC/manifest.json`.
Tgas, EOS, geometria, fluidos, fechamento, inicialização e controles numéricos
permanecem idênticos, inclusive `tol_res=tol_inc=1e-8` e perturbação inicial
0,01. Cada condição mantém seu contorno constante durante toda a simulação.

| Variante | Ps (Pa absoluto) | jg0 padrão (m/s) | jl (m/s) | mul (kg/m²/s) | mug (kg/m²/s) |
|---|---:|---:|---:|---:|---:|
| Ps− | 201000 | 1,16 | 0,547 | 545,38635 | 1,396675899947762 |
| Ps+ | 211000 | 1,16 | 0,547 | 545,38635 | 1,396675899947762 |
| jg− | 206000 | 1,10 | 0,547 | 545,38635 | 1,3244340430539123 |
| jg+ | 206000 | 1,22 | 0,547 | 545,38635 | 1,4689177568416119 |
| jl− | 206000 | 1,16 | 0,544 | 542,3952 | 1,396675899947762 |
| jl+ | 206000 | 1,16 | 0,550 | 548,3775 | 1,396675899947762 |

Os overrides de massa são explícitos:
`mul=rho_l0*jl` e `mug=P_standard*jg/(Rg*T_standard)`, com Pstandard=101300 Pa
e Tstandard=293,15 K. Variar Ps não muda mug porque a base do jg é padrão.
Alterar jg/jl também altera apenas seu override correspondente; não se
mantém um fluxo de massa antigo incompatível com a vazão pretendida.

[Gerador](../scripts/prepare_static_bc_sensitivity.py) e
[quatro testes](../tests/test_static_bc_sensitivity.py) executados: **PASS**.
O gerador confere hashes, controles, metadados efetivos e conjunto exato de
campos alterados; não inicia solver. Plano/decks/selagem ficam externamente
em `/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119/bc_static/`.
O executor imutável `ea973346fe0c9bec8a96df88abc98bd205d54b708668b717a2e47de4857e7a24`
foi lançado com um worker, em sequência, respeitando quatro solvers globais.

O [watcher](../scripts/watch_static_bc_sensitivity.py) usa as mesmas cópias
congeladas dos analisadores de ciclos e conservação, lê somente execuções
terminadas e não lê fields.csv. O critério físico de balanço vem dos passos
aceitos P2, para seis combinações fase/trecho, sem o piso absoluto legado do
analisador de ciclos. O resumo distingue término, disponibilidade de dez
intervalos e acomodação. Resultados não acomodados permanecem identificados;
eventual extensão deve mudar somente a duração em uma nova coorte, sem
sobrescrever 3600 s. Ausência de ciclos não prova estabilidade física e não
autoriza ajustar parâmetros para produzir oscilação.

## P2: diferença de densidade e perímetro na fonte primária

Foi relida visualmente a Eq.1/p.6 (PDF10) do relatório dos próprios autores,
[Taitel, Lee e Dukler (1977), DOE/ID/01571--T1](https://www.osti.gov/servlets/purl/5902200).
Seu SHA256 é `c953722d3abb7e22369b4c89b6fd15c6a28482d896a3bbdfc8f654f41f9854fa`.
A fonte resume a teoria estacionária de 1976, mas **não** alegamos leitura
integral do artigo de 1976: o endpoint editorial forneceu somente resumo.
Para tubo horizontal e candidato estratificado, a Eq.1 imprime:

```text
Ucrit_TLD77² = (1-h/D)² * g * A_g/Si * (rho_l/rho_g - 1)
Si = dA_l/dh = D*sin(phi)                 [Eq.11/p.9]
```

Burr15/p.3 usa `Sl=D*phi` e a razão `(rho_l-rho_g)/rho_l`. As Eqs.53–57 de
Burr mudam a densidade, mas ainda não se tornam a Eq.1 de 1977 por manterem
Sl. São diferenças de modelo, não conversões de unidades. Não implementamos
troca de fórmula ou fechamento. A fonte separa perda de estratificação,
intermitente/anular e suave/ondulado; uma desigualdade não identifica sozinha
golfada. A hipótese de gás muito mais rápido que líquido é explícita na p.9.

O [calculador independente](../scripts/compare_stratified_primary_criteria.py)
reutiliza exatamente os envelopes experimentais da Fig.9 de Carraretto2020
já registrados. Não encontra densidades por ponto nessa fonte; por isso
expressa a resposta de 1977 como condição sobre a razão de densidades:

| Controle | Burr15: margem mínima ug−Ucrit (m/s) | Eq.1/1977 não aponta perda em todo envelope se rho_l/rho_g > |
|---|---:|---:|
| Fig.9(a), terceiro marcador | +0,346475 | 9,974269 |
| Fig.9(b), terceiro marcador | +0,446875 | 12,523962 |

Para alpha>0,5 desses envelopes, `(1-h/D)²*A_g/Si` cresce com alpha.
O limite usa esse fator no alpha mínimo e ug no máximo JG/alpha mínimo.
As razões mínimas ug/ul no envelope são 4,9110 e 4,7529; reportá-las não
valida a hipótese aproximada de gás muito mais rápido que o líquido. Os valores de
rho_l/rho_g acima são **limiares condicionais derivados**, não densidades
experimentais atribuídas. A discordância de Burr15 continua robusta; a
alternativa não recebe uma aprovação experimental incondicional.
[Três testes geométricos/de limite](../tests/test_stratified_primary_criteria.py)
passaram, e o [JSON](../validation/diagnostics/static_bc_primary_20261119/primary_criteria.json)
registra as equações, versões e lacunas.

Baliño2010/caso14 é não estratificado SF, mas alpha médio=0,74 e jg0 padrão
não estabelecem um estado local com interface plana. A comparação registra
`NOT_APPLICABLE_GEOMETRY_NOT_ESTABLISHED`, sem calcular altura ou limiar.
O par local antes/depois de perda de estratificação com holdup e JG continua
**BLOCKED dataset**. `AUTO_BURR15` permanece diagnóstico, restrito ao candidato
estratificado; a guarda com fechamento drift e a checagem `cid==1` permanecem
como documentadas na auditoria anterior. Nenhuma lei de retorno foi criada.
