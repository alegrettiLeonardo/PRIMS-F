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

A conversão STANDARD da alimentação e a EOS local são operações distintas.
No freeze 69813ee, main.f90/315 usa a densidade ideal de referência, enquanto
module_find.f90/121 usa rho_g=p/Cg². O [cálculo independente Decimal60](../validation/diagnostics/static_bc_primary_20261119/standard_vs_barotropic.json)
confirma rho_std_ideal=1,20403094823 kg/m³ e mug=1,39667589995 kg/m²/s.
A EOS BAROTROPIC com Cg=343 m/s daria rho_g=0,861035792909 kg/m³ em Pstd;
a razão de densidades é 0,715127625394. Portanto a mesma massa corresponde
a jg_model(Pstd)=1,62208808443 m/s (razão 1,39835179692), embora a prescrição
experimental STANDARD seja 1,16 m/s. Isso limita a equivalência dos observáveis
com a fonte ideal e não autoriza reescalar a alimentação. Os seis casos e as
extensões preservam EOS e massa prescrita. ISOTHERMAL já é uma opção do
programa, mas não participa deste freeze/refinamento.

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
sobrescrever 3600 s. A nova coorte de 7200 s foi autorizada com até dois
workers externos, um thread por solver e teto global de quatro. Os seis
casos básicos permaneceram sequenciais; max_workers=1 do deck e os demais
controles numéricos continuam idênticos. O novo plano declarou o teto antes
da selagem. Ausência de ciclos não prova estabilidade física e não autoriza
ajustar parâmetros para produzir oscilação.

Os [seis resultados básicos de 3600 s](../validation/diagnostics/static_bc_3600_20261119/README.md)
foram concluídos e auditados. O pacote conserva o resultado provisório Ps+
e a ausência de ciclos resolvidos jg+, além das quatro comparações periódicas
que satisfizeram os critérios numéricos. As extensões de 7200 s têm plano e
selagem separados; o pacote de 3600 s não será substituído.

### Compatibilidade do observável alpha_p e da pressão

Baliño2010/p.652 chama alpha_p exp de fração de vazio experimental no
oleoduto e informa médias, salvo entradas marcadas Max. O texto não define
janela nem ponderação espacial/volumétrica. O densitômetro gama é descrito
na base do riser; a Fig.9 mostra G1 no trecho final do oleoduto junto a P4,
antes da curva da base, sem coordenada métrica ou extensão do feixe.
A associação de G1 à coluna alpha_p é uma inferência, não um mapeamento
instrumental explicitamente fornecido. Na p.655, alpha_p sim é identificado
como o valor do estado estacionário; a p.645 assume alpha_p constante no
transiente. Assim, alpha_p exp=0,75 do caso13 não valida diretamente uma
média espaço-tempo de todo o oleoduto no modelo distribuído. Pressão de base,
fração local e médias espaciais devem permanecer observáveis distintos.

A Fig.18 identifica pressão na base do riser, mas não associa explicitamente
essa série a P4 ou P5. A Fig.9 mostra P4 no trecho final do oleoduto e P5 na
curva/base do riser, sem fornecer coordenada métrica ou extensão da tomada.
O traço de pressão na junção é um proxy numérico fixo; a pressão da primeira
célula é outro observável. Nenhum dos dois recebe equivalência instrumental
exata apenas a partir desse desenho.

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
