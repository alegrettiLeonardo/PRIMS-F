# Diagnósticos P0/P1 — domínio conservativo único

Os programas deste diretório são ferramentas isoladas. Não ativam MUSCL/TVD,
SDIRK2, Rusanov ou outro fechamento no transiente padrão. As evidências estão
identificadas pelos operadores `1365094`, `18aa7f7` e `994736f`.
O relatório de decisão é [QUALIFICATION_P0P1](../../../docs/QUALIFICATION_P0P1.md).

## Compilar e reproduzir

GNU Fortran 13.3.0; Fortran 2018; release `-O2`; BLAS/LAPACK 3.12.0.
Não reutilizar objetos de outro operador. Os dois manifestos `frozen_*_BUILD.json`
registram hashes dos fontes/binários; o primeiro também registra o Makefile.
Executáveis e PDFs ficam fora do Git. Use diretórios de saída novos.

```sh
make -j2 FC=gfortran LIBS='-llapack -lblas'
make check_all FC=gfortran LIBS='-llapack -lblas'
make test FC=gfortran LIBS='-llapack -lblas'
PYTHONPATH=scripts python3 -m unittest discover -s tests -p 'test_*.py' -v
python3 tests/test_boundary_telemetry_cli.py --binary bin/prims_f_check_all --output tmp/boundary_novo
python3 tests/test_feed_contract_cli.py --binary bin/prims_f --output tmp/feed_novo
```

`run_p0_campaign.py --frozen DIR_BINARIOS --output SAIDA_NOVA` executa os
casos declarados no protocolo, com limites de custo explícitos. `DIR_BINARIOS`
deve conter `prims_f` e seu `build_manifest.json`; nunca inferir o operador
pelo nome do arquivo. As entradas efetivas são gravadas por execução.

Para compilar os drivers contra os objetos release do operador escolhido:

```sh
objs=$(find tmp/release -name '*.o' ! -name main.o -print)
gfortran -std=f2018 -O2 -ffree-line-length-none -fopenmp -Itmp/release \
  validation/diagnostics/p0p1_20261119/restart_refinement.f90 $objs \
  -llapack -lblas -o /tmp/primsf_refine
# criar SAIDA_NOVA antes da execução; U_INICIAL contém uma linha U1 U2 U3 por volume
/tmp/primsf_refine SAIDA_NOVA U_INICIAL N_POR_TRECHO DT_MAX_S DURACAO_S
```

O driver fixa explicitamente as propriedades do caso OSC BAROTROPIC e o mesmo
solver. Para a campanha registrada, foi ligado aos objetos imutáveis `18aa7f7`.
`run_restart_refinement.py --initial-fields CAMPANHA_A/fields.csv --output NOVA`
reproduz a projeção conservativa e as sete combinações do protocolo. O SHA
registrado por uma reprodução será o do seu checkout, não automaticamente
`18aa7f7`. A condição inicial comum original está em `common_initial_U.txt`.
Ela é reconstruída dos campos de `1365094` em 2400 s: não é uma continuação
bit a bit nem uma nova integração de 2400 s em todas as malhas.

Extensões usam o U aceito do checkpoint, com o mesmo executável e controles;
o driver reinicia seu dt_start em 0,001 s. O offset físico e essa reinicialização
são registrados. Não concatenar os acumuladores de massa sem contabilizar
os inventários de cada segmento. O exportador de checkpoint contém U,
geometria, escalas e passo de trabalho; não constitui um restart geral de
Chexal com múltiplos ramos/histórico.

```sh
python3 validation/diagnostics/p0p1_20261119/analyze_refinement.py \
  --runs RESULTADOS_BASE --extensions RESULTADOS_EXTENSAO --output ANALISE_NOVA
```

`reduced_wave.f90` exercita ondas pequenas e reflexão em pressão fixa, para
líquido e gás puros. Declare a fase do reservatório explicitamente; o domínio
atual não fornece sua composição. `replay_p0.f90` é restrito à malha 51+51
dos snapshots e depende de `qualification2_20261119/globalization_flux/offline_api.f90`.
Seus argumentos são pasta, Cg SI, base líquida, método de contorno e fluxo
(1 Roe/2 Rusanov). Recompila-se com os módulos debug e o offline_api antes
do programa. `probe_parameters.dat` contém omega_u, Ps, jl, jg, Gl, Gg, reservado;
Gl/Gg são fluxos mássicos por área, não vazões totais.

## Protótipo linear P1

`wave_api.f90` é ligado como biblioteca PIC com todos os módulos do mesmo
operador e `offline_api.f90`. `linear_waves.py --library waves.so --fixture
PASTA_SS1 --output NOVA` usa o Jacobiano real de um estado admissível, seus
três autovalores reais e a recuperação não linear de primitivas para qualificar
os traços reconstruídos. Ondas características periódicas são avançadas
separadamente: isto não comprova estabilidade do sistema não linear acoplado.

A comparação espacial usa FV1 ou MUSCL/MC, ambos com SSPRK2. MC é escrito
explicitamente como minmod(2 dl, (dl+dr)/2, 2 dr). Se um traço sai do domínio,
reduz-se a inclinação conservando a média; se a média for inadmissível, a
proposta é recusada. A comparação temporal mantém o símbolo FV1 e verifica
BE ou SDIRK2 contra a solução semidiscreta exata. O tableau adotado é
`A=[[gamma,0],[1-gamma,gamma]], b=[1-gamma,gamma]`,
`gamma=1-1/sqrt(2)`. Por derivação própria, `b^T 1=1`,
`b^T c=2 gamma-gamma^2=1/2` e
`R(z)=[1+(1-2 gamma)z]/(1-gamma z)^2`. Os passos do protótipo linear são
adimensionais; velocidades físicas e escalas são fornecidas separadamente.
Não se emprega gamma na EOS isotérmica.

As famílias de métodos são documentadas por van Leer (1979),
[DOI](https://doi.org/10.1016/0021-9991(79)90145-1), e Alexander (1977),
[DOI](https://doi.org/10.1137/0714068). O texto integral editorial de van Leer
retornou HTTP403, e o de Alexander não estava aberto; não atribuímos páginas
ou equações específicas dessas fontes às fórmulas próprias acima.

## Escopo das evidências

`P0_EVIDENCE.json` contém os seis balanços, os erros de imposição física e a
comparação curta `18aa7f7`/`994736f`. Normalização por massa processada é
informada separadamente daquela pelo inventário inicial. Fase exatamente
ausente usa erro absoluto, com normalização indefinida; nenhum piso de massa.
`EOS_CAMPAIGN_1365094.json` preserva a campanha EOS corretamente parametrizada.
Os pilotos anteriores com dt_max=0,2 s continuam fora do Git, identificados
como protocolo inválido para a comparação pretendida de 0,125 s.

Os PNGs são novos gráficos científicos, não cópias dos PDFs. A sobreposição
experimental usa a curva digitalizada da Fig.18 de Baliño et al., com origens
temporais independentes; a banda é incerteza gráfica, não intervalo estatístico.
Os limiares de frentes são cruzamentos geométricos de alfa, testados pelas
regressões existentes; não são prova automática de interface material,
penetração, blowdown, SS1/SS3/OSC ou severe slugging validado.
