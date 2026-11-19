# Continuação da qualificação numérica e física

Base de desenvolvimento `12eaa07`, operador de referência `69813ee`.
Os resultados desta campanha usam o binário release já verificado,
SHA256 `7442291fd6cdc311627e79500dc5a4804813b42c7bf4daecf88b4b8b2322288a`.
O compilador é GNU Fortran13.3.0, Fortran2018, `-O2` e LAPACK/BLAS.
A raiz externa é `PRIMS-F-qualification2-results-20261119`.

## Campanhas independentes

O estudo espacial mantém teto de passo0,125s e usa 51,101,201 e401 volumes
por trecho. Cada trecho tem volumes uniformes em comprimento de arco,
oleoduto57,4m e riser12,546254m. O estudo temporal usa51 volumes por trecho
e tetos0,125,0,0625 e0,03125s. As listas completas dos passos aceitos estão
em `steps.csv`. Os casos começam da mesma definição física de estado
estacionário, integrada na malha própria de cada caso.

O nominal51/0,125 em3600s é a execução congelada da campanha numérica
anterior, com o mesmo operador/binário/controles. Novas execuções ficam
identificadas por manifesto. Não se misturam dados de outro SHA ou de outro
controlador. As sensibilidades estáticas usam o mesmo nominal e variam
uma entrada por vez dentro dos desvios publicados do caso13. As mudanças
de vazão atualizam somente a massa da fase correspondente.

| Eixo de sensibilidade | Nominal | Extremos |
|---|---:|---:|
| Pressão absoluta de separador (Pa) |206000|201000/211000|
| Velocidade superficial padrão de gás (m/s) |1,16|1,10/1,22|
| Velocidade superficial líquida de referência (m/s) |0,547|0,544/0,550|

Esses desvios são os publicados a partir dos sinais. Não se atribui uma
distribuição estatística ou frequência de oscilação. São ensaios de
sensibilidade estática ao ponto de operação. O efeito de forçamento
externo temporal permanece sem identificação quando faltam os sinais
sincronizados da alimentação e do separador.

## Critérios fixados

1. Um caso só termina normalmente se o processo retorna0, status0 e o tempo
   físico solicitado foi atingido. Falhas e interrupções ficam identificadas.
2. Iniciar com3600s. Onze picos de pressão cercados por mínimos delimitam
   dez intervalos completos. Excursão mínima100Pa e RMS mínimo1Pa evitam
   inventar ciclos em sinais planos.
3. Acomodação exige até1% de variação da semi-amplitude entre os extremos,
   e até1% nas médias/excursões não nulas dos13canais entre as primeiras
   e últimas cinco janelas. Extender os casos que não passam.
4. Convergência compara pares consecutivos por eixo nas janelas acomodadas.
   O alvo é1% em período e semi-amplitude. Reportar também médias,
   formas de onda, inventários, fluxos, defasagens e tendência de refinamento.
   Uma sequência de amplitudes crescente não determina ordem assintótica.
5. As formas de onda usam256fases por ciclo com origem no pico de pressão.
   A origem é comum aos canais. Não se otimiza uma translação por canal
   para reduzir o erro. Registrar erroRMS absoluto, relativo à excursão
   de referência e a fase do máximo dos canais.
6. Conservar o critério1e-6 por inventário inicial para cada fase e trecho.
   Registrar também erros absolutos, massa absoluta processada e por ciclo,
   a diferença dos dois fluxos da junção e o residual aceito acumulado.
7. A auditoria de contorno reconstrói os estados a partir de cada quadro
   aceito armazenado. Confere espectros reais, número de ondas entrantes,
   velocidades assinadas e posto das restrições projetadas nos modos.
   Os estados Newton recusados têm um replay separado.

## Comparação experimental

A pressão na junção e na primeira célula do riser são reportadas
separadamente. A referência é o caso13, Fig.18 e Tabela3 de Baliño,
Burr e Nemoto2010. A digitalização cobre1200–1500s, em544colunas,
com pressão absoluta e envelope de espessura/calibração.

Comparar excursão pico a pico da mesma duração e com a mesma definição.
O período publicado19,3–33,8s tem precedência sobre uma estimativa feita
do pequeno recorte raster. O espectro próprio do recorte precisa informar
resolução de frequência, janela e incerteza gráfica. Reportar a razão
entre os dois picos e os possíveis harmônicos, sem assumir que todo pico
seja uma oscilação independente. A comparação temporal usa origens
independentes e não pressupõe sincronização.

SS1/caso2 tem85,6s e SS3/caso11 tem76s. A presença de gás na base é
considerada junto às frentes, holdup, inventários e fluxos. A Eq.15
continua diagnóstica, com a divergência de densidades das Eqs.15/54
visível. Critérios alternativos exigem teste independente e regra
reversa antes de entrar no operador.
