# Sensibilidade independente de tol_res:1e-8→1e-10

O operador permaneceu `a179e51d9d63e8aeed91a94d829b5b483d107864`,
binário SHA256 `02ce51462ec3c954cec9364b17444b7279be72ec150f87e307f1b5fc8dd0a4c3`.
A comparação textual dos decks verificou que somente `tol_res` mudou de
1e-8 para 1e-10, além da pasta de saída. `tol_inc=1e-8`, física, fluxos,
discretização, controlador temporal e globalização são idênticos. Não
houve recompilação. O timeout externo de 900s wall foi declarado antes dos
dois smokes. Esse aperto testa erro de aceitação do Newton, não ajusta
coeficientes físicos nem relaxa o domínio.

| Caso | Término observado | CPU(s) | Wall(s) | Último steps.csv | Último conservation.csv |
|---|---|---:|---:|---|---|
|D_SS1|EXTERNAL_TIMEOUT|899.262720|900.001755|step2138,t=78.577569191712726s|step2138,t=78.577569191712726s|
|D_SS3|NEWTON_GLOBALIZATION_FAILURE|47.377847|47.472056|step810,t=97.131852191077598s|step810,t=97.131852191077598s|

Em TIMEOUT, cada buffer fornece somente seu último registro completo,
podendo terminar em passos diferentes. A contagem é um limite inferior,
não um estado terminal reconstruído. O snapshot remanescente é uma
rejeição anterior. Em falha interna, o motivo e a proposta estão nos
sidecars preservados; não se interpreta NaN/huge como norma calculada.

## Máximos de conservação no intervalo observado

Para cada combinação, o máximo do módulo do balanço bruto em todo o
ledger é dividido separadamente por M0 e pela massa absoluta processada
pelas duas faces até o fim da janela. Ambos devem ser <=1e-6. O resíduo
Newton acumulado é registrado, sem ser subtraído para corrigir o erro.
São janelas iniciadas em t0, não ciclos acomodados. Denominadores, máximos,
inventários e SHA de cada ledger estão em `balance_normalizations.json`.

D_SS1, até78.577569191712726s:

| Fase/região | Máx.absoluto(kg) | Máx./M0 | Máx./Mprocessada | Ambos<=1e-6 |
|---|---:|---:|---:|---|
|l_system|2.215077e-07|2.764134e-09|1.201607e-09|PASS|
|l_pipe|4.482969e-07|8.207018e-09|2.489547e-09|PASS|
|l_riser|2.268206e-07|8.890474e-09|1.230278e-09|PASS|
|g_system|1.260328e-09|7.036216e-09|2.670837e-08|PASS|
|g_pipe|1.358497e-09|7.736752e-09|2.433312e-08|PASS|
|g_riser|9.816955e-11|2.781019e-08|1.695753e-09|PASS|

D_SS3, até97.131852191077598s:

| Fase/região | Máx.absoluto(kg) | Máx./M0 | Máx./Mprocessada | Ambos<=1e-6 |
|---|---:|---:|---:|---|
|l_system|9.237056e-14|1.734282e-15|7.325752e-16|PASS|
|l_pipe|1.961098e-12|5.720080e-14|1.667645e-14|PASS|
|l_riser|1.918465e-12|1.010936e-13|1.521310e-14|PASS|
|g_system|4.718448e-16|2.137307e-15|1.857257e-15|PASS|
|g_pipe|4.825029e-13|2.362346e-12|1.795508e-12|PASS|
|g_riser|4.823086e-13|2.919770e-11|1.726408e-12|PASS|

O ensaio tol8 anterior falhou na conservação normalizada SS1: pior
razão M0=4.762852e-6 e pior razão Mprocessada=4.576636e-6. Esse resultado
permanece identificado em `../v3_smokes/`; não foi substituído por tol10.
Aperto de tolerância e igualdade dos fluxos comuns são condições
numéricas; não demonstram, por si, avanço sustentado ou validação física.
As fontes permanecem congeladas e `ACTIVE_SET` é opção experimental.

A auditoria dos contadores após t78.5s encontrou zero passos aceitos
sem nova solução linear nos dois casos. O avanço lento observado não
resultou de aceitar Uold sem atualização pelo teste inicial de resíduo.

D_SS1 não completou seu tempo solicitado de100s.

D_SS3 não completou seu tempo solicitado de110s.
