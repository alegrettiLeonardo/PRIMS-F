# Envelopes do caso OSC e refinamento temporal

Gerado em 2026-11-18T18:51:59.619966+00:00.

A classificação do sinal muda com o passo de tempo. Os resultados abaixo são observáveis numéricos da pressão na junção; não demonstram autovalores do operador, convergência global, um ciclo limite assintótico ou validação física de severe slugging.

## Método e janelas

Os máximos e mínimos locais são refinados por uma parábola de três amostras, sem extrapolação. Cada pico usa os dois mínimos vizinhos para interpolar a linha de base. A semi-amplitude é metade da diferença entre o pico e essa linha de base. O ajuste principal é a regressão linear de log(amplitude) em **40–400 s**, usando amplitudes ≥ 1 Pa. A separação mínima entre extremos do mesmo tipo é 10 s e somente pares de mínimos separados por 8–40 s são associados a um ciclo. Esses parâmetros medem o sinal observado; não definem um regime físico.

O diagnóstico de amortecimento pode usar sinais pequenos. Para rotular dez ciclos finais como úteis, exige-se excursão pico–mínimo ≥ 100 Pa em todos eles e RMS da pressão em torno da média ≥ 1 Pa. Onze picos consecutivos definem exatamente dez intervalos. O RMS é integrado no tempo sobre a interpolação linear das amostras. Dados ainda em execução são explicitamente identificados.

Os refinados 0,0625 e 0,03125 s também gravam o monitor com intervalo menor. O JSON registra essa cadência e repete o ajuste inicial em uma grade comum de amostragem de 0,125 s, por interpolação linear. Isso permite distinguir o efeito da gravação do efeito do passo do solver; a tabela principal usa as amostras nativas.

## Crescimento inicial observado

| N | Δt máximo (s) | Fim observado (s) | Execução completa | Taxa (s⁻¹) | R² | Picos | Intervalo médio entre picos (s) |
|---:|---:|---:|:---:|---:|---:|---:|---:|
| 51 | 0.5 | 2400.000 | sim | -0.010982774 | 1.000000000 | 15 | 22.999262 |
| 51 | 0.25 | 2400.000 | sim | -0.002002954 | 1.000000000 | 15 | 22.939933 |
| 51 | 0.125 | 2400.000 | sim | +0.002513812 | 0.999999911 | 16 | 22.937111 |
| 51 | 0.0625 | 2400.000 | sim | +0.004769835 | 0.999999408 | 16 | 22.943745 |
| 51 | 0.03125 | 1329.811 | não | +0.005894235 | 0.999998428 | 16 | 22.949669 |
| 101 | 0.5 | 2400.000 | sim | -0.009550168 | 1.000000000 | 15 | 23.039848 |
| 101 | 0.125 | 2400.000 | sim | +0.004150573 | 0.999999776 | 16 | 22.987533 |
| 201 | 0.5 | 2400.000 | sim | -0.008819474 | 1.000000000 | 15 | 23.060673 |
| 201 | 0.125 | 1351.557 | não | +0.004985521 | 0.999999997 | 15 | 23.014067 |

O alto R² descreve a regularidade exponencial desta janela. O erro padrão da regressão não é uma incerteza de modelo: os ciclos vêm de uma única trajetória numérica, e não de réplicas independentes.

## Acomodação nas execuções completas

| N | Δt (s) | RMS 900–1400 / 1400–1900 / 1900–2400 (Pa) | Últimos 10 ciclos úteis | Variação de amplitude | Período médio (s) | CV do período |
|---:|---:|---|:---:|---:|---:|---:|
| 51 | 0.5 | 0.00175339 / 0 / 5.82077e-11 | não | — | — | — |
| 51 | 0.25 | 12.5753 / 4.62835 / 1.68953 | não | — | — | — |
| 51 | 0.125 | 2202.86 / 5057.91 / 6703.66 | sim | +2.3986% | 23.103011 | 8.00491e-05 |
| 51 | 0.0625 | 8765.01 / 9681.93 / 9644.51 | sim | +0.0005% | 23.265630 | 1.28259e-06 |
| 101 | 0.5 | 0.00374788 / 5.82077e-11 / 5.82077e-11 | não | — | — | — |
| 101 | 0.125 | 5819.35 / 8801.26 / 8956.87 | sim | +0.0213% | 23.266396 | 3.40431e-06 |
| 201 | 0.5 | 0.00401923 / 5.82077e-11 / 5.82077e-11 | não | — | — | — |

A estabilidade do período não garante a acomodação da amplitude. O JSON registra separadamente se a amplitude dos últimos dez ciclos varia menos de 1% ou de 5%; essas duas verificações são descrições, não critérios universais de validação ou demonstração de um ciclo limite.

## Por que Euler implícito pode ocultar o crescimento

Para o problema escalar z′ = (σ + iω)z, a amplificação de Euler implícito é

```text
q(h) = 1 / [1 − h(σ + iω)]
g_BE(h) = log|q| / h = −log[(1 − hσ)² + (hω)²] / (2h)
ω_BE(h) = atan2(hω, 1 − hσ) / h
```

Mesmo para σ > 0, há crescimento somente se h < 2σ/(σ² + ω²). A estabilidade absoluta no semiplano esquerdo não garante a reprodução de instabilidades no semiplano direito. Para h pequeno, g_BE ≈ σ + h(σ² − ω²)/2: em um modo oscilatório, o termo temporal pode acrescentar amortecimento suficiente para inverter o sinal observado.

O ajuste descritivo dos envelopes N = 51 com h = 0,5; 0,25; 0,125 s produz σ efetivo = 0.0070326971 s⁻¹ e ω efetivo = 0.26879921 rad/s, com passo de mudança de sinal ≈ 0.194535 s. **Esses parâmetros não são autovalores demonstrados.** O modo LAGGED congela a dissipação de Roe do passo anterior; sua linearização completa é uma separação de operadores, e não necessariamente a fórmula escalar de um único operador fixo.

| Δt (s) | Usado no ajuste | Taxa observada (s⁻¹) | Previsão BE (s⁻¹) | Erro na previsão do período |
|---:|:---:|---:|---:|---:|
| 0.5 | sim | -0.010982774 | -0.010982394 | +1.888% |
| 0.25 | sim | -0.002002954 | -0.002004108 | +1.871% |
| 0.125 | sim | +0.002513812 | +0.002514583 | +1.858% |
| 0.0625 | não | +0.004769835 | +0.004774669 | +1.844% |
| 0.03125 | não | +0.005894235 | +0.005904060 | +1.833% |

Nos passos 0,125; 0,0625; 0,03125 s, a razão entre os incrementos da taxa é 2.006425, equivalente à ordem observada 1.004627 para **esta grandeza e janela**. A extrapolação linear dos dois passos menores dá 0.007018634 s⁻¹. Isso é um indício de erro temporal de primeira ordem no crescimento inicial, sem constituir GCI global ou convergência da amplitude final, do período, do campo transiente ou do regime físico.

## Reprodutibilidade e limites

O JSON contém hashes dos bytes lidos de monitor.csv, run_metadata.txt e run_status.txt, janelas completas/incompletas, R², pisos, valores de cada pico e resultados de todos os ajustes. O CSV osc_envelopes.csv permite reconstruir as curvas de amplitude; osc_envelopes_fits.csv lista as regressões. Os PNG e PDF são figuras independentes para exportação.

As execuções refinadas em curso podem mudar suas estatísticas tardias. Reexecute este script após a conclusão. Nenhum alvo experimental é usado para ajustar as taxas, o modelo descritivo ou os limiares. Captar uma oscilação de pressão não comprova a sequência física de bloqueio, expulsão e retorno de líquido requerida para caracterizar severe slugging.
