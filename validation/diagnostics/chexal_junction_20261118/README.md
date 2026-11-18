# Junção MATCHED e traço gasoso negativo — 18/11/2026

**Resultado:** Chexal não elimina o limite do traço horizontal MATCHED com gás
negativo. Não foi encontrada uma inconsistência algébrica de balanço/pressão nos
estados admissíveis exercitados. Nenhuma projeção, clipping, barreira ou alteração
do operador de junção foi introduzida nesta etapa.

`test_chexal_junction.f90` verifica ID1 no oleoduto e IDs5/6/7 no riser: os dois
fluxos de massa são comuns, pressões reconstruídas separadamente a partir dos
momentos de ambos os traços coincidem, holdups distintos são preservados e o
balanço de duas células com larguras diferentes telescopa. Também verifica
limites de fases ausentes nos contornos de entrada/saída, sem confundir vazão
nula com inventário ausente. Os 108 checks incluem saída vertical contracorrente
nos ramos HIGH/LOW, jl=-0,1 e jg=0,1 m/s, com Pout/Pinterior=0,9999/1/1,0001:
pressão imposta, momento retido, sinais das fases, residual do ramo e ida/volta
na pressão inalterada. A suíte executada passou 10/10 programas.
São verificações de operadores, não demonstração
de condições características bem postas em todo regime.

## Duas evidências distintas

1. O teste unitário reaplica o traço histórico preservado em
   `bg_ss1_stall/summary.json`: jl=0,5296144723724787 m/s e
   jg=-6,800580025319073e-11 m/s. Os três seletores retornam `CHEXAL_DOMAIN=61`
   e zero raízes no traço horizontal, com alfa indefinido.
2. `reproduce_snapshot.py` recompila as fontes e refaz a linearização Newton do
   snapshot completo Bhagwat SS1, mantendo essa lei nas células. Só os traços
   propostos são avaliados adicionalmente em Chexal. Nas cinco propostas com
   fatores 1, 1/2, 1/4, 1/8 e 1/16, todas as 48 células permanecem admissíveis,
   a junção Bhagwat retorna 1012 e os 15 replays Chexal retornam 61. O passo completo
   reexecutado produz jg=-6,7683662081e-11 m/s e vazão gasosa
   -3,7291426317e-13 kg/s. A pequena diferença em relação ao valor histórico é
   preservada: estas são execuções distintas e seus hashes constam dos resumos.

No passo completo, as parcelas central e dissipativa gasosas são respectivamente
1,85455794147e-4 e 1,85455794520e-4 kg/s. O sinal negativo vem da subtração entre
essas parcelas do fluxo numérico. As velocidades gasosas celulares continuam
positivas; não há demonstração de retorno físico do gás. Reproduzir:

```
python3 validation/diagnostics/chexal_junction_20261118/reproduce_snapshot.py --output tmp/replay_new_label
```

## Domínio e limite de modelagem

MATCHED resolve o traço do riser em theta=0 rad desde a horizontal. Chexal1992
p.264 permite apenas cocorrente horizontal; a contracorrente citada na p.262
é vertical, gás subindo e líquido descendo. A escolha HIGH/LOW da p.268 não
transforma o traço horizontal negativo em estado admissível.

O acoplamento atual constrói pressão por uma aproximação central/dissipativa e
em seguida inverte o fechamento. O teste não o transforma numa solução exata
do problema de Riemann da junção. Uma mudança para seguir esse estado exigiria
uma formulação de ondas/condições características ou frente móvel compatível
com as fontes e seus domínios, informação adicional de regime/ramo e validação
própria. Não há justificativa aqui para zerar o fluxo gasoso, trocar espectros,
alargar o domínio Chexal ou alterar automaticamente a viscosidade numérica.
