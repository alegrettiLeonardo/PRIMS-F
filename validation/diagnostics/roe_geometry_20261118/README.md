# Caminho Roe com inclinação variável — 18/11/2026

Verificação isolada sobre `47433f8` (alpha-tail), sem Chexal. O método agora integra
`F_U(U(t), theta(t))`, com interpolação linear tanto de U quanto da inclinação
desde a horizontal. O argumento opcional `theta_right` preserva o caminho antigo
quando omitido. MATCHED já transforma ambos os estados para a lei e inclinação
do oleoduto; esse caminho continua com inclinação constante.

A identidade correta é `F_R-F_L = integral(F_U dt) DeltaU + integral(F_theta dt) DeltaTheta`.
O termo geométrico foi verificado independentemente por diferenças de fluxos a U
fixo e integração Simpson. A matriz 3x3 é uma aproximação dissipativa: não se
afirma uma identidade exata Roe `DeltaF=A DeltaU` com ângulos diferentes.
O fluxo conserva os fluxos físicos das extremidades nas respectivas geometrias
e usa um fluxo comum por face. Nenhum espectro complexo é substituído nem há
clipping de estados.

`test_roe_geometry.f90` verifica a identidade aumentada, detecta a contribuição
geométrica não nula, reversão e reparametrização do caminho, preservação exata
da API de ângulo constante, derivadas completas do fluxo e espectro real dos
estados exercitados. A suíte passou sete programas sobre os fontes identificados
no resumo. Isso não estabelece hiperbolicidade global nem valida slugging.

Regressão congelada SS3 N51+51/dtmax0,25 s: a reconstrução dos endpoints usa os
valores de alfa, pressão e U3 da tabela histórica `roe_path_nodes.csv`, junto
com EOS/escalas do respectivo `run_metadata.txt`. Reproduz a falha14 do caminho
na inclinação média; o caminho aumentado e o fluxo de face retornam sucesso,
com autovalores reais e sem modificar as células ou a geometria do caso.
A correção local não elimina limites celulares ul=0 nem traços gasosos negativos
de MATCHED. Novos transientes e refinamentos são necessários.
