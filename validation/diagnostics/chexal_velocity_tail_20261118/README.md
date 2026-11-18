# Correção da derivada com velocidade positiva pequena

A versão congelada 7321e67 limitava a redução do stencil centrado a16 passos.
Quando ul se aproximava de2,2586e-8 m/s, as perturbações cruzavam ul=0 e podiam
recusar estados admissíveis. A nova avaliação usa stencil unilateral de ordem4,
assinaturas de peça e convergência, sem modificar correlação ou domínio físico.
A revisão também corrigiu uma recusa espúria da guarda bilateral fixa junto a
Rel=Reg: seu raio passa a diminuir antes de rejeitar o estado.

Executados: build release; make test FC=gfortran LIBS='-llapack -lblas'; e os
19 testes Python da CI. Resultado: **3439 verificações em12 programas Fortran**,
incluindo87 de cauda/guarda e21 contra referências independentes mpmath60.
O resumo identifica hashes finais, compilador, flags e commit pai. Não é
validação experimental. A evidência anterior3331/10 permanece preservada em
chexal_integration_20261118.

Os replays em chexal_velocity_tail_replay_20261118 e CSVs independentes em
chexal_velocity_reference_20261118 identificam por hash a primeira versão
intermediária, antes do ajuste da guarda. Não são atribuídos retroativamente
ao commit final; novas campanhas terão outra versão congelada.
