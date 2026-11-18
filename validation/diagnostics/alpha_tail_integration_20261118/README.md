# Integração alpha-tail, 18/11/2026

Base congelada: `05ad799`. O candidato foi revisado e aplicado sem alterações
na física em `scr/module_steady_state.f90`; o hash final é
`76311d3a412b6699ccb25a93c74062b2f993181bbecd9aff0928530a02c01624`.
A busca geométrica só entra após zero raízes da busca antiga, Bhagwat, ângulo
não negativo e fluxos cocorrentes. A margem positiva de C0 justifica unicidade
na cauda; os extremos devem representar um intervalo real. Não há piso de alfa,
clipping nem extensão de contracorrente.

Executado antes e depois da aplicação: O0/O3, 55 guards, 398 testes estacionários
por versão, 101 registros antigos idênticos byte a byte. Revisão O2: 129 contratos
(108 raízes novas, 18 sucessos e três exclusões preservados). Os seis casos antigos
com tolerância frouxa permanecem uma limitação de precisão, não uma aprovação.
`make test FC=gfortran LIBS='-llapack -lblas'` aprovou seis programas; Python CI
aprovou 16 testes. Comandos, hashes, erros e logs estão neste diretório.

Lacuna de empacotamento: `prior_evidence.json` não existe no commit base embora
o manifesto o mencione. `candidate_manifest_historical.json` preserva o manifesto
original; o runner verifica todos os insumos presentes e relata o arquivo
histórico como `NOT_DISTRIBUTED`, sem fabricar seu conteúdo ou validar seu hash.
As notas README/REVIEW do pacote descrevem a revisão histórica e ficam preservadas.
Os resumos desta pasta são novas execuções independentes.

Não foram repetidos transientes neste commit. Os resultados anteriores continuam
associados à baseline sem alpha-tail; os limites de ul=0, MATCHED gasoso negativo,
Roe com inclinação variável e a validação experimental permanecem trabalhos distintos.
