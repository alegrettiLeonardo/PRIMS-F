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

## Recuperação posterior da evidência histórica

Depois da integração, `prior_evidence.json` foi recuperado do ZIP local fornecido
`PRIMSF_Evidencias_Implementacao.zip`. O SHA-256 coincide exatamente com o manifesto
original: `f88afc885178f1c3456ed0d03c9177566939b7ab90edf77ce77f411ae169af89`.
O arquivo foi restaurado como entrada obrigatória de checksum. Os registros
`NOT_DISTRIBUTED` acima documentam o primeiro checkout e foram preservados;
a lacuna de distribuição está resolvida, sem regeneração da evidência histórica.
O manifesto atualizado contém o hash do ZIP e a proveniência da recuperação.
