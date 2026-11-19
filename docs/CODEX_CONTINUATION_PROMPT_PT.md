> Nota histórica da base `05ad799`, anterior à integração. Para o estado atual,
> leia [CODEX_HANDOFF_STATE.md](CODEX_HANDOFF_STATE.md) e
> [CHEXAL_1992_IMPLEMENTATION.md](CHEXAL_1992_IMPLEMENTATION.md).

# Prompt de continuação — PRIMS-F / Chexal

Copie o bloco abaixo para o Codex aberto no projeto do WSL. Os PDFs anexados à
conversa precisam estar acessíveis nesse ambiente; o clone não os contém.

```text
Continue a auditoria e implementação Fortran 2018 do PRIMS-F em:
/home/leonardoab/Documentos/PRIMS-F
GitHub: alegrettiLeonardo/PRIMS-F; PR #1 já incorporado ao main.
Base de continuação: origin/audit/chexal-continuation-20261118.

1. Leia AGENTS.md e git status. Preserve alterações locais. Confira/fetch a branch
   acima; se houver divergência, use worktree separado, sem reset ou checkout
   destrutivo. Leia primeiro docs/CODEX_HANDOFF_STATE.md; consulte outros arquivos
   com rg conforme a tarefa. Não recomece toda a auditoria nem despeje PDFs/logs
   completos no contexto. Registre decisões e resultados em arquivos concisos.

2. Antes de mudar física, revise e integre em commit separado o candidato
   validation/diagnostics/alpha_tail_candidate/ (ainda NÃO aplicado): corrige
   NO_ROOT para alfa positivo minúsculo, sem estender o domínio de Bhagwat.
   Reproduza seus testes e mantenha a evidência anterior identificada por versão.

3. Objetivo: implementar CHEXAL_1992 como fechamento selecionável, preservando
   Bhagwat para regressão. Use chexal1992.pdf, a dissertação, sistema_oleoduto_riser
   e Baliño2010; localize os PDFs, cite equações/páginas e não publique seus
   conteúdos no GitHub. Se faltar fonte necessária, registre a lacuna e avance
   nas tarefas independentes, sem inventar coeficientes.

4. Respeite o domínio da fonte: Chexal1992 p.262 admite contracorrente VERTICAL,
   jg>0 e jl<0; p.264 restringe horizontal a cocorrente; p.268 admite duas raízes
   e exige informação adicional para seleção. Não misture versões 1992/1997 nem
   habilite contracorrente inclinada sem referência primária específica. A
   catenária atual não é vertical. Use parâmetros ar–água, Reynolds assinados e
   conversão explícita do ângulo (Chexal: graus desde a vertical).

5. Formule R=alpha*[C0*(jl+jg)+Vgj]-jg, sem divisão por j quando j=0. Implemente
   seleção/continuação explícita de ramo, CCFL, limites monofásicos e status de
   domínio/singularidade. Corrija conjuntamente recuperação de primitivas, fluxo,
   derivadas/Jacobianos, estacionário, contornos e junção. Não use clipping,
   valores absolutos indiscriminados ou raízes/espectros artificiais. Examine
   separadamente o traço gasoso negativo MATCHED/Roe: Chexal não o resolve
   automaticamente. Revise também o caminho de Roe com inclinação variável.

6. Execute make test FC=gfortran LIBS='-llapack -lblas' e os testes Python da CI.
   Acrescente referências independentes para Chexal: cocorrente/contracorrente
   vertical, duas raízes/CCFL, j=0, fases ausentes, sinais, unidades, invariância
   de escala e derivadas no mesmo ramo. Teste ida/volta U↔primitivas, conservação
   e hiperbolicidade; não trate passagem unitária como validação experimental.

7. Repita os dois refinamentos interrompidos em worktree separado, com a
   baseline congelada, sem alpha-tail/Chexal, para fechar a campanha anterior. Após os testes locais, execute controle estável e casos instáveis, inicialmente
   2400 s e depois até cerca de dez ciclos acomodados. Refine malha e passo
   separadamente; preserve EOS, vazões e geometria na comparação de fechamentos.
   Registre pressão, fração líquida, fluxos por fase, inventários e balanço integral.
   Analise acumulação/compressão, produção, penetração e blowdown conforme o regime;
   bloqueio completo de gás não é requisito universal de OSC/SS3. Documente a
   divergência SS3: dissertação 26 s; Baliño2010 tabela3/p.653, caso11, 76 s.
Entregue commits pequenos, testes realmente executados, métricas, fontes e limites
restantes. Mantenha o PR draft até a validação; não ajuste parâmetros para forçar
slugging nem declare sucesso por oscilações de pressão isoladas. Prossiga com
autonomia nas correções locais autorizadas e faça perguntas apenas se indispensável.
```
