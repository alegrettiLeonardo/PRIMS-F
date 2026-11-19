# Chexal: CCFL e propagação de status — 18/11/2026

Correção sobre `b0d33cb25ab2480b3eb0f5d37c18edc340b0f122`.
`summary.json` registra o build release, 3447 verificações/12 programas Fortran,
19 testes Python e os hashes realmente testados. O teste core tem 147 checks.

`probe.f90` consulta a mesma referência ar–água do teste core. Os arquivos
`core_probe_before.csv` e `core_probe_after.csv` mantêm os resultados de CCFL;
células vazias indicam ausência de saída física em status não nulo.
`core_recovery_before.csv` e `core_recovery_after.csv` registram a recuperação
sem dicas exatamente na dobra de `jg=0,1 m/s`, em ambos os ramos. Os hashes
do núcleo antes/depois e do probe estão em `core_probe_source_hashes.json`.

Antes, CCFL com `jg=1e-4 m/s` tinha resíduo relativo 0,225 apesar de status OK;
para `jg<=1e-6 m/s` nas amostras, alpha=0 e resíduo relativo -1. A busca não
resolve esses exemplos e agora retorna NUMERIC. As referências `jg=0,1/0,01`
preservam resultados com resíduo relativo menor que `8e-16`.

Na dobra de `jg=0,1`, HIGH e LOW antes retornavam OK recuperando cocorrente
descendente. Agora propagam SINGULAR, sem trocar de ramo após descartar a
singularidade. As verificações locais de interioridade, resíduo e máximo
constam em `docs/CHEXAL_REMAINING_LIMITS.md`; não são prova global das raízes.

Reprodução do probe atual a partir da raiz do repositório:

```sh
mkdir -p tmp/chexal_ccfl_probe
gfortran -std=f2018 -O0 -ffree-line-length-none -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -Jtmp/chexal_ccfl_probe \
  -Itmp/chexal_ccfl_probe scr/mod_kinds.f90 scr/module_chexal_1992.f90 \
  validation/diagnostics/chexal_ccfl_status_20261118/probe.f90 \
  -o tmp/chexal_ccfl_probe/probe
tmp/chexal_ccfl_probe/probe
```

Para o antes, extraia o núcleo do commit identificado acima para um diretório
temporário e compile o mesmo probe contra essa cópia. Nenhum PDF ou extrato
documental integra esta evidência. Isto é verificação numérica, não validação
experimental.
