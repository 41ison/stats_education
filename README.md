# Sinal ou ruído? Significância estatística em múltiplos testes

Dashboard educativo em R/Shiny que mostra, de forma interativa, por que o
p-valor ajustado deve ser usado no lugar do p-valor nominal quando milhares
de testes são feitos ao mesmo tempo, como em proteômica e transcriptômica.

A ideia é uma experiência de "revelação": o volcano plot começa cheio de
candidatos em vermelho e, ao ativar o ajuste FDR, os pontos que não resistem
à correção caem para baixo da linha de corte e ficam cinza.

## O que o app mostra

- **Volcano plot animado:** `log2(Fold Change)` no eixo X e `-log10(p-valor)`
  no eixo Y, com linha de corte em p = 0,05.
- **Interruptor "Aplicar Filtro de Acasos (Ajuste FDR)":** alterna o eixo Y
  entre o p-valor nominal e o p-valor ajustado por Benjamini-Hochberg.
- **Tooltips educativos:** ao passar o mouse sobre um ponto, o app mostra o
  p-valor nominal e o ajustado daquela proteína e o que aconteceu com ela.
- **Caixas de resumo:** candidatos iniciais, sobreviventes ao ajuste e
  quantos dos pontos em destaque não têm efeito real.
- **Box de justificativa:** texto explicativo com os números da simulação
  atual e uma tabela comparando os critérios nominal, FDR e Bonferroni.

## Sobre os dados

Os dados são **simulados** a cada execução. O app gera um experimento com
dois grupos (controle e tratado), no qual uma fração das proteínas tem efeito
real e o restante é apenas ruído, e aplica um teste t de Welch por proteína.
Como a verdade é conhecida, é possível contar falsos positivos e efeitos
reais encontrados por cada critério.

Na barra lateral é possível ajustar:

| Controle | Intervalo | Padrão |
|---|---|---|
| Proteínas testadas | 500 a 5000 | 3000 |
| % que realmente mudam | 0 a 20 | 5 |
| Réplicas por grupo | 3 a 10 | 6 |

O botão "Nova simulação" sorteia um novo experimento com os mesmos
parâmetros.

## Como executar

Requer R (4.1 ou superior) e os pacotes `shiny`, `bslib` e `plotly`:

```r
install.packages(c("shiny", "bslib", "plotly"))
```

Clone o repositório e rode o app a partir da pasta do projeto:

```r
shiny::runApp("app.R")
```

No Positron ou no RStudio, também é possível abrir `app.R` e usar o botão de
execução do app.

## Sugestões de uso em aula

1. Comece com o interruptor desligado e pergunte quantas proteínas "mudaram".
2. Ative o ajuste FDR e observe quantos candidatos caem.
3. Coloque "% que realmente mudam" em 0 e mostre que os candidatos nominais
   continuam aparecendo, mesmo sem nenhum efeito real.
4. Aumente o número de réplicas e veja mais efeitos reais resistirem ao
   ajuste.
5. Compare os três critérios na tabela do box de justificativa.

## Licença

Distribuído sob a licença MIT. Veja o arquivo [LICENSE](LICENSE).
