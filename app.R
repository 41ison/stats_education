library(shiny)
library(bslib)
library(plotly)

alpha <- 0.05

col_red <- "rgba(228, 26, 28, 0.85)"
col_blue <- "rgba(55, 126, 184, 0.85)"
col_small <- "rgba(230, 159, 0, 0.75)"
# Conventional effect-size cutoff: at least a doubling or a halving
fc_cut <- 1
col_fallen <- "rgba(110, 110, 110, 0.30)"
col_never <- "rgba(170, 180, 195, 0.35)"
# Faint red outline marks greys that used to be red
col_fallen_line <- "rgba(228, 26, 28, 0.45)"
col_no_line <- "rgba(255, 255, 255, 0)"

# Simulated two-group experiment: a few proteins truly change, the rest are noise
simulate_experiment <- function(n, prop_true, reps) {
  is_true <- seq_len(n) %in% sample.int(n, round(n * prop_true))
  effect <- ifelse(
    is_true,
    sample(c(-1, 1), n, replace = TRUE) * runif(n, 1, 2.5),
    0
  )
  sd <- runif(n, 0.3, 0.8)

  control <- matrix(rnorm(n * reps, 0, sd), nrow = n)
  treated <- matrix(rnorm(n * reps, effect, sd), nrow = n)

  # Vectorised Welch t-test, one per protein
  m1 <- rowMeans(control)
  m2 <- rowMeans(treated)
  se1 <- rowSums((control - m1)^2) / (reps - 1) / reps
  se2 <- rowSums((treated - m2)^2) / (reps - 1) / reps
  t_stat <- (m2 - m1) / sqrt(se1 + se2)
  df <- (se1 + se2)^2 / (se1^2 / (reps - 1) + se2^2 / (reps - 1))
  p <- 2 * pt(-abs(t_stat), df)

  data.frame(
    protein = sprintf("P%04d", seq_len(n)),
    is_true = is_true,
    log2fc = m2 - m1,
    p_nominal = p,
    p_adjusted = p.adjust(p, method = "BH")
  )
}

fmt_p <- function(p) {
  ifelse(
    p < 0.001,
    formatC(p, format = "e", digits = 1),
    formatC(p, format = "f", digits = 3)
  )
}

# Everything about the main trace that differs between the two views
view_state <- function(d, adjusted) {
  candidate <- d$p_nominal < alpha
  survivor <- d$p_adjusted < alpha
  head <- sprintf("<b>Proteína %s</b><br>", d$protein)

  if (!adjusted) {
    return(list(
      y = -log10(d$p_nominal),
      color = ifelse(candidate, col_red, col_never),
      line = rep(col_no_line, nrow(d)),
      text = paste0(
        head,
        ifelse(
          candidate,
          sprintf(
            "Candidata inicial (p-nominal = %s).<br>Ainda sem correção para múltiplos testes.",
            fmt_p(d$p_nominal)
          ),
          sprintf(
            "Não estatisticamente significativa (p-nominal = %s).",
            fmt_p(d$p_nominal)
          )
        )
      ),
      y_title = "-log10(p-valor nominal)",
      legend_red = "Candidatos Iniciais (p-nominal < 0,05)"
    ))
  }

  fallen <- candidate & !survivor
  list(
    y = -log10(d$p_adjusted),
    color = ifelse(survivor, col_red, ifelse(fallen, col_fallen, col_never)),
    line = ifelse(fallen, col_fallen_line, col_no_line),
    text = paste0(
      head,
      ifelse(
        survivor,
        sprintf(
          "Sobreviveu ao ajuste (p-nominal = %s, p-ajustado = %s).",
          fmt_p(d$p_nominal),
          fmt_p(d$p_adjusted)
        ),
        ifelse(
          fallen,
          sprintf(
            "Parecia promissora (p-nominal = %s),<br>mas não se distingue do ruído estatístico (p-ajustado = %s).",
            fmt_p(d$p_nominal),
            fmt_p(d$p_adjusted)
          ),
          sprintf(
            "Nunca foi candidata (p-nominal = %s, p-ajustado = %s).",
            fmt_p(d$p_nominal),
            fmt_p(d$p_adjusted)
          )
        )
      )
    ),
    y_title = "-log10(p-valor ajustado, FDR)",
    legend_red = "Sobreviveram ao ajuste (p ajustado < 0,05)"
  )
}

ref_table <- function(header, rows) {
  tags$table(
    class = "table table-sm align-middle w-auto",
    tags$thead(tags$tr(lapply(header, tags$th))),
    tags$tbody(lapply(rows, \(r) tags$tr(lapply(r, tags$td))))
  )
}

# Static text of the "Anatomia do volcano plot" tab
anatomy_md <- list(
  reading = "
Cada ponto é uma proteína, e o gráfico responde a duas perguntas ao mesmo tempo:

- **Eixo horizontal: quanto mudou?** É o tamanho do efeito, em log2(fold change).
  Pontos à direita são mais abundantes no grupo tratado; à esquerda, menos
  abundantes; no centro, praticamente iguais nos dois grupos.
- **Eixo vertical: quão difícil é atribuir essa mudança ao acaso?** É a evidência
  estatística, em −log10(p-valor). Quanto mais alto o ponto, menor o p-valor.

O nome vem da forma. A maioria das proteínas muda pouco e tem p-valor alto,
formando a base larga do vulcão. Efeitos grandes tendem a vir com p-valores
pequenos, e por isso os pontos \"espirram\" para cima e para os lados.

- **Cantos superiores:** mudança grande *e* estatisticamente significativa. São os
  candidatos mais interessantes.
- **Topo, perto do centro:** estatisticamente significativa, mas pequena. Pode ser
  real e ainda assim irrelevante do ponto de vista biológico.
- **Laterais, embaixo:** mudança aparentemente grande, mas com tanta variação entre
  as réplicas que não dá para separá-la do acaso.
- **Base central:** nem mudança, nem evidência.

Os dois eixos não são independentes: o p-valor é calculado a partir da própria
diferença entre os grupos, dividida pela sua incerteza. É por isso que o centro
do gráfico fica vazio no alto.
",
  double_cut = "
Na prática, uma proteína só costuma entrar na lista de resultados quando passa
por **dois cortes ao mesmo tempo**:

- **Linha horizontal, evidência:** p-valor (de preferência o ajustado) abaixo de
  0,05. A mudança é estatisticamente significativa.
- **Linhas verticais, tamanho do efeito:** |log2FC| ≥ 1. A abundância pelo
  menos dobrou ou caiu pelo menos pela metade.

Cada corte responde a uma pergunta diferente, e nenhum substitui o outro. Com
muitas réplicas ou pouca variação, diferenças minúsculas chegam a p-valores
pequenos: são estatisticamente significativas, mas talvez não importem. Com
poucas réplicas, fold changes grandes aparecem por puro acaso: impressionam,
mas não se sustentam.

No gráfico acima, os cortes dividem o plano em regiões:

- **Vermelho e azul:** passam nos dois cortes.
- **Laranja:** estatisticamente significativos, mas com efeito menor que o
  corte.
- **Cinza:** não passam no corte de p-valor, qualquer que seja o fold change.

**Cuidados.** O ±1 é uma convenção, não uma lei: para algumas proteínas, uma
mudança de 30% pode ter consequência biológica, e o corte adequado depende do
sistema e da técnica de medida. Os cortes devem ser definidos antes de olhar os
dados. E a garantia do FDR vale para o conjunto selecionado pelo p-valor
ajustado; ao filtrar depois pelo fold change, o conjunto muda e a garantia não
se transfere automaticamente.
",
  p_value = "
**O que é.** O p-valor é a probabilidade de observar uma diferença tão grande
quanto a observada, ou maior, *supondo que não exista diferença real* entre os
grupos (a hipótese nula). Um p-valor pequeno diz que os dados seriam
surpreendentes se tudo fosse acaso. Ele **não** é a probabilidade de a hipótese
nula ser verdadeira, nem mede o tamanho ou a importância do efeito.

**De onde vem a palavra.** Em 1885, **Francis Edgeworth** usou \"significant\"
para dizer que uma diferença *significava* (indicava) algo além de uma
flutuação acidental. Era \"significar\" no sentido de sinalizar, não de ser
importante.

**De onde vem o P.** Em 1900, **Karl Pearson** publicou o teste qui-quadrado e,
com ele, o cálculo de *P*: a probabilidade de desvios tão grandes ou maiores
surgirem por amostragem aleatória.

**De onde vem o 0,05.** Em 1925, em *Statistical Methods for Research Workers*,
**Ronald Fisher** popularizou os testes de significância e sugeriu 1 em 20
(0,05) como um limite conveniente. No seu famoso experimento da senhora que
prova chá, o resultado também *significava* algo no sentido de ajudar a
interpretar o que se observou. Para Fisher, o p-valor era uma medida
contínua de evidência contra a hipótese nula, a ser pesada pelo pesquisador, e
não um veredito automático.

**De onde vem a decisão \"significativo / não significativo\".** Entre 1928 e
1933, **Jerzy Neyman** e **Egon Pearson** (filho de Karl) propuseram outra
lógica: fixar *antes* do experimento a taxa de erro aceitável (α), considerar
uma hipótese alternativa e distinguir dois erros, o falso positivo (tipo I) e o
falso negativo (tipo II). Daí vem o conceito de poder do teste. Fisher discordou
publicamente dessa abordagem até o fim da vida.

**O que usamos hoje** é um híbrido das duas escolas: calcula-se o p-valor de
Fisher e compara-se com o α de Neyman-Pearson. A linha tracejada do volcano plot
é esse α; a altura de cada ponto é o p-valor.

**Um problema de linguagem.** No sentido de Edgeworth, \"significativo\" quer
dizer apenas \"sinaliza algo que provavelmente não é acaso\". Em *\"Why say
'statistically significant' rather than just 'significant'? A plea to rid the
medical literature of linguistic ambiguity\"*, Fingerhut e colaboradores (2022)
contam que, com o uso disseminado por Karl Pearson e seus discípulos, esse
sentido foi se deslocando aos poucos, até a palavra e o próprio p-valor
passarem, incorretamente, por sinônimos de \"importante\".

Os autores dão duas razões para nunca escrever \"significativo\" sozinho em
textos científicos. Primeira: o p-valor é o resultado de um teste *estatístico*,
e não uma propriedade do efeito ou da população estudada. Segunda: sem o
modificador, quem lê guarda, sem perceber, o sentido de \"importante\" ou
\"comprovado\". Daí a recomendação: diga-se sempre **estatisticamente
significativo**. No volcano plot essa distinção é visível: a significância
estatística está no eixo vertical, e a relevância do efeito, no horizontal.
",
  p_adjusted = "
O raciocínio de Fisher e de Neyman-Pearson foi pensado para **um** teste. Um
volcano plot resume milhares. Com α = 0,05, 5% das proteínas sem efeito algum
ficam acima da linha de corte só por acaso.

O **p-valor ajustado** corrige cada p-valor pelo número de testes realizados.
Com o método de Benjamini-Hochberg (1995), o corte de 0,05 passa a controlar a
*taxa de falsas descobertas* (FDR): entre os pontos acima da linha, espera-se
que no máximo 5% sejam falsos positivos.

No gráfico, o ajuste nunca diminui um p-valor, então os pontos descem ou ficam
onde estão. A ordem entre eles se mantém, e a posição horizontal não muda. Use o
interruptor na barra lateral para ver isso acontecer; a aba Simulação detalha os
números.
",
  fold_change = "
O **fold change (FC)** é a razão entre as médias dos dois grupos: tratado
dividido por controle. FC = 2 significa o dobro; FC = 0,5, a metade.

Na escala original, essa razão é assimétrica. Todas as reduções ficam espremidas
entre 0 e 1, enquanto os aumentos vão de 1 ao infinito. Dobrar e cair pela
metade são mudanças de mesma magnitude, mas ficariam a distâncias diferentes do
centro.

O logaritmo resolve isso:

- **Simetria.** log2(2) = +1 e log2(0,5) = −1. Aumentos e reduções equivalentes
  ficam à mesma distância do zero, e \"sem mudança\" fica exatamente no centro.
- **Leitura direta.** Com base 2, cada unidade é uma duplicação: +1 é o dobro,
  +2 é o quádruplo, +3 é oito vezes.
- **Razão vira diferença.** log2(tratado / controle) = log2(tratado) −
  log2(controle). Como abundâncias costumam ser analisadas em escala log, onde a
  variação é mais estável, o log2FC é simplesmente a diferença entre as médias.
",
  neg_log10 = "
Os p-valores interessantes são os pequenos, e na escala original todos eles se
amontoam junto ao zero: 0,01 e 0,000001 seriam indistinguíveis no gráfico,
embora representem evidências muito diferentes.

- **O log10 separa as ordens de grandeza.** Cada unidade no eixo corresponde a
  um p-valor 10 vezes menor.
- **O sinal de menos inverte a escala.** O logaritmo de um número entre 0 e 1 é
  negativo, então sem a inversão os resultados mais fortes ficariam embaixo. Com
  ela, mais evidência fica mais alto, que é como lemos gráficos naturalmente.
",
  references = "
- Fingerhut A, Wexner S, Behrns K, et al. Why say \"statistically significant\"
  rather than just \"significant\"? A plea to rid the medical literature of
  linguistic ambiguity. *Surgery*. 2022;172:1039-1040.
  doi:10.1016/j.surg.2022.08.019. Publicação secundária, em acesso aberto:
  *Annals of Laparoscopic and Endoscopic Surgery*. 2023;8:8.
  doi:10.21037/ales-22-78.
- Edgeworth FY. Methods of statistics. *Journal of the Statistical Society of
  London*, Jubilee Volume. 1885:181-217.
- Pearson K. On the criterion that a given system of deviations from the
  probable in the case of a correlated system of variables is such that it can
  be reasonably supposed to have arisen from random sampling. *Philosophical
  Magazine*. 1900;50:157-175.
- Fisher RA. *Statistical Methods for Research Workers*. Edinburgh: Oliver and
  Boyd; 1925.
- Neyman J, Pearson ES. On the problem of the most efficient tests of
  statistical hypotheses. *Philosophical Transactions of the Royal Society A*.
  1933;231:289-337.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and
  powerful approach to multiple testing. *Journal of the Royal Statistical
  Society B*. 1995;57:289-300.
"
)

ui <- page_sidebar(
  title = "Sinal ou ruído? Significância estatística em múltiplos testes",
  theme = bs_theme(version = 5),
  # Scrolling page, so the plot card keeps its fixed height above the text card
  fillable = FALSE,
  sidebar = sidebar(
    width = 330,
    div(
      class = "p-3 rounded border border-danger-subtle bg-danger-subtle fw-semibold",
      input_switch(
        "fdr",
        "Aplicar Filtro de Acasos (Ajuste FDR)",
        value = FALSE
      )
    ),
    p(
      class = "small text-muted",
      "Desligado: cada proteína é julgada pelo seu p-valor nominal, como se fosse o único teste do experimento.",
      "Ligado: os p-valores são ajustados (Benjamini-Hochberg) para o fato de que milhares de testes foram feitos ao mesmo tempo."
    ),
    accordion(
      open = FALSE,
      accordion_panel(
        "Configurar a simulação",
        sliderInput(
          "n",
          "Proteínas testadas",
          min = 500,
          max = 5000,
          value = 3000,
          step = 500
        ),
        sliderInput(
          "prop_true",
          "% que realmente mudam",
          min = 0,
          max = 20,
          value = 5,
          step = 1
        ),
        sliderInput(
          "reps",
          "Réplicas por grupo",
          min = 3,
          max = 10,
          value = 6,
          step = 1
        )
      )
    ),
    actionButton("resample", "Nova simulação", class = "btn-outline-secondary"),
    p(
      class = "small text-muted",
      "Os dados são simulados, então sabemos quais proteínas mudaram de verdade.",
      "Experimente 0% de mudanças reais: os candidatos iniciais continuam aparecendo."
    )
  ),
  navset_underline(
    id = "tab",
    nav_panel(
      "Simulação",
      class = "pt-3",
      layout_column_wrap(
        width = 1 / 3,
        fill = FALSE,
        value_box(
          title = "Candidatos iniciais (p nominal < 0,05)",
          value = textOutput("n_candidates"),
          p(textOutput("expected_chance", inline = TRUE)),
          theme = "danger"
        ),
        value_box(
          title = "Sobrevivem ao ajuste FDR",
          value = textOutput("n_survivors"),
          p(textOutput("n_fallen", inline = TRUE)),
          theme = "success"
        ),
        value_box(
          title = "Em destaque agora, mas sem efeito real",
          value = textOutput("n_false"),
          p(textOutput("false_share", inline = TRUE)),
          theme = "secondary"
        )
      ),
      card(
        full_screen = TRUE,
        height = "750px",
        card_header(textOutput("plot_title")),
        plotlyOutput("volcano"),
        card_footer(class = "small", textOutput("explanation"))
      ),
      card(
        card_header("Por que usar o p-valor ajustado em vez do nominal?"),
        uiOutput("justification")
      )
    ),
    nav_panel(
      "Anatomia do volcano plot",
      class = "pt-3",
      card(
        full_screen = TRUE,
        height = "600px",
        card_header("Como ler o gráfico"),
        plotlyOutput("anatomy"),
        card_footer(
          class = "small",
          "Mesmos dados da aba Simulação. O eixo vertical acompanha o interruptor da barra lateral (p-valor nominal ou ajustado)."
        )
      ),
      card(
        card_header("Interpretação gráfica"),
        markdown(anatomy_md$reading)
      ),
      card(
        card_header("O critério duplo: p-valor e fold change"),
        markdown(anatomy_md$double_cut)
      ),
      card(
        card_header("p-valor: uma breve história da palavra \"significativo\""),
        markdown(anatomy_md$p_value)
      ),
      card(
        card_header("p-valor ajustado"),
        markdown(anatomy_md$p_adjusted)
      ),
      card(
        card_header("Fold change e por que usar log2"),
        markdown(anatomy_md$fold_change),
        ref_table(
          c("Mudança no tratado", "Fold change", "log2FC"),
          list(
            c("8 vezes maior", "8", "+3"),
            c("2 vezes maior", "2", "+1"),
            c("sem mudança", "1", "0"),
            c("metade", "0,5", "−1"),
            c("um oitavo", "0,125", "−3")
          )
        )
      ),
      card(
        card_header("Por que −log10 no eixo do p-valor"),
        markdown(anatomy_md$neg_log10),
        ref_table(
          c("p-valor", "−log10(p)", "Leitura"),
          list(
            c("1", "0", "nenhuma evidência contra o acaso"),
            c("0,05", "1,30", "linha de corte convencional"),
            c("0,01", "2", ""),
            c("0,001", "3", ""),
            c("0,000001", "6", "evidência muito forte")
          )
        )
      ),
      card(
        card_header("Referências"),
        markdown(anatomy_md$references)
      )
    )
  )
)

server <- function(input, output, session) {
  experiment <- reactive({
    input$resample
    simulate_experiment(input$n, input$prop_true / 100, input$reps)
  })

  counts <- reactive({
    d <- experiment()
    candidate <- d$p_nominal < alpha
    survivor <- d$p_adjusted < alpha
    highlighted <- if (input$fdr) survivor else candidate
    list(
      n = nrow(d),
      candidates = sum(candidate),
      survivors = sum(survivor),
      highlighted = sum(highlighted),
      false_highlighted = sum(highlighted & !d$is_true),
      true_total = sum(d$is_true),
      true_found = sum(highlighted & d$is_true)
    )
  })

  # Independent of the switch: compares all three decision rules on the same data
  methods <- reactive({
    d <- experiment()
    rules <- list(
      "Nominal (sem ajuste)" = d$p_nominal < alpha,
      "FDR (Benjamini-Hochberg)" = d$p_adjusted < alpha,
      "Bonferroni" = p.adjust(d$p_nominal, method = "bonferroni") < alpha
    )
    data.frame(
      method = names(rules),
      highlighted = vapply(rules, sum, integer(1)),
      true_found = vapply(rules, \(r) sum(r & d$is_true), integer(1)),
      false_pos = vapply(rules, \(r) sum(r & !d$is_true), integer(1)),
      row.names = NULL
    )
  })

  output$justification <- renderUI({
    m <- nrow(experiment())
    true_total <- sum(experiment()$is_true)
    tab <- methods()
    nominal <- tab[1, ]
    bh <- tab[2, ]
    bonf <- tab[3, ]

    p_any <- 1 - (1 - alpha)^m
    p_any_txt <- if (p_any > 0.999) {
      "mais de 99,9%"
    } else {
      sprintf("%s%%", format(round(100 * p_any, 1), decimal.mark = ","))
    }
    share_false <- function(row) {
      if (row$highlighted == 0) {
        "–"
      } else {
        sprintf("%.0f%%", 100 * row$false_pos / row$highlighted)
      }
    }

    intro <- sprintf(
      "
**O p-valor nominal responde a uma pergunta que não fizemos.** Ele mede a
chance de ver um resultado tão extremo *se esta fosse a única proteína
testada* e nada tivesse mudado. Com o corte em 0,05, aceitamos 5%% de chance
de alarme falso em **um** teste.

**Em muitos testes, o alarme falso deixa de ser exceção e vira rotina.**
Se nenhuma proteína mudasse, 5%% delas ainda teriam p < 0,05: nos **%d**
testes desta simulação, cerca de **%.0f** \"descobertas\" surgiriam por puro
acaso. A chance de haver ao menos um falso positivo, 1 − 0,95<sup>m</sup>
para testes independentes, seria de **%s**. Aqui, **%d** proteínas têm
p nominal < 0,05, mas **%d** delas não têm efeito real algum.

**O p-valor ajustado leva em conta quantos testes foram feitos.** O método de
Benjamini-Hochberg ordena os p-valores do menor para o maior e multiplica
cada um por m / posição (número de testes dividido pela posição na lista).
O corte de 0,05 passa a controlar a *taxa de falsas descobertas* (FDR):
entre as proteínas declaradas estatisticamente significativas, espera-se que no máximo 5%%
sejam falsos positivos.

**Por que FDR e não uma correção mais dura?** A correção de Bonferroni
multiplica todo p-valor por m (aqui, por %d) e busca evitar *qualquer* falso
positivo. O preço é perder efeitos reais. O FDR aceita uma pequena fração de
erros em troca de mais poder, o que combina com estudos exploratórios cujos
resultados serão validados depois. Nesta simulação, com **%d** proteínas que
realmente mudaram:
",
      m,
      m * alpha,
      p_any_txt,
      nominal$highlighted,
      nominal$false_pos,
      m,
      true_total
    )

    comparison <- tags$table(
      class = "table table-sm align-middle w-auto",
      tags$thead(tags$tr(
        tags$th("Critério (corte em 0,05)"),
        tags$th(class = "text-end", "Pontos destacados"),
        tags$th(
          class = "text-end",
          sprintf("Efeitos reais encontrados (de %d)", true_total)
        ),
        tags$th(class = "text-end", "Falsos positivos"),
        tags$th(class = "text-end", "% de falsos entre os destacados")
      )),
      tags$tbody(lapply(seq_len(nrow(tab)), function(i) {
        tags$tr(
          tags$td(tab$method[i]),
          tags$td(class = "text-end", tab$highlighted[i]),
          tags$td(class = "text-end", tab$true_found[i]),
          tags$td(class = "text-end", tab$false_pos[i]),
          tags$td(class = "text-end", share_false(tab[i, ]))
        )
      }))
    )

    extra_true <- bh$true_found - bonf$true_found
    verdict <- if (extra_true > 0) {
      sprintf(
        "Os dois removem quase todo o ruído do critério nominal, mas o FDR preserva **%d** a mais dos efeitos reais.",
        extra_true
      )
    } else {
      "Nesta simulação, os dois métodos encontram o mesmo número de efeitos reais."
    }

    closing <- sprintf(
      "
Pontos destacados: **%d** com Bonferroni, **%d** com FDR e **%d** sem ajuste.
%s

**O que o ajuste não faz.** Ele não garante que cada ponto destacado seja
real (5%% de falsas descobertas esperadas não é zero), e efeitos verdadeiros
porém fracos podem não sobreviver. Ele também não substitui réplicas:
aumente o número de réplicas na simulação e veja mais efeitos reais
resistirem ao ajuste.
",
      bonf$highlighted,
      bh$highlighted,
      nominal$highlighted,
      verdict
    )

    tagList(markdown(intro), comparison, markdown(closing))
  })

  output$n_candidates <- renderText(counts()$candidates)
  output$expected_chance <- renderText({
    sprintf(
      "Só por acaso, esperaríamos cerca de %.0f em %d testes mesmo que nada tivesse mudado.",
      counts()$n * alpha,
      counts()$n
    )
  })
  output$n_survivors <- renderText(counts()$survivors)
  output$n_fallen <- renderText({
    sprintf(
      "%d candidatos não resistem ao ajuste.",
      counts()$candidates - counts()$survivors
    )
  })
  output$n_false <- renderText(counts()$false_highlighted)
  output$false_share <- renderText({
    k <- counts()
    if (k$highlighted == 0) {
      return("Nenhum ponto em destaque.")
    }
    sprintf(
      "%.0f%% dos %d pontos vermelhos são falsos positivos.",
      100 * k$false_highlighted / k$highlighted,
      k$highlighted
    )
  })

  output$plot_title <- renderText({
    if (input$fdr) {
      "Visão ajustada: o que resiste ao filtro de acasos"
    } else {
      "Visão nominal: ilusão de sucesso"
    }
  })

  output$explanation <- renderText({
    k <- counts()
    if (input$fdr) {
      sprintf(
        "Após o ajuste, restam %d pontos vermelhos; %d deles são efeitos reais (de %d existentes na simulação). O ajuste controla a proporção esperada de falsas descobertas em 5%%: não a zera, e efeitos reais fracos também podem cair. Passe o mouse sobre os pontos cinza com contorno vermelho para ver quem caiu.",
        k$highlighted,
        k$true_found,
        k$true_total
      )
    } else {
      sprintf(
        "%d proteínas parecem ter mudado, mas apenas %d mudaram de verdade nesta simulação. Com %d testes a 5%%, dezenas de p-valores pequenos surgem por puro acaso. Ative o filtro na barra lateral.",
        k$candidates,
        k$true_total,
        k$n
      )
    }
  })

  # Rendered once per simulation; the switch only animates the existing plot
  output$volcano <- renderPlotly({
    d <- experiment()
    s <- view_state(d, isolate(input$fdr))
    x_max <- max(abs(d$log2fc)) * 1.05
    y_max <- max(-log10(d$p_nominal)) * 1.05

    # Legend entries are dummy points parked outside the fixed axis range,
    # so the real points can stay in one trace and animate with stable identity
    legend_trace <- function(p, name, color, line = col_no_line) {
      add_trace(
        p,
        x = 0,
        y = -10,
        type = "scatter",
        mode = "markers",
        name = name,
        hoverinfo = "skip",
        marker = list(
          size = 9,
          color = color,
          line = list(color = line, width = 1)
        )
      )
    }

    plot_ly() |>
      add_trace(
        x = d$log2fc,
        y = s$y,
        type = "scatter",
        mode = "markers",
        text = s$text,
        hoverinfo = "text",
        showlegend = FALSE,
        marker = list(
          size = 6,
          color = s$color,
          line = list(color = s$line, width = 1)
        )
      ) |>
      legend_trace(s$legend_red, col_red) |>
      legend_trace(
        "Caíram com o ajuste (antes vermelhos)",
        col_fallen,
        col_fallen_line
      ) |>
      legend_trace("Não estatisticamente significativos", col_never) |>
      layout(
        xaxis = list(title = "log2(Fold Change)", range = c(-x_max, x_max)),
        yaxis = list(
          title = list(text = s$y_title),
          range = c(-0.05, y_max),
          zeroline = FALSE
        ),
        legend = list(orientation = "h", y = 1.08),
        shapes = list(list(
          type = "line",
          xref = "paper",
          x0 = 0,
          x1 = 1,
          y0 = -log10(alpha),
          y1 = -log10(alpha),
          line = list(dash = "dash", width = 1.5, color = "#444444")
        )),
        annotations = list(list(
          xref = "paper",
          x = 1,
          xanchor = "right",
          y = -log10(alpha),
          yanchor = "bottom",
          text = "corte: p = 0,05",
          showarrow = FALSE,
          font = list(size = 12, color = "#444444")
        ))
      ) |>
      config(displayModeBar = FALSE)
  })

  output$anatomy <- renderPlotly({
    d <- experiment()
    p <- if (input$fdr) d$p_adjusted else d$p_nominal
    y <- -log10(p)
    sig <- p < alpha
    big <- abs(d$log2fc) >= fc_cut
    # Floors keep the cutoff lines and labels visible in any simulation
    x_max <- max(abs(d$log2fc), 2 * fc_cut) * 1.05
    y_max <- max(y, 3) * 1.2
    cut <- -log10(alpha)

    fc_line <- function(x) {
      list(
        type = "line",
        yref = "paper",
        x0 = x,
        x1 = x,
        y0 = 0,
        y1 = 1,
        line = list(dash = "dash", width = 1.5, color = "#444444")
      )
    }

    label <- function(x, y, text, color = "#444444") {
      list(
        x = x,
        y = y,
        text = text,
        showarrow = FALSE,
        font = list(size = 12, color = color),
        bgcolor = "rgba(255, 255, 255, 0.8)"
      )
    }

    plot_ly(
      x = d$log2fc,
      y = y,
      type = "scatter",
      mode = "markers",
      hoverinfo = "text",
      text = sprintf(
        "<b>Proteína %s</b><br>log2FC = %.2f (fold change = %.2f)<br>p = %s, −log10(p) = %.2f",
        d$protein,
        d$log2fc,
        2^d$log2fc,
        fmt_p(p),
        y
      ),
      marker = list(
        size = 6,
        color = ifelse(
          !sig,
          col_never,
          ifelse(!big, col_small, ifelse(d$log2fc > 0, col_red, col_blue))
        )
      )
    ) |>
      layout(
        xaxis = list(
          title = "← menor no tratado &nbsp;&nbsp; log2(Fold Change) &nbsp;&nbsp; maior no tratado →",
          range = c(-x_max, x_max),
          zeroline = FALSE
        ),
        yaxis = list(
          title = if (input$fdr) {
            "-log10(p-valor ajustado) &nbsp; mais evidência →"
          } else {
            "-log10(p-valor) &nbsp; mais evidência →"
          },
          range = c(-0.05, y_max),
          zeroline = FALSE
        ),
        shapes = list(
          list(
            type = "line",
            xref = "paper",
            x0 = 0,
            x1 = 1,
            y0 = cut,
            y1 = cut,
            line = list(dash = "dash", width = 1.5, color = "#444444")
          ),
          fc_line(-fc_cut),
          fc_line(fc_cut)
        ),
        annotations = list(
          label(
            -(fc_cut + x_max) / 2,
            y_max * 0.93,
            "<b>Diminuem no tratado</b><br>passam nos dois cortes",
            "rgb(55, 126, 184)"
          ),
          label(
            (fc_cut + x_max) / 2,
            y_max * 0.93,
            "<b>Aumentam no tratado</b><br>passam nos dois cortes",
            "rgb(228, 26, 28)"
          ),
          label(
            0,
            y_max * 0.93,
            "<b>Efeito pequeno</b><br>estatisticamente<br>significativo,<br>|log2FC| &lt; 1",
            "rgb(170, 110, 0)"
          ),
          label(0, cut / 2, "abaixo do corte: não se distingue do acaso"),
          list(
            x = fc_cut,
            xanchor = "left",
            yref = "paper",
            y = 0,
            yanchor = "bottom",
            text = " log2FC = +1 (dobro)",
            showarrow = FALSE,
            font = list(size = 12, color = "#444444")
          ),
          list(
            x = -fc_cut,
            xanchor = "right",
            yref = "paper",
            y = 0,
            yanchor = "bottom",
            text = "log2FC = −1 (metade) ",
            showarrow = FALSE,
            font = list(size = 12, color = "#444444")
          ),
          list(
            xref = "paper",
            x = 1,
            xanchor = "right",
            y = cut,
            yanchor = "bottom",
            text = "corte: p = 0,05, ou −log10(p) = 1,30",
            showarrow = FALSE,
            font = list(size = 12, color = "#444444")
          )
        )
      ) |>
      config(displayModeBar = FALSE)
  })

  observeEvent(input$fdr, ignoreInit = TRUE, {
    s <- view_state(experiment(), input$fdr)
    plotlyProxy("volcano", session) |>
      plotlyProxyInvoke(
        "animate",
        list(
          data = list(
            list(
              y = s$y,
              text = s$text,
              marker = list(color = s$color, line = list(color = s$line))
            ),
            list(name = s$legend_red)
          ),
          traces = list(0, 1),
          layout = list(yaxis = list(title = list(text = s$y_title)))
        ),
        list(
          transition = list(duration = 1500, easing = "cubic-in-out"),
          # redraw applies the non-animatable changes (hover text, titles, legend)
          frame = list(duration = 1500, redraw = TRUE)
        )
      )
  })
}

shinyApp(ui, server)
