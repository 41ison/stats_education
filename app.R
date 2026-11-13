library(shiny)
library(bslib)
library(plotly)

alpha <- 0.05

col_red <- "rgba(228, 26, 28, 0.85)"
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
          sprintf("Não significativa (p-nominal = %s).", fmt_p(d$p_nominal))
        )
      ),
      y_title = "-log10(p-valor nominal)",
      legend_red = "Candidatos Iniciais (p nominal < 0,05)"
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
entre as proteínas declaradas significativas, espera-se que no máximo 5%%
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
      "Visão nominal: a ilusão"
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
      legend_trace("Não significativos", col_never) |>
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
