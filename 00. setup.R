# 1. Preparar pacotes ----
# Execute este arquivo uma vez antes da aula. O roteiro principal também o carrega.
# As funções ficam aqui para que a aula acompanhe os passos da análise.
# Documentação: blocos #' no formato roxygen2, com @param, @return e @examples.
# São comentários de consulta no código; não é necessário instalar roxygen2.
# Exemplos com \dontrun usam os objetos preparados pelo roteiro da aula.
dir.create("biblioteca_R", showWarnings = FALSE)
versao_r = paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
lib_usuario = file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", versao_r)
.libPaths(c("biblioteca_R", lib_usuario, .libPaths()))
pacotes = c("igraph", "ggplot2", "ggrepel", "ggpubr", "vegan", "bipartite", "jsonlite")
faltantes = pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltantes)) install.packages(faltantes, lib = "biblioteca_R", repos = "https://cloud.r-project.org")
if (!requireNamespace("ggsankey", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = "biblioteca_R", repos = "https://cloud.r-project.org")
  remotes::install_github("davidsjoberg/ggsankey", lib = "biblioteca_R", upgrade = "never")
}
library(igraph)
library(ggplot2)
library(ggsankey)
# 2. Conferir matriz e categorias ----
#' Conferir a matriz bipartida e as categorias das espécies
#'
#' Interrompe a análise quando os dados violam as condições de entrada da aula.
#' @param matriz Matriz numérica de inteiros não negativos, com pelo menos duas
#'   linhas e duas colunas, nomes únicos entre ambas as guildas e margens positivas.
#' @param categorias Data frame com id, especie, guilda, grupo e cor. Os ids devem
#'   corresponder aos nomes da matriz; cada grupo deve ter uma única cor hexadecimal.
#' @param tipo Texto: "contagem" para frequências ou "binaria" para valores 0/1.
#' @return TRUE invisível quando as verificações passam; caso contrário, um erro.
#' @details Não corrige nem remove registros. Cada guilda deve corresponder a um
#'   lado da matriz para os gráficos. Confira também a identificação taxonômica.
#' @examples
#' \dontrun{
#' conferir_dados(matriz, categorias, tipo = "contagem")
#' }
conferir_dados = function(matriz, categorias, tipo) {
  stopifnot(is.matrix(matriz), is.numeric(matriz), all(is.finite(matriz)),
            all(matriz >= 0), all(matriz == round(matriz)),
            nrow(matriz) >= 2, ncol(matriz) >= 2,
            all(rowSums(matriz) > 0), all(colSums(matriz) > 0))
  stopifnot(tipo %in% c("contagem", "binaria"))
  if (tipo == "binaria") stopifnot(all(matriz %in% c(0, 1)))
  ids = c(rownames(matriz), colnames(matriz))
  stopifnot(!anyDuplicated(ids), !anyDuplicated(categorias$id),
            setequal(ids, categorias$id), all(nzchar(categorias$grupo)),
            all(grepl("^#[0-9A-Fa-f]{6}$", categorias$cor)))
  mapa = unique(categorias[, c("grupo", "cor")])
  stopifnot(!anyDuplicated(mapa$grupo))
  invisible(TRUE)
}
# 3. Preparar os pares de interação ----
#' Converter a matriz em uma lista de pares interagentes
#'
#' Mantém uma linha para cada célula com valor maior que zero.
#' @param matriz Matriz de interações com nomes nas linhas e colunas.
#' @return Data frame com Plantas (linhas), Animais (colunas) e frequencia (peso).
#' @details Os nomes Plantas e Animais são rótulos de colunas; a função não infere
#'   taxonomia. Zeros não são exportados como arestas.
#' @examples
#' m = matrix(c(2, 0, 1, 3), 2, dimnames = list(c("P1", "P2"), c("A1", "A2")))
#' listar_interacoes(m)
listar_interacoes = function(matriz) {
  pos = which(matriz > 0, arr.ind = TRUE)
  data.frame(Plantas = rownames(matriz)[pos[, 1]], Animais = colnames(matriz)[pos[, 2]],
             frequencia = as.numeric(matriz[pos]))
}
# 4. Desenhar o Sankey com categorias e nomes alinhados ----
# cores recebe um vetor nomeado, por exemplo: c(Plantas = "#009E73").
# Categorias não informadas mantêm suas cores; nomes de outras redes são ignorados.
#' Substituir cores de uma paleta nomeada
#'
#' Atualiza somente os grupos informados e mantém os demais valores da paleta.
#' @param padrao Vetor de cores nomeado pelas categorias existentes.
#' @param cores Vetor de cores nomeado pelas categorias a alterar, ou NULL para
#'   manter o padrão. Aceita nomes de cores do R e códigos hexadecimais.
#' @return Vetor nomeado com as cores atualizadas.
#' @details Nomes de grupos ausentes são ignorados, mas ocorre erro se nenhum nome
#'   informado corresponder à paleta. Cores inválidas também produzem erro.
#' @examples
#' personalizar_cores(c(Plantas = "green", Animais = "orange"), c(Plantas = "#009E73"))
personalizar_cores = function(padrao, cores = NULL) {
  if (is.null(cores)) return(padrao)
  if (!is.character(cores) || is.null(names(cores)) || anyNA(cores) ||
      anyNA(names(cores)) || any(!nzchar(names(cores))) || anyDuplicated(names(cores))) {
    stop("Informe cores como vetor nomeado, por exemplo: c(Plantas = '#009E73').", call. = FALSE)
  }
  grDevices::col2rgb(cores)
  usados = intersect(names(padrao), names(cores))
  if (!length(usados)) stop("Nenhum nome da paleta corresponde a: ", paste(names(padrao), collapse = ", "), call. = FALSE)
  padrao[usados] = cores[usados]
  padrao
}
# Limita somente os nomes, priorizando a soma dos pesos; empates seguem a ordem dos dados.
#' Selecionar os nós que receberão nomes
#'
#' Ordena os candidatos pela soma dos pesos e limita somente a quantidade de textos.
#' @param pesos Vetor numérico de somas dos pesos por nó.
#' @param min_rotulo Número não negativo: soma mínima para exibir o nome; padrão 10.
#' @param n_rotulos Inteiro não negativo: máximo de nomes; padrão 20. Inf permite
#'   todos os elegíveis e 0 oculta todos. Empates seguem a ordem de entrada.
#' @return Vetor de índices dos nós selecionados, em ordem decrescente de peso.
#' @examples
#' selecionar_rotulos(c(4, 20, 10, 20), min_rotulo = 0, n_rotulos = 2)
selecionar_rotulos = function(pesos, min_rotulo = 10, n_rotulos = 20) {
  stopifnot(length(min_rotulo) == 1, is.finite(min_rotulo), min_rotulo >= 0,
            is.numeric(n_rotulos), length(n_rotulos) == 1, !is.na(n_rotulos), n_rotulos >= 0,
            n_rotulos == Inf || n_rotulos == floor(n_rotulos))
  candidatos = which(pesos >= min_rotulo)
  candidatos = candidatos[order(-pesos[candidatos], candidatos)]
  head(candidatos, min(n_rotulos, length(candidatos)))
}
#' Desenhar o Sankey completo da rede de interação
#'
#' Representa os pesos como larguras e centraliza cada nome no respectivo nó.
#' @param matriz Matriz de interações validada por conferir_dados().
#' @param categorias Data frame com id, especie, guilda, grupo e cor.
#' @param info Metadados da análise. Argumento mantido por compatibilidade; o
#'   desenho atual obtém as guildas e categorias de categorias, sem consultar info.
#' @param min_rotulo Soma mínima de pesos para um nome ser elegível; padrão 10.
#' @param cores Paleta nomeada por grupo, ou NULL para usar categorias$cor.
#' @param cor_fluxos Cor única dos fluxos, ou NULL para a cor da categoria de origem.
#' @param cor_borda Cor do contorno dos nós; NA, padrão, remove o contorno.
#' @param cor_texto Cor dos títulos, eixos e legendas; padrão "black".
#' @param n_rotulos Máximo de nomes, priorizando maiores somas; padrão 20.
#'   Inf inclui todos os elegíveis e 0 oculta os nomes.
#' @return Objeto ggplot, sem exportação automática.
#' @details Todas as interações são mantidas. Os textos são coloridos por grupo e
#'   escritos em itálico; nomes sobrepostos podem ser omitidos por check_overlap.
#' @examples
#' \dontrun{
#' desenhar_rede(matriz, categorias, info, n_rotulos = 15, cores = c(Plantas = "darkgreen"))
#' }
desenhar_rede = function(matriz, categorias, info, min_rotulo = 10, cores = NULL,
                         cor_fluxos = NULL, cor_borda = NA, cor_texto = "black", n_rotulos = 20) {
  stopifnot(length(min_rotulo) == 1, is.finite(min_rotulo), min_rotulo >= 0)
  pares = listar_interacoes(matriz)
  plot_data = make_long(pares, Plantas, Animais, value = frequencia)
  plot_data$grupo = factor(categorias$grupo[match(as.character(plot_data$node), categorias$id)],
                           levels = unique(categorias$grupo))
  paleta = unique(categorias[, c("grupo", "cor")])
  cores = personalizar_cores(setNames(paleta$cor, paleta$grupo), cores)
  guilda_linhas = unique(categorias$guilda[match(rownames(matriz), categorias$id)])
  guilda_colunas = unique(categorias$guilda[match(colnames(matriz), categorias$id)])
  stopifnot(length(guilda_linhas) == 1, length(guilda_colunas) == 1)
  parametros_sankey = list(flow.alpha = 0.55, flow.color = NA, node.color = cor_borda, width = 0.12, space = 2)
  # NULL mantém os fluxos com a cor da categoria de origem; uma cor fixa os uniformiza.
  if (!is.null(cor_fluxos)) parametros_sankey$flow.fill = cor_fluxos
  base = ggplot(plot_data, aes(x = x, next_x = next_x, node = node, next_node = next_node, value = value, fill = grupo)) +
    do.call(geom_sankey, parametros_sankey) +
    scale_fill_manual(values = cores, drop = FALSE) +
    scale_x_discrete(labels = c(guilda_linhas, guilda_colunas), expand = expansion(add = c(0.7, 0.7))) +
    labs(x = NULL, y = NULL, fill = "Categoria") +
    theme_minimal(base_family = "sans") +
    theme(text = element_text(color = cor_texto), panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(size = 14, color = cor_texto, face = "bold"),
          legend.position = "right", legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10), legend.background = element_blank(),
          plot.title = element_text(size = 16, face = "bold"), plot.subtitle = element_text(size = 11),
          plot.caption = element_text(size = 10, hjust = 0), plot.margin = margin(15, 20, 15, 20))
  construido = ggplot_build(base)
  camada_nos = which(vapply(base$layers, function(z) inherits(z$geom, "GeomRect"), logical(1)))
  camada_fluxos = which(vapply(base$layers, function(z) inherits(z$geom, "GeomPolygon"), logical(1)))
  stopifnot(length(camada_nos) == 1, length(camada_fluxos) == 1)
  nos = construido$data[[camada_nos]]
  fluxos = construido$data[[camada_fluxos]]
  stopifnot(nrow(nos) == sum(dim(matriz)), length(unique(fluxos$group)) == sum(matriz > 0))
  pesos = c(rowSums(matriz), colSums(matriz))
  stopifnot(all(nos$freq == pesos[match(as.character(nos$node), names(pesos))]))
  # Filtrar apenas nomes. Os retângulos e fluxos permanecem completos.
  rotulos = nos[selecionar_rotulos(nos$freq, min_rotulo, n_rotulos), ]
  rotulos$nome = categorias$especie[match(as.character(rotulos$node), categorias$id)]
  rotulos$categoria = categorias$grupo[match(as.character(rotulos$node), categorias$id)]
  rotulos$y_texto = (rotulos$ymin + rotulos$ymax) / 2
  esquerda = rotulos$x == min(nos$x)
  rotulos$x_texto = ifelse(esquerda, rotulos$xmin - 0.02, rotulos$xmax + 0.02)
  rotulos$alinhamento = ifelse(esquerda, 1, 0)
  grafico = base +
    geom_text(data = rotulos, aes(x = x_texto, y = y_texto, label = nome, hjust = alinhamento, color = categoria),
              inherit.aes = FALSE, vjust = 0.5, family = "sans", fontface = "italic", size = 3, check_overlap = TRUE) +
    scale_color_manual(values = cores, guide = "none")
  grafico
}
# 4.1. Grafo do igraph desenhado com ggplot2 ----
# Padrão: livre (Fruchterman–Reingold), com todos os nós distribuídos pelo plano.
# A estrutura continua bipartida: o layout não cria interações dentro das guildas.
# O layout muda somente a posição dos nós; arestas, categorias e pesos são mantidos.
# Referências: https://r.igraph.org/reference/layout_with_fr.html
# https://ggrepel.slowkow.com/reference/geom_text_repel.html
#' Desenhar um grafo do igraph com ggplot2 e nomes repelidos
#'
#' Mostra todas as espécies e interações, com cores por categoria e área dos nós
#' proporcional ao grau. Os nomes ficam em itálico e têm linhas-guia tracejadas.
#' @param grafo Grafo igraph não direcionado, com atributos name e type nos
#'   vértices e weight positivo nas arestas; use graph_from_biadjacency_matrix().
#' @param categorias Data frame com id, especie, guilda, grupo e cor para os nós.
#' @param min_rotulo Soma mínima dos pesos para selecionar um nome; padrão 10.
#' @param disposicao "livre" para Fruchterman–Reingold (padrão) ou "bipartido"
#'   para duas colunas. A estrutura de interações permanece a mesma.
#' @param semente Inteiro usado no layout e no ggrepel; padrão 123.
#' @param curvatura Número entre -1 e 1; padrão 0.12. Zero gera arestas retas.
#'   Na disposição bipartida, as arestas são sempre retas.
#' @param ponderar_layout Lógico; FALSE atribui atração igual às arestas. TRUE
#'   usa os pesos da interação para a atração no layout livre.
#' @param cores Vetor de cores nomeado por grupo, ou NULL para categorias$cor.
#' @param cor_arestas Cor das ligações ecológicas; padrão "#607D8B".
#' @param cor_borda Cor do contorno dos nós; padrão "grey20".
#' @param cor_texto Cor das legendas e títulos; padrão "black". Os nomes das
#'   espécies recebem a cor da categoria, independentemente desse argumento.
#' @param cor_guias Cor das linhas entre nomes e nós; padrão "grey55".
#' @param posicao_rotulos "proxima" (padrão) desloca nomes localmente; "externa"
#'   usa as margens. "livre" é um sinônimo de "proxima". Só afeta o layout livre.
#' @param n_rotulos Máximo de nomes; padrão 20. Inf inclui todos os elegíveis e
#'   0 oculta os nomes. A prioridade é a soma dos pesos, não o grau.
#' @param atracao_rotulos Número positivo, padrão 0.5; maior atração aproxima
#'   nomes dos nós no modo "proxima".
#' @param repulsao_rotulos Número positivo, padrão 2; maior repulsão abre espaço
#'   entre os textos e os nós no modo "proxima".
#' @return Objeto ggplot, sem exportação automática.
#' @details O layout não representa distância ecológica ou geográfica. A função
#'   define a semente global do R. O ggrepel reduz sobreposições, mas a leitura
#'   também depende do número de nomes e das dimensões físicas da figura.
#' @seealso https://r.igraph.org/reference/layout_with_fr.html
#'   https://ggrepel.slowkow.com/reference/geom_text_repel.html
#' @examples
#' \dontrun{
#' g = igraph::graph_from_biadjacency_matrix(matriz, weighted = TRUE)
#' desenhar_grafo_ggplot(g, categorias, min_rotulo = 0, n_rotulos = 12)
#' }
desenhar_grafo_ggplot = function(grafo, categorias, min_rotulo = 10, disposicao = "livre",
                                semente = 123, curvatura = 0.12, ponderar_layout = FALSE,
                                cores = NULL, cor_arestas = "#607D8B", cor_borda = "grey20",
                                cor_texto = "black", cor_guias = "grey55", posicao_rotulos = "proxima",
                                n_rotulos = 20, atracao_rotulos = 0.5, repulsao_rotulos = 2) {
  stopifnot(inherits(grafo, "igraph"), !igraph::is_directed(grafo),
            disposicao %in% c("bipartido", "livre"),
            posicao_rotulos %in% c("proxima", "externa", "livre"),
            length(atracao_rotulos) == 1, is.finite(atracao_rotulos), atracao_rotulos > 0,
            length(repulsao_rotulos) == 1, is.finite(repulsao_rotulos), repulsao_rotulos > 0,
            length(min_rotulo) == 1, is.finite(min_rotulo), min_rotulo >= 0,
            length(curvatura) == 1, is.finite(curvatura), abs(curvatura) <= 1,
            is.logical(ponderar_layout), length(ponderar_layout) == 1, !is.na(ponderar_layout))
  ids = igraph::V(grafo)$name
  tipos = igraph::V(grafo)$type
  stopifnot(length(tipos) == length(ids), is.logical(tipos), any(tipos), any(!tipos),
            !anyDuplicated(ids), !anyDuplicated(categorias$id), all(ids %in% categorias$id))
  nos = categorias[match(ids, categorias$id), ]
  nos$grupo = factor(nos$grupo, levels = unique(categorias$grupo))
  nos$grau = as.numeric(igraph::degree(grafo))
  nos$forca = as.numeric(igraph::strength(grafo, weights = igraph::E(grafo)$weight))
  # Nas duas colunas, usar a ordem calculada pelo igraph para reduzir cruzamentos.
  if (disposicao == "bipartido") {
    posicoes = igraph::layout_as_bipartite(grafo, types = tipos)
    nos$x = as.numeric(tipos)
    nos$y = 0
    for (tipo in c(FALSE, TRUE)) {
      indices = which(tipos == tipo)
      ordem = indices[order(posicoes[indices, 1])]
      nos$y[ordem] = seq(1, -1, length.out = length(ordem))
    }
    nos$x_rotulo = nos$x + ifelse(tipos, 0.04, -0.04)
    nos$y_rotulo = nos$y
    nos$alinhamento = ifelse(tipos, 0, 1)
  } else {
    set.seed(semente)
    # Por padrão, as posições dependem dos parceiros; os pesos aparecem nas linhas.
    # TRUE aproxima mais os pares de maior peso. Isso altera apenas o desenho.
    pesos_layout = if (ponderar_layout) igraph::E(grafo)$weight else rep(1, igraph::ecount(grafo))
    posicoes = igraph::layout_with_fr(grafo, weights = pesos_layout, niter = 2000, grid = "nogrid")
    # Centralizar e escalar por um único fator preserva as proporções do layout.
    posicoes = sweep(posicoes, 2, apply(posicoes, 2, function(x) mean(range(x))), "-")
    posicoes = posicoes / max(abs(posicoes))
    nos$x = posicoes[, 1]
    nos$y = posicoes[, 2]
    nos$x_rotulo = nos$x
    nos$y_rotulo = nos$y + 0.055
    nos$alinhamento = 0.5
  }
  arestas = igraph::as_data_frame(grafo, what = "edges")
  stopifnot(nrow(arestas) > 0, all(is.finite(arestas$weight)), all(arestas$weight > 0))
  origem = match(arestas$from, nos$id)
  destino = match(arestas$to, nos$id)
  arestas$x = nos$x[origem]
  arestas$y = nos$y[origem]
  arestas$xend = nos$x[destino]
  arestas$yend = nos$y[destino]
  paleta = unique(categorias[categorias$id %in% ids, c("grupo", "cor")])
  stopifnot(!anyDuplicated(paleta$grupo))
  cores = personalizar_cores(setNames(paleta$cor, paleta$grupo), cores)
  indices_rotulos = selecionar_rotulos(nos$forca, min_rotulo, n_rotulos)
  rotulos = nos[indices_rotulos, ]
  nos$rotulo = ""
  nos$rotulo[indices_rotulos] = as.character(nos$especie[indices_rotulos])
  curvatura_usada = if (disposicao == "bipartido") 0 else curvatura
  figura = ggplot(nos, aes(x = x, y = y))
  if (length(unique(arestas$weight)) == 1) {
    figura = figura +
      geom_curve(data = arestas, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, color = cor_arestas, alpha = 0.4, linewidth = 0.4,
                 curvature = curvatura_usada, lineend = "round")
  } else {
    figura = figura +
      geom_curve(data = arestas, aes(x = x, y = y, xend = xend, yend = yend, linewidth = weight),
                 inherit.aes = FALSE, color = cor_arestas, alpha = 0.4,
                 curvature = curvatura_usada, lineend = "round") +
      scale_linewidth_continuous(name = "Peso da interação", range = c(0.2, 1.6))
  }
  figura = figura +
    geom_point(aes(fill = grupo, size = grau), shape = 21, color = cor_borda, stroke = 0.3) +
    scale_fill_manual(name = "Categoria", values = cores) +
    scale_color_manual(values = cores, guide = "none") +
    scale_size_area(name = "Número de parceiros", max_size = 9, aesthetics = c("size", "point.size")) +
    guides(fill = guide_legend(order = 1, override.aes = list(size = 4)),
           size = guide_legend(order = 2, override.aes = list(fill = "grey45", color = "grey45")),
           linewidth = guide_legend(order = 3)) +
    theme_void(base_family = "sans") +
    theme(text = element_text(color = cor_texto), legend.position = "right", legend.title = element_text(size = 12, face = "bold", color = cor_texto),
          legend.text = element_text(size = 10, color = cor_texto), legend.background = element_blank(),
          plot.margin = margin(15, 20, 15, 20))
  if (disposicao == "bipartido") {
    guildas = c(as.character(nos$guilda[which(!tipos)[1]]), as.character(nos$guilda[which(tipos)[1]]))
    figura = figura +
      geom_text(data = rotulos, aes(x = x_rotulo, y = y_rotulo, label = especie, hjust = alinhamento, color = grupo),
                inherit.aes = FALSE, family = "sans", fontface = "italic", size = 3, check_overlap = TRUE) +
      annotate("text", x = c(0, 1), y = -1.12, label = guildas, family = "sans", size = 4.5, fontface = "bold", color = cor_texto) +
      coord_cartesian(xlim = c(-0.6, 1.6), ylim = c(-1.18, 1.05), expand = FALSE, clip = "off")
  } else if (posicao_rotulos == "externa") {
    # A posição dos nós permanece a do Fruchterman–Reingold.
    # Somente os textos vão para as margens; os lados não representam guildas.
    rotulos$lado = ifelse(rotulos$x < 0, -1, 1)
    for (lado in c(-1, 1)) {
      texto_lado = rotulos[rotulos$lado == lado, ]
      if (!nrow(texto_lado)) next
      texto_lado = texto_lado[order(texto_lado$y), ]
      altura_textos = if (nrow(texto_lado) == 1) texto_lado$y else seq(-1.05, 1.05, length.out = nrow(texto_lado))
      # position_nudge_repel conserva a origem da linha no nó correto.
      figura = figura +
        ggrepel::geom_text_repel(data = texto_lado, aes(x = x, y = y, label = especie, point.size = grau, color = grupo),
                                inherit.aes = FALSE, family = "sans", fontface = "italic", size = 3.2,
                                hjust = if (lado == -1) 1 else 0, direction = "y",
                                position = ggrepel::position_nudge_repel(x = lado * 1.2 - texto_lado$x,
                                                                        y = altura_textos - texto_lado$y),
                                seed = semente, box.padding = 0.35, point.padding = 0.15,
                                min.segment.length = 0, segment.color = cor_guias, segment.size = 0.3,
                                segment.linetype = "dashed", force = 1, force_pull = 1, max.overlaps = Inf,
                                ylim = c(-1.15, 1.15), max.time = 2, max.iter = 100000, show.legend = FALSE)
    }
    figura = figura + coord_fixed(xlim = c(-2.15, 2.15), ylim = c(-1.2, 1.2), clip = "off")
  } else {
    # Todos os nós funcionam como obstáculos, mesmo quando seu nome não é exibido.
    # A atração mantém cada nome perto da espécie; a repulsão abre espaço localmente.
    # O halo branco melhora a leitura nos cruzamentos com as linhas da rede.
    angulos = seq_len(nrow(nos)) * pi * (3 - sqrt(5))
    figura = figura +
      ggrepel::geom_text_repel(aes(label = rotulo, point.size = grau, color = grupo), family = "sans", fontface = "italic", size = 3.2,
                              seed = semente, box.padding = 0.2, point.padding = 0.1,
                              position = ggrepel::position_nudge_repel(x = 0.05 * cos(angulos), y = 0.05 * sin(angulos)),
                              min.segment.length = 0, segment.color = cor_guias, segment.size = 0.25, segment.linetype = "dashed",
                              force = repulsao_rotulos, force_pull = atracao_rotulos, max.overlaps = Inf,
                              bg.color = "white", bg.r = 0.08,
                              max.time = 5, max.iter = 1000000, show.legend = FALSE) +
      coord_fixed(xlim = c(-1.25, 1.25), ylim = c(-1.2, 1.2), clip = "off")
  }
  figura
}
# 4.2. Exportar a adaptação do grafo ----
#' Exportar o grafo em PNG e TIFF
#'
#' Salva versões de apresentação e impressão da mesma figura.
#' @param figura Objeto ggplot produzido pela função de desenho.
#' @param pasta Diretório de destino; padrão "figures". Criado se necessário.
#' @param nome Nome-base sem extensão; padrão "grafo_livre".
#' @return A figura, invisivelmente. Grava PNG a 180 dpi e TIFF/LZW a 600 dpi,
#'   ambos com 14 × 10 polegadas e fundo branco. Arquivos homônimos são substituídos.
#' @examples
#' \dontrun{
#' salvar_grafo_ggplot(figura_grafo, pasta = pasta_figuras)
#' }
salvar_grafo_ggplot = function(figura, pasta = "figures", nome = "grafo_livre") {
  dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(pasta, paste0(nome, ".png")), plot = figura, width = 14, height = 10, units = "in", dpi = 180, bg = "white")
  ggsave(file.path(pasta, paste0(nome, ".tiff")), plot = figura, width = 14, height = 10, units = "in", dpi = 600, compression = "lzw", bg = "white")
  invisible(figura)
}
# 5. Detectar módulos: mesmo esforço no observado e nas simulações ----
#' Procurar a partição de maior modularidade entre buscas Louvain
#'
#' Repete a busca estocástica e mantém a partição com maior Q ponderado.
#' @param grafo Grafo igraph não direcionado com pesos positivos em weight.
#' @param n_inicios Inteiro positivo de buscas independentes; padrão 10.
#' @param semente Inteiro para reproduzir as buscas; padrão 123.
#' @return Vetor de pertencimento aos módulos, na ordem dos vértices do grafo.
#' @details Define a semente global. O melhor resultado é aproximado; Q é a
#'   modularidade geral do igraph, não uma formulação específica bipartida.
#' @examples
#' \dontrun{
#' modulos = detectar_modulos(grafo, n_inicios = 10)
#' igraph::modularity(grafo, modulos, weights = igraph::E(grafo)$weight)
#' }
detectar_modulos = function(grafo, n_inicios = 10, semente = 123) {
  stopifnot(n_inicios >= 1, n_inicios == as.integer(n_inicios))
  set.seed(semente)
  grupos = replicate(n_inicios, membership(cluster_louvain(grafo, weights = E(grafo)$weight)), simplify = FALSE)
  q = vapply(grupos, function(x) modularity(grafo, x, weights = E(grafo)$weight), numeric(1))
  grupos[[which.max(q)]]
}
# 6. Repetir as métricas do roteiro nas redes nulas ----
#' Calcular as cinco métricas de uma matriz de interação
#'
#' Calcula conectância, NODF binário, Q ponderado, ligações por espécie e H2'.
#' @param matriz Matriz numérica de interações com margens positivas.
#' @param n_inicios Número de buscas Louvain para Q; padrão 10.
#' @param semente Semente da detecção dos módulos; padrão 123.
#' @param retornar_modulos Lógico; FALSE (padrão) retorna somente as métricas.
#'   TRUE também retorna a partição Louvain usada no cálculo de Q.
#' @return Vetor numérico nomeado: conectancia, NODF, Q, ligacoes_por_especie e H2.
#'   Com retornar_modulos = TRUE, retorna uma lista com metricas e modulos.
#' @details C = L/(P*A) e L/S = L/(P+A) usam o número de células positivas.
#'   NODF usa presença/ausência e varia de 0 a 100. Q e H2' usam os pesos.
#'   Em dados binários, H2' é descritivo e não estima preferência por frequência.
#'   A função usa a semente global por meio de detectar_modulos().
#' @examples
#' \dontrun{
#' observado = medir_rede(matriz)
#' }
medir_rede = function(matriz, n_inicios = 10, semente = 123, retornar_modulos = FALSE) {
  stopifnot(is.logical(retornar_modulos), length(retornar_modulos) == 1, !is.na(retornar_modulos))
  binaria = 1L * (matriz > 0)
  grafo = graph_from_biadjacency_matrix(matriz, directed = FALSE, weighted = TRUE)
  modulos = detectar_modulos(grafo, n_inicios, semente)
  valores = c(conectancia = sum(binaria) / length(binaria),
    NODF = unname(vegan::nestednodf(binaria, order = TRUE, weighted = FALSE)$statistic["NODF"]),
    Q = modularity(grafo, modulos, weights = E(grafo)$weight),
    ligacoes_por_especie = sum(binaria) / sum(dim(binaria)),
    H2 = as.numeric(bipartite::networklevel(matriz, index = "H2", H2_integer = TRUE)))
  if (retornar_modulos) return(list(metricas = valores, modulos = modulos))
  valores
}
# 7. Gerar modelos nulos adequados ao tipo de dado ----
# Contagens: swap.web para NODF/Q/H2'; r2dtable para conectância e L/S.
# Binária: quasiswap mantém graus e valores 0/1. C, L/S e H2' ficam fixados.
# Não transformar presença/ausência em contagens pelo uso de r2dtable.
#' Gerar redes nulas com restrições adequadas ao tipo de matriz
#'
#' Usa swap.web e r2dtable para contagens, ou quasiswap para presença/ausência.
#' @param matriz Matriz de inteiros não negativos, validada por conferir_dados().
#' @param tipo "contagem" ou "binaria"; deve corresponder aos dados observados.
#' @param B Número inteiro de redes por modelo, no mínimo 19; padrão 999.
#' @param semente Inteiro usado nos sorteios; padrão 1102. r2dtable usa semente + 1.
#' @return Lista nomeada por modelo, contendo B matrizes completas em cada elemento.
#' @details Todos os modelos mantêm margens. swap.web também mantém L; quasiswap
#'   mantém L e entradas 0/1. r2dtable deixa L variar e não é usado em dados binários.
#'   Verifica essas restrições após os sorteios e altera o estado aleatório global.
#'   Não impõe zeros estruturais. A inferência depende das restrições do modelo.
#' @examples
#' \dontrun{
#' nulas = gerar_nulas(matriz, tipo = info$tipo, B = 199)
#' }
gerar_nulas = function(matriz, tipo, B = 999, semente = 1102) {
  stopifnot(B >= 19, B == as.integer(B), tipo %in% c("contagem", "binaria"))
  message("Gerando ", B, " redes por modelo (", tipo, ")...")
  set.seed(semente)
  if (tipo == "contagem") {
    nulas = list(swap.web = bipartite::swap.web(N = B, web = matriz, verbose = FALSE))
    set.seed(semente + 1)
    nulas$r2dtable = stats::r2dtable(n = B, r = rowSums(matriz), c = colSums(matriz))
  } else {
    objeto = vegan::nullmodel(matriz, method = "quasiswap")
    cubo = stats::simulate(objeto, nsim = B, seed = semente)
    nulas = list(quasiswap = lapply(seq_len(B), function(i) cubo[, , i]))
  }
  for (modelo in names(nulas)) {
    for (m in nulas[[modelo]]) {
      stopifnot(identical(dim(m), dim(matriz)), all(is.finite(m)), all(m >= 0), all(m == round(m)),
                all(rowSums(m) == rowSums(matriz)), all(colSums(m) == colSums(matriz)))
      if (modelo != "r2dtable") stopifnot(sum(m > 0) == sum(matriz > 0))
      if (modelo == "quasiswap") stopifnot(all(m %in% c(0, 1)))
    }
  }
  nulas
}
# 8. Construir as distribuições nulas ----
#' Calcular as métricas nas redes nulas
#'
#' Usa o mesmo esforço de busca de modularidade em todas as réplicas.
#' @param nulas Lista de matrizes por modelo, produzida por gerar_nulas().
#' @param n_inicios Número de buscas Louvain em cada réplica; padrão 10.
#' @return Lista de matrizes numéricas; linhas são réplicas e colunas são métricas.
#'   Para r2dtable, calcula somente conectancia e ligacoes_por_especie.
#' @details Para os demais modelos, calcula as cinco métricas; testar_metricas()
#'   remove as constantes. A modularidade da réplica i usa a semente 40000 + i.
#' @examples
#' \dontrun{
#' simulacoes = medir_nulas(nulas, n_inicios = 10)
#' }
medir_nulas = function(nulas, n_inicios = 10) {
  simulacoes = list()
  for (modelo in names(nulas)) {
    lista = nulas[[modelo]]
    if (modelo == "r2dtable") {
      L = vapply(lista, function(m) sum(m > 0), numeric(1))
      simulacoes[[modelo]] = cbind(conectancia = L / length(lista[[1]]),
                                   ligacoes_por_especie = L / sum(dim(lista[[1]])))
    } else {
      valores = matrix(NA_real_, length(lista), 5,
                       dimnames = list(NULL, c("conectancia", "NODF", "Q", "ligacoes_por_especie", "H2")))
      for (i in seq_along(lista)) {
        valores[i, ] = medir_rede(lista[[i]], n_inicios, semente = 40000 + i)
        if (i %% 100 == 0) message(modelo, ": ", i, "/", length(lista))
      }
      simulacoes[[modelo]] = valores
    }
  }
  simulacoes
}
# 9. SES, p empírico bilateral e correção de Holm ----
nomes_metricas = c(conectancia = "Conectância", NODF = "Aninhamento (NODF)",
                   Q = "Modularidade (Q)", ligacoes_por_especie = "Ligações por espécie (L/S)",
                   H2 = "Especialização (H2')")
#' Comparar o observado com as distribuições nulas
#'
#' Calcula SES, caudas empíricas, p bilateral e ajuste de Holm.
#' @param observado Vetor nomeado de métricas, como o retorno de medir_rede().
#' @param simulacoes Lista de distribuições retornada por medir_nulas().
#' @return Data frame com modelo, código/nome da métrica, observado, média/desvio
#'   nulos, quantis, SES, p_inferior, p_superior, p_bilateral, B, status e p_holm.
#'   Métricas fixadas pelo modelo ou sem variação simulada são omitidas.
#' @details SES = (observado - média nula)/desvio nulo. Cada cauda usa
#'   (1 + número de extremos)/(B + 1), incluindo empates; o p bilateral é o dobro
#'   da menor cauda, limitado a 1. Holm é aplicado aos testes não redundantes;
#'   conectância e L/S compartilham o ajuste. Quantis nulos não são ICs do observado.
#' @examples
#' \dontrun{
#' testes = testar_metricas(observado, simulacoes)
#' testes[, c("metrica", "SES", "p_bilateral", "p_holm")]
#' }
testar_metricas = function(observado, simulacoes) {
  tabelas = list()
  for (modelo in names(simulacoes)) {
    sim = simulacoes[[modelo]]
    tabelas[[modelo]] = do.call(rbind, lapply(colnames(sim), function(codigo) {
      x = sim[, codigo]
      obs = observado[[codigo]]
      stopifnot(all(is.finite(x)), is.finite(obs))
      fixa = (modelo != "r2dtable" && codigo %in% c("conectancia", "ligacoes_por_especie")) ||
        (modelo == "quasiswap" && codigo == "H2")
      if (fixa) stopifnot(all(abs(x - obs) < 1e-8))
      desvio = sd(x)
      variavel = !fixa && desvio > 1e-12
      tolerancia = 1e-10 * max(1, abs(obs))
      p_inf = if (variavel) (1 + sum(x <= obs + tolerancia)) / (length(x) + 1) else NA_real_
      p_sup = if (variavel) (1 + sum(x >= obs - tolerancia)) / (length(x) + 1) else NA_real_
      data.frame(modelo = modelo, codigo = codigo, metrica = nomes_metricas[[codigo]],
                 observado = obs, media_nula = mean(x), dp_nulo = if (fixa) 0 else desvio,
                 q025 = unname(quantile(x, 0.025)), q975 = unname(quantile(x, 0.975)),
                 SES = if (variavel) (obs - mean(x)) / desvio else NA_real_,
                 p_inferior = p_inf, p_superior = p_sup,
                 p_bilateral = if (variavel) min(1, 2 * min(p_inf, p_sup)) else NA_real_, B = length(x),
                 status = if (fixa) "Fixada pelo modelo" else if (!variavel) "Sem variação simulada" else "Comparação calculada")
    }))
  }
  resultado = do.call(rbind, tabelas)
  # A tabela final contém apenas métricas com distribuição nula variável.
  resultado = resultado[resultado$status == "Comparação calculada", ]
  rownames(resultado) = NULL
  # C e L/S são proporcionais: contar essa evidência apenas uma vez.
  selecionados = is.finite(resultado$p_bilateral) & resultado$codigo != "ligacoes_por_especie"
  resultado$p_holm = NA_real_
  resultado$p_holm[selecionados] = p.adjust(resultado$p_bilateral[selecionados], method = "holm")
  i_c = which(resultado$modelo == "r2dtable" & resultado$codigo == "conectancia")
  i_ls = which(resultado$modelo == "r2dtable" & resultado$codigo == "ligacoes_por_especie")
  if (length(i_c)) {
    stopifnot(isTRUE(all.equal(resultado$p_bilateral[i_c], resultado$p_bilateral[i_ls])),
              isTRUE(all.equal(resultado$SES[i_c], resultado$SES[i_ls], tolerance = 1e-8)))
    resultado$p_holm[i_ls] = resultado$p_holm[i_c]
  }
  resultado
}
# 10. Visualizar o observado e as simulações ----
#' Organizar as distribuições nulas em formato longo
#'
#' Converte matrizes de simulações em uma tabela para gráficos e exportação.
#' @param simulacoes Lista de matrizes numéricas produzida por medir_nulas().
#' @return Data frame com modelo, codigo, metrica, simulacao e valor.
#' @examples
#' \dontrun{
#' distribuicoes = organizar_distribuicoes(simulacoes)
#' head(distribuicoes)
#' }
organizar_distribuicoes = function(simulacoes) {
  do.call(rbind, lapply(names(simulacoes), function(modelo) {
    sim = simulacoes[[modelo]]
    do.call(rbind, lapply(colnames(sim), function(codigo) {
      data.frame(modelo = modelo, codigo = codigo, metrica = nomes_metricas[[codigo]],
                 simulacao = seq_len(nrow(sim)), valor = sim[, codigo])
    }))
  }))
}
#' Desenhar histogramas nulos com a referência observada
#'
#' Cria um gráfico por modelo, com painéis somente para as métricas testáveis.
#' @param simulacoes Lista de distribuições retornada por medir_nulas().
#' @param testes Data frame retornado por testar_metricas().
#' @param cores Paleta nomeada por swap.web, r2dtable e/ou quasiswap; NULL usa o padrão.
#' @param cor_observado Cor da linha vertical observada; padrão "#D55E00".
#' @param cor_borda Cor do contorno das barras; padrão "white".
#' @param cor_texto Cor dos textos; padrão "black".
#' @return Lista nomeada de objetos ggplot. Os gráficos também são impressos no
#'   dispositivo ativo; não grava arquivos.
#' @examples
#' \dontrun{
#' graficos_testes = desenhar_testes(simulacoes, testes)
#' }
desenhar_testes = function(simulacoes, testes, cores = NULL, cor_observado = "#D55E00",
                           cor_borda = "white", cor_texto = "black") {
  distribuicoes = organizar_distribuicoes(simulacoes)
  cores = personalizar_cores(c(swap.web = "#009E73", r2dtable = "#56B4E9", quasiswap = "#0072B2"), cores)
  graficos = list()
  for (modelo in names(simulacoes)) {
    referencia = testes[testes$modelo == modelo & is.finite(testes$p_bilateral), ]
    dados = distribuicoes[distribuicoes$modelo == modelo & distribuicoes$codigo %in% referencia$codigo, ]
    if (!nrow(dados)) next
    niveis = unname(nomes_metricas[nomes_metricas %in% referencia$metrica])
    dados$metrica = factor(dados$metrica, levels = niveis)
    referencia$metrica = factor(referencia$metrica, levels = niveis)
    grafico = ggplot(dados, aes(x = valor)) +
      geom_histogram(bins = 25, fill = cores[[modelo]], color = cor_borda) +
      geom_vline(data = referencia, aes(xintercept = observado), color = cor_observado, linewidth = 0.9) +
      facet_wrap(~metrica, scales = "free", ncol = min(3, nrow(referencia))) +
      labs(x = "Valor da métrica", y = "Número de simulações", title = paste("Modelo nulo:", modelo),
           subtitle = "Linha vertical: observado | Histogramas: redes aleatorizadas",
           caption = if (modelo == "r2dtable") "Conectância e L/S representam a mesma evidência." else "A comparação depende das restrições do modelo nulo.") +
      theme_minimal(base_family = "sans") +
      theme(text = element_text(color = cor_texto), axis.text = element_text(size = 12, color = cor_texto), axis.title = element_text(size = 14, face = "bold"),
            strip.text = element_text(size = 12, face = "bold"), panel.grid.minor = element_blank(),
            plot.title = element_text(size = 16, face = "bold"), legend.position = "right")
    graficos[[modelo]] = grafico
    print(grafico)
  }
  graficos
}
# 11. Exportar figuras, tabelas e objetos ----
#' Salvar tabelas, figuras e objetos da análise completa
#'
#' Exporta os resultados da aula para pastas explicitamente informadas.
#' @param matriz Matriz de interação analisada.
#' @param categorias Tabela de identificação e classificação das espécies.
#' @param info Metadados da análise, incluindo tipo; preservados no arquivo RDS.
#' @param figura Objeto ggplot da rede em Sankey.
#' @param grafo Grafo igraph construído a partir de matriz.
#' @param modulos Vetor de módulos na mesma ordem dos vértices do grafo.
#' @param observado Vetor nomeado com as cinco métricas observadas.
#' @param nulas Lista das redes nulas geradas, por modelo.
#' @param simulacoes Lista das matrizes de métricas das redes nulas.
#' @param testes Tabela de inferência produzida por testar_metricas().
#' @param graficos Lista de ggplots retornada por desenhar_testes().
#' @param B Número de simulações usado na análise, registrado nos metadados.
#' @param n_inicios Número de buscas Louvain, registrado nos metadados.
#' @param min_rotulo Limite de peso para nomes, registrado nos metadados.
#' @param pasta_resultados Diretório de tabelas e RDS; padrão "results".
#' @param pasta_figuras Diretório das figuras; padrão "figures".
#' @param n_rotulos Máximo de nomes configurado na análise; padrão 20.
#' @return Caminho de pasta_resultados, invisivelmente. Grava CSVs, analise.rds,
#'   sessionInfo.txt, PNGs e TIFF/LZW a 600 dpi; substitui arquivos homônimos.
#' @details Cria pastas quando necessário. Métricas e testes usam CSV2 (ponto e
#'   vírgula); espécies e interações usam CSV. As sementes gravadas correspondem
#'   aos padrões do roteiro (123, 1102 e 40000). O grafo livre é salvo separadamente
#'   por salvar_grafo_ggplot().
#' @examples
#' \dontrun{
#' salvar_aula(matriz, categorias, info, figura_rede, grafo, modulos, observado,
#'             nulas, simulacoes, testes, graficos_testes, B, n_inicios, min_rotulo,
#'             pasta_resultados = pasta_resultados, pasta_figuras = pasta_figuras)
#' }
salvar_aula = function(matriz, categorias, info, figura, grafo, modulos, observado,
                       nulas, simulacoes, testes, graficos, B, n_inicios, min_rotulo,
                       pasta_resultados = "results", pasta_figuras = "figures", n_rotulos = 20) {
  dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
  dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)
  write.csv(listar_interacoes(matriz), file.path(pasta_resultados, "interacoes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  metricas = data.frame(metrica = unname(nomes_metricas[names(observado)]), valor = unname(observado))
  write.csv2(metricas, file.path(pasta_resultados, "metricas.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv2(testes, file.path(pasta_resultados, "testes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv2(organizar_distribuicoes(simulacoes), file.path(pasta_resultados, "distribuicoes_nulas.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  nos = categorias[match(V(grafo)$name, categorias$id), ]
  nos$grau = as.numeric(degree(grafo))
  nos$soma_interacoes = as.numeric(strength(grafo, weights = E(grafo)$weight))
  nos$modulo = as.integer(modulos)
  write.csv(nos, file.path(pasta_resultados, "especies_modulos.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  ggsave(file.path(pasta_figuras, "rede.png"), plot = figura, width = 14, height = 10, units = "in", dpi = 180, bg = "white")
  ggsave(file.path(pasta_figuras, "rede.tiff"), plot = figura, width = 14, height = 10, units = "in", dpi = 600, compression = "lzw", bg = "white")
  for (modelo in names(graficos)) {
    ggsave(file.path(pasta_figuras, paste0("nulo_", modelo, ".png")), plot = graficos[[modelo]],
           width = 13, height = 5.5, units = "in", dpi = 180, bg = "white")
  }
  saveRDS(list(matriz = matriz, categorias = categorias, info = info, observado = observado,
               nulas = nulas, simulacoes = simulacoes, testes = testes, B = B,
               n_inicios = n_inicios, min_rotulo = min_rotulo, n_rotulos = n_rotulos, semente_observado = 123,
               semente_nulas = 1102, semente_metricas_nulas = 40000), file.path(pasta_resultados, "analise.rds"))
  capture.output(sessionInfo(), file = file.path(pasta_resultados, "sessionInfo.txt"))
  message("Concluído. Figuras em ", pasta_figuras, "; tabelas em ", pasta_resultados, ".")
  invisible(pasta_resultados)
}
# Fontes dos algoritmos: documentação de bipartite::swap.web, stats::r2dtable,
# vegan::nullmodel, vegan::nestednodf e igraph::cluster_louvain.
# 12. Demonstração de alocação sob margens fixas ----
# Amostragem hipergeométrica condicional: mesma distribuição de tabelas que
# r2dtable, com outra implementação para tornar a construção observável.
# Não é uma gravação dos passos internos em C de stats::r2dtable.
#' Mostrar a construção didática de uma tabela com margens fixas
#'
#' Aloca contagens célula a célula por sorteios hipergeométricos condicionais.
#' @param matriz Matriz de contagens inteiras com nomes e margens positivas.
#' @param semente Inteiro para reproduzir a construção; padrão 5101.
#' @return Lista com nome, quadros (matrizes), textos de cada etapa, explicacao e regra.
#' @details A distribuição da tabela final equivale à de r2dtable; os quadros não
#'   registram seus passos internos em C. Estados parciais não são réplicas nulas
#'   completas. A função altera a semente global.
#' @examples
#' \dontrun{
#' trajeto = trajetoria_r2d(matriz)
#' tail(trajeto$quadros, 1)
#' }
trajetoria_r2d = function(matriz, semente = 5101) {
  set.seed(semente)
  m = matriz * 0
  linhas = rowSums(matriz)
  colunas = colSums(matriz)
  quadros = list(m)
  textos = "Início: todas as interações ainda serão alocadas."
  for (i in seq_len(nrow(matriz) - 1)) {
    faltam_linha = linhas[i]
    populacao = sum(colunas)
    for (j in seq_len(ncol(matriz) - 1)) {
      valor = stats::rhyper(1, colunas[j], populacao - colunas[j], faltam_linha)
      m[i, j] = valor
      faltam_linha = faltam_linha - valor
      populacao = populacao - colunas[j]
      colunas[j] = colunas[j] - valor
      if (valor > 0) {
        quadros[[length(quadros) + 1]] = m
        textos = c(textos, paste("Alocar", valor, "em", rownames(m)[i], "×", colnames(m)[j]))
      }
    }
    m[i, ncol(m)] = faltam_linha
    colunas[ncol(m)] = colunas[ncol(m)] - faltam_linha
    if (faltam_linha > 0) {
      quadros[[length(quadros) + 1]] = m
      textos = c(textos, paste("Completar o total de", rownames(m)[i]))
    }
  }
  m[nrow(m), ] = colunas
  quadros[[length(quadros) + 1]] = m
  textos = c(textos, "Matriz completa: totais de todas as linhas e colunas preservados.")
  stopifnot(all(rowSums(m) == rowSums(matriz)), all(colSums(m) == colSums(matriz)))
  list(nome = "r2dtable — construção didática", quadros = quadros, textos = textos,
       explicacao = "Sorteios hipergeométricos condicionais reproduzem a distribuição de r2dtable. Esta versão didática mostra as alocações; não registra o código C do pacote. Quadros intermediários são matrizes incompletas.",
       regra = "Ao completar: totais por espécie fixos; número de ligações pode variar.")
}
# 13. Capturar estados reais da implementação instalada de swap.web ----
# Instrumentar uma cópia local da função, sem alterar o pacote nem seu sorteio.
# Registrar first após inicializações e trocas; o registrador não usa números aleatórios.
#' Registrar uma trajetória real de swap.web
#'
#' Instrumenta uma cópia local da função instalada para registrar estados da matriz.
#' @param matriz Matriz de contagens inteiras com margens positivas.
#' @param semente Inteiro para reproduzir a chamada; padrão 5102.
#' @return Lista com nome, quadros, textos, explicacao e regra da demonstração.
#' @details Não modifica o pacote instalado. As margens permanecem fixas, e L
#'   coincide com o observado ao terminar. Depende do nome interno first; uma
#'   mudança na implementação do pacote pode interromper a captura. Altera a
#'   semente global. Quadros sucessivos não são réplicas independentes de um teste.
#' @examples
#' \dontrun{
#' trajeto = trajetoria_swap(matriz)
#' length(trajeto$quadros)
#' }
trajetoria_swap = function(matriz, semente = 5102) {
  quadros = list()
  #' Registrar estados distintos no ambiente da trajetória
  #' @param m Matriz corrente capturada de swap.web.
  #' @return m invisível; acrescenta o estado à lista quadros do ambiente externo.
  #' @details Auxiliar local de trajetoria_swap(); não fica disponível após source().
  #' @keywords internal
  #' @noRd
  registrar = function(m) {
    if (!length(quadros) || !identical(m, quadros[[length(quadros)]])) quadros[[length(quadros) + 1]] <<- m
    invisible(m)
  }
  #' Inserir chamadas ao registrador após atribuições à matriz first
  #' @param expr Expressão ou chamada R percorrida recursivamente.
  #' @return Expressão R instrumentada, sem executar o sorteio.
  #' @details Auxiliar local de trajetoria_swap(); não modifica o pacote instalado.
  #' @keywords internal
  #' @noRd
  instrumentar = function(expr) {
    if (missing(expr)) return(quote(expr = ))
    if (!is.call(expr)) return(expr)
    expr = as.call(lapply(as.list(expr), instrumentar))
    operador = as.character(expr[[1]])[1]
    if (operador %in% c("<-", "=")) {
      esquerda = expr[[2]]
      alvo = identical(esquerda, as.name("first")) ||
        (is.call(esquerda) && identical(esquerda[[1]], as.name("[")) && identical(esquerda[[2]], as.name("first")))
      if (alvo) return(as.call(list(as.name("{"), expr, quote(registrar(first)))))
    }
    expr
  }
  funcao = bipartite::swap.web
  ambiente = new.env(parent = environment(funcao))
  ambiente$registrar = registrar
  environment(funcao) = ambiente
  body(funcao) = instrumentar(body(funcao))
  set.seed(semente)
  final = funcao(N = 1, web = matriz, verbose = FALSE)[[1]]
  # Uma atualização do pacote que mudar os nomes internos deve ser detectada.
  stopifnot(length(quadros) > 0, isTRUE(all.equal(unname(tail(quadros, 1)[[1]]), unname(final))))
  stopifnot(all(vapply(quadros, function(m) all(rowSums(m) == rowSums(matriz)) && all(colSums(m) == colSums(matriz)), logical(1))),
            sum(final > 0) == sum(matriz > 0))
  textos = c("Sorteio inicial por r2dtable: margens já preservadas.", rep("Troca de pesos em um bloco 2 × 2: margens preservadas.", length(quadros) - 1))
  for (i in seq_along(quadros)[-1]) {
    if (sum(quadros[[i]] != quadros[[i - 1]]) > 4) textos[i] = "Reinício previsto pelo algoritmo: nova matriz com as mesmas margens."
  }
  textos[length(textos)] = paste(textos[length(textos)], "Número de ligações atingiu a meta observada.")
  list(nome = "swap.web — trajetória real", quadros = quadros, textos = textos,
       explicacao = "Estados capturados durante uma chamada real de bipartite::swap.web. O algoritmo parte de r2dtable e altera blocos 2 × 2. O número de ligações coincide com o observado somente ao final.",
       regra = "A cada troca: totais por espécie fixos. Ao terminar: número de ligações também fixo.")
}
# 14. Comparar sorteios binários independentes ----
# O código C do quasiswap não expõe estados internos. Mostrar saídas completas,
# identificadas como sorteios independentes, sem inventar uma trajetória interna.
#' Preparar uma sequência de redes binárias aleatórias completas
#'
#' Sorteia matrizes pelo quasiswap, preservando os graus das duas guildas.
#' @param matriz Matriz binária 0/1 com margens positivas.
#' @param semente Inteiro para reproduzir os sorteios; padrão 5103.
#' @param n_redes Número de matrizes completas a exibir; padrão 20.
#' @return Lista com nome, quadros, textos, explicacao e regra da demonstração.
#' @details Cada quadro é uma saída completa; não mostra os passos internos do
#'   código C. A semente é fornecida a simulate() para reproduzir a sequência.
#' @examples
#' \dontrun{
#' trajeto = trajetoria_quasiswap(1L * (matriz > 0), n_redes = 10)
#' }
trajetoria_quasiswap = function(matriz, semente = 5103, n_redes = 20) {
  cubo = stats::simulate(vegan::nullmodel(matriz, "quasiswap"), nsim = n_redes, seed = semente)
  quadros = lapply(seq_len(n_redes), function(i) cubo[, , i])
  stopifnot(all(vapply(quadros, function(m) all(m %in% c(0, 1)) && all(rowSums(m) == rowSums(matriz)) && all(colSums(m) == colSums(matriz)), logical(1))))
  list(nome = "quasiswap — sorteios completos", quadros = quadros,
       textos = paste("Rede binária completa", seq_len(n_redes), "de", n_redes),
       explicacao = "Cada quadro é uma matriz binária completa gerada por vegan::quasiswap. São sorteios independentes; as transições não representam os passos internos do algoritmo em C.",
       regra = "Em cada rede completa: graus das espécies, número de ligações e valores 0/1 preservados.")
}
# 15. Montar a demonstração interativa em HTML local ----
#' Criar a apresentação HTML dos modelos nulos
#'
#' Carrega um exemplo da disciplina, prepara os quadros e preenche o template local.
#' @param rede Nome do exemplo: "polinizacao" (padrão), "dispersao" ou "herbivoria".
#' @param cores Paleta nomeada por categoria, ou NULL para as cores dos CSVs.
#' @param n_rotulos Máximo de nomes elegíveis na apresentação; padrão 20.
#' @param min_rotulo Soma mínima dos pesos observados para exibir um nome; padrão 10.
#' @return Caminho absoluto do HTML criado em figures/<rede>/; substitui o arquivo
#'   homônimo. O HTML contém os dados e pode ser aberto sem internet.
#' @details Requer data/<rede>_matriz.csv, data/<rede>_categorias.csv, data/redes.csv
#'   e recursos/animacao_nulos.html. A seleção de nomes é fixa entre os quadros;
#'   nomes próximos podem ser omitidos para manter a leitura. O arquivo não abre
#'   automaticamente. As trajetórias usam as sementes de suas funções auxiliares.
#' @examples
#' \dontrun{
#' arquivo = criar_demonstracao("polinizacao", n_rotulos = 12)
#' abrir_demonstracao(arquivo)
#' }
criar_demonstracao = function(rede = "polinizacao", cores = NULL, n_rotulos = 20, min_rotulo = 10) {
  matriz = as.matrix(read.csv(paste0("data/", rede, "_matriz.csv"), row.names = 1, check.names = FALSE, fileEncoding = "UTF-8"))
  categorias = read.csv(paste0("data/", rede, "_categorias.csv"), fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  catalogo = read.csv("data/redes.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
  info = catalogo[catalogo$rede == rede, ]
  conferir_dados(matriz, categorias, info$tipo)
  if (info$tipo == "contagem") {
    trajetorias = list(trajetoria_r2d(matriz), trajetoria_swap(matriz))
  } else {
    trajetorias = list(trajetoria_quasiswap(matriz))
  }
  #' Serializar uma matriz para o desenho da demonstração
  #' @param m Matriz de um quadro, completo ou parcial.
  #' @return Lista com arestas (índices iniciando em zero), L, total e margens.
  #' @details Auxiliar local de criar_demonstracao(), usado na exportação para JSON.
  #' @keywords internal
  #' @noRd
  empacotar = function(m) {
    pos = which(m > 0, arr.ind = TRUE)
    arestas = lapply(seq_len(nrow(pos)), function(k) c(pos[k, 1] - 1, pos[k, 2] - 1, m[pos[k, 1], pos[k, 2]]))
    list(arestas = arestas, L = sum(m > 0), total = sum(m), linhas = unname(rowSums(m)), colunas = unname(colSums(m)))
  }
  nos = categorias[match(c(rownames(matriz), colnames(matriz)), categorias$id), ]
  paleta = unique(nos[, c("grupo", "cor")])
  cores = personalizar_cores(setNames(paleta$cor, paleta$grupo), cores)
  nos$cor = unname(cores[nos$grupo])
  modelos = lapply(trajetorias, function(x) list(nome = x$nome, explicacao = x$explicacao,
                                                 regra = x$regra, textos = x$textos, quadros = lapply(x$quadros, empacotar)))
  dados = list(titulo = info$titulo, tipo = info$tipo, animais = info$animal,
               plantas = rownames(matriz), visitantes = colnames(matriz),
               rotulos = as.list(as.integer(selecionar_rotulos(c(rowSums(matriz), colSums(matriz)), min_rotulo, n_rotulos) - 1)),
               cores = nos$cor, grupos = unique(nos[, c("grupo", "cor")]),
               observado = empacotar(matriz), modelos = modelos)
  template = paste(readLines("recursos/animacao_nulos.html", encoding = "UTF-8", warn = FALSE), collapse = "\n")
  partes = strsplit(template, "/*DADOS*/", fixed = TRUE)[[1]]
  stopifnot(length(partes) == 2)
  json = jsonlite::toJSON(dados, auto_unbox = TRUE, dataframe = "rows", digits = 12)
  pasta = file.path("figures", rede)
  dir.create(pasta, recursive = TRUE, showWarnings = FALSE)
  arquivo = file.path(pasta, "demonstracao_modelos_nulos.html")
  writeLines(enc2utf8(paste0(partes[1], json, partes[2])), arquivo, useBytes = TRUE)
  message("Demonstração salva em ", arquivo)
  normalizePath(arquivo, winslash = "/")
}
#' Abrir a demonstração no navegador
#'
#' Abre o HTML em sessões interativas e informa seu caminho ao usar Rscript.
#' @param arquivo Caminho do HTML retornado por criar_demonstracao().
#' @return O caminho recebido, invisivelmente.
#' @examples
#' \dontrun{
#' abrir_demonstracao("figures/polinizacao/demonstracao_modelos_nulos.html")
#' }
abrir_demonstracao = function(arquivo) {
  if (interactive()) utils::browseURL(arquivo) else message("Abra este HTML no navegador: ", arquivo)
  invisible(arquivo)
}
# 16. Comparar redes observadas lado a lado ----
#' Comparar os grafos e as métricas observadas de várias redes
#'
#' Recalcula as cinco métricas diretamente das matrizes, desenha um grafo por
#' rede e escreve seus valores no painel. Não executa modelos nulos.
#' @param matrizes Lista nomeada de matrizes de interação, com pelo menos duas redes.
#' @param categorias Lista nomeada de tabelas de espécies, com os mesmos nomes de
#'   matrizes. Cada tabela contém id, especie, guilda, grupo e cor.
#' @param tipos Vetor nomeado por rede com "contagem" ou "binaria".
#' @param titulos Vetor de títulos na ordem das matrizes; padrão names(matrizes).
#' @param n_inicios Número de buscas Louvain por rede; padrão 10, igual ao roteiro.
#' @param semente Semente dos módulos e layouts, aplicada a cada rede; padrão 123.
#' @param min_rotulo Soma mínima dos pesos para exibir nomes; padrão 0.
#' @param n_rotulos Máximo de nomes por painel; padrão 8. Inf inclui todos os
#'   elegíveis e 0 oculta os nomes. Não exclui nós ou interações.
#' @param cores Paleta nomeada por grupo, ou NULL para as cores das categorias.
#'   Quando um grupo aparece em mais de uma rede, sua primeira cor é a referência.
#' @param cor_arestas Cor das interações; padrão "#607D8B".
#' @param cor_guias Cor das linhas tracejadas entre nomes e nós; padrão "grey55".
#' @param curvatura Curvatura das interações, entre -1 e 1; padrão 0.12.
#' @param atracao_rotulos Atração positiva dos nomes aos nós; padrão 0.5.
#' @param repulsao_rotulos Repulsão positiva dos nomes; padrão 2.
#' @return Lista com figura (ggplot montado com ggpubr), paineis (ggplots),
#'   metricas (data frame, uma linha por rede), grafos (lista igraph) e modulos
#'   (partições Louvain usadas no cálculo do Q de cada rede).
#' @details A legenda de categorias e a escala da área dos nós são comuns.
#'   As espessuras são relativas aos pesos dentro de cada rede: as unidades
#'   amostrais diferem, portanto não há uma escala quantitativa comum de pesos.
#'   H2' em dados binários é apenas descritivo. O layout usa atração igual nas
#'   arestas e não representa distância ecológica. Define a semente global.
#' @seealso https://rpkgs.datanovia.com/ggpubr/reference/ggarrange.html
#' @examples
#' \dontrun{
#' # Objetos importados no script 03. comparacao_redes.R:
#' comparacao = comparar_redes(matrizes, categorias, tipos, titulos, n_rotulos = 8)
#' print(comparacao$figura)
#' comparacao$metricas
#' }
comparar_redes = function(matrizes, categorias, tipos, titulos = names(matrizes),
                          n_inicios = 10, semente = 123, min_rotulo = 0, n_rotulos = 8,
                          cores = NULL, cor_arestas = "#607D8B", cor_guias = "grey55",
                          curvatura = 0.12, atracao_rotulos = 0.5, repulsao_rotulos = 2) {
  ids = names(matrizes)
  stopifnot(is.list(matrizes), length(matrizes) >= 2, length(matrizes) <= 26,
            length(ids) == length(matrizes), !anyNA(ids), all(nzchar(ids)), !anyDuplicated(ids),
            is.list(categorias), !anyDuplicated(names(categorias)), setequal(names(categorias), ids),
            !anyDuplicated(names(tipos)), setequal(names(tipos), ids),
            length(titulos) == length(ids), !anyNA(titulos), all(nzchar(titulos)))
  categorias = categorias[ids]
  tipos = tipos[ids]
  for (id in ids) conferir_dados(matrizes[[id]], categorias[[id]], tipos[[id]])
  grafos = lapply(matrizes, igraph::graph_from_biadjacency_matrix, directed = FALSE, weighted = TRUE)
  analises = lapply(matrizes, medir_rede, n_inicios = n_inicios, semente = semente, retornar_modulos = TRUE)
  valores = vapply(analises, function(x) x$metricas, setNames(numeric(5), names(nomes_metricas)))
  modulos = lapply(analises, function(x) x$modulos)
  metricas = data.frame(rede = ids, titulo = unname(titulos), tipo = unname(tipos),
                        linhas = vapply(matrizes, nrow, integer(1)),
                        colunas = vapply(matrizes, ncol, integer(1)),
                        nos = vapply(grafos, igraph::vcount, numeric(1)),
                        ligacoes = vapply(grafos, igraph::ecount, numeric(1)),
                        t(valores), row.names = NULL, check.names = FALSE)
  paletas = lapply(categorias, function(x) x[, c("grupo", "cor")])
  paleta = do.call(rbind, paletas)
  paleta = paleta[!duplicated(paleta$grupo), ]
  cores = personalizar_cores(setNames(paleta$cor, paleta$grupo), cores)
  grau_maximo = max(unlist(lapply(grafos, igraph::degree)))
  marcas_grau = pretty(c(0, grau_maximo), n = 3)
  marcas_grau = marcas_grau[marcas_grau > 0 & marcas_grau <= grau_maximo]
  paineis = setNames(vector("list", length(ids)), ids)
  for (i in seq_along(ids)) {
    id = ids[i]
    figura = desenhar_grafo_ggplot(grafos[[id]], categorias[[id]], min_rotulo = min_rotulo,
                                  n_rotulos = n_rotulos, semente = semente, cores = cores,
                                  cor_arestas = cor_arestas, cor_guias = cor_guias,
                                  curvatura = curvatura, atracao_rotulos = atracao_rotulos,
                                  repulsao_rotulos = repulsao_rotulos)
    # Substituir as escalas locais por escalas idênticas antes de reunir as legendas.
    figura$scales$scales = Filter(function(x) !any(x$aesthetics %in% c("fill", "colour", "size", "point.size")),
                                  figura$scales$scales)
    camada_nos = which(vapply(figura$layers, function(x) inherits(x$geom, "GeomPoint"), logical(1)))
    figura$layers[[camada_nos]]$show.legend = TRUE
    valores_texto = formatC(valores[, i], format = "f", digits = 3, decimal.mark = ",")
    texto_metricas = paste(paste(unname(nomes_metricas), valores_texto, sep = " = "), collapse = "\n")
    tipo_texto = if (tipos[[id]] == "binaria") "Presença/ausência" else "Contagens"
    paineis[[id]] = figura +
      scale_fill_manual(name = "Categoria", values = cores, limits = names(cores), drop = FALSE) +
      scale_color_manual(values = cores, limits = names(cores), guide = "none") +
      scale_size_area(name = "Número de\nparceiros", max_size = 7, limits = c(0, grau_maximo),
                      breaks = marcas_grau, aesthetics = c("size", "point.size")) +
      guides(linewidth = "none", point.size = "none") +
      labs(title = paste0("(", letters[i], ") ", titulos[i]),
           subtitle = paste0(metricas$nos[i], " nós | ", metricas$ligacoes[i], " ligações | ", tipo_texto),
           caption = texto_metricas) +
      theme(plot.title = element_text(size = 15, face = "bold", hjust = 0),
            plot.subtitle = element_text(size = 10, hjust = 0),
            plot.caption = element_text(size = 11, hjust = 0, lineheight = 1.15),
            plot.caption.position = "plot", plot.title.position = "plot",
            plot.margin = margin(12, 18, 10, 18))
  }
  figura = ggpubr::ggarrange(plotlist = paineis, ncol = length(paineis), nrow = 1,
                            align = "hv", common.legend = TRUE, legend = "right")
  nota = paste("Área dos nós: grau, com escala comum. Espessura: peso relativo dentro de cada rede.",
               "Valores descritivos; tamanho e amostragem diferem. H2' binário não estima preferência por frequência.", sep = "\n")
  figura = ggpubr::annotate_figure(figura,
                                 top = ggpubr::text_grob("Redes de interações — métricas observadas", family = "sans", face = "bold", size = 17),
                                 bottom = ggpubr::text_grob(nota, family = "sans", size = 10))
  list(figura = figura, paineis = paineis, metricas = metricas, grafos = grafos, modulos = modulos)
}
#' Exportar a comparação das redes e os valores que aparecem nos painéis
#'
#' Salva uma figura conjunta e uma tabela com as métricas sem arredondamento.
#' @param comparacao Lista retornada por comparar_redes() ou comparar_modulos().
#' @param pasta_figuras Diretório de destino das figuras; padrão "figures/comparacao".
#' @param pasta_resultados Diretório da tabela; padrão "results/comparacao".
#' @param nome Nome-base dos arquivos, sem extensão; padrão "comparacao_redes".
#' @param largura Largura da figura em polegadas; padrão 20, adequado a três painéis.
#' @param altura Altura da figura em polegadas; padrão 9.
#' @return Vetor nomeado com os caminhos dos arquivos, invisivelmente. Grava PNG
#'   a 180 dpi, TIFF/LZW a 600 dpi, PDF vetorial e CSV2; substitui arquivos homônimos.
#'   Quando a lista contém especies_modulos, exporta também essa tabela em CSV.
#' @details Cria os diretórios quando necessário. Não repete os cálculos ou sorteios.
#' @examples
#' \dontrun{
#' salvar_comparacao(comparacao, largura = 20, altura = 9)
#' }
salvar_comparacao = function(comparacao, pasta_figuras = "figures/comparacao",
                             pasta_resultados = "results/comparacao", nome = "comparacao_redes",
                             largura = 20, altura = 9) {
  dir.create(pasta_figuras, recursive = TRUE, showWarnings = FALSE)
  dir.create(pasta_resultados, recursive = TRUE, showWarnings = FALSE)
  arquivos = c(png = file.path(pasta_figuras, paste0(nome, ".png")),
                tiff = file.path(pasta_figuras, paste0(nome, ".tiff")),
                pdf = file.path(pasta_figuras, paste0(nome, ".pdf")),
                csv = file.path(pasta_resultados, paste0(nome, ".csv")))
  ggsave(arquivos[["png"]], plot = comparacao$figura, width = largura, height = altura, units = "in", dpi = 180, bg = "white")
  ggsave(arquivos[["tiff"]], plot = comparacao$figura, width = largura, height = altura, units = "in", dpi = 600, compression = "lzw", bg = "white")
  ggsave(arquivos[["pdf"]], plot = comparacao$figura, width = largura, height = altura, units = "in", device = grDevices::cairo_pdf, bg = "white")
  write.csv2(comparacao$metricas, arquivos[["csv"]], row.names = FALSE, fileEncoding = "UTF-8")
  if (!is.null(comparacao$especies_modulos)) {
    arquivos = c(arquivos, especies_modulos = file.path(pasta_resultados, paste0(nome, "_especies_modulos.csv")))
    write.csv(comparacao$especies_modulos, arquivos[["especies_modulos"]], row.names = FALSE, fileEncoding = "UTF-8")
  }
  invisible(arquivos)
}
# 17. Evidenciar os módulos de Louvain ----
#' Construir contornos convexos ao redor dos módulos
#'
#' Expande as posições dos nós por uma pequena margem e calcula um polígono
#' convexo para cada módulo. Também representa módulos com um ou dois nós.
#' @param nos Data frame com x, y (coordenadas) e modulo (identificador do módulo).
#' @param margem Número positivo na escala do layout, padrão 0.045; aumenta a
#'   folga ao redor dos nós. Não representa incerteza ou distância ecológica.
#' @param n_pontos Inteiro de pontos por círculo de expansão; padrão 32, mínimo 8.
#' @return Data frame com x, y e modulo, contendo os vértices ordenados dos polígonos.
#' @details Usa grDevices::chull(). Os polígonos podem se sobrepor e conter nós de
#'   outros módulos em seu interior; o pertencimento é definido pelo Louvain e
#'   pela cor do nó. Não são regiões geográficas nem intervalos de confiança.
#' @seealso https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/chull.html
#' @examples
#' pontos = data.frame(x = c(0, 1, 0.5), y = c(0, 0, 1), modulo = c(1, 1, 2))
#' contornar_modulos(pontos)
contornar_modulos = function(nos, margem = 0.045, n_pontos = 32) {
  stopifnot(all(c("x", "y", "modulo") %in% names(nos)), nrow(nos) > 0,
            all(is.finite(nos$x)), all(is.finite(nos$y)), !anyNA(nos$modulo),
            length(margem) == 1, is.finite(margem), margem > 0,
            length(n_pontos) == 1, is.finite(n_pontos), n_pontos >= 8, n_pontos == floor(n_pontos))
  angulos = seq(0, 2 * pi, length.out = n_pontos + 1)[-1]
  poligonos = lapply(split(nos, nos$modulo, drop = TRUE), function(pontos) {
    x = as.vector(outer(pontos$x, margem * cos(angulos), "+"))
    y = as.vector(outer(pontos$y, margem * sin(angulos), "+"))
    borda = grDevices::chull(x, y)
    data.frame(x = x[borda], y = y[borda], modulo = pontos$modulo[1])
  })
  resultado = do.call(rbind, poligonos)
  rownames(resultado) = NULL
  resultado
}
#' Comparar as redes coloridas pelos módulos de Louvain
#'
#' Reutiliza os grafos, coordenadas, nomes e métricas de uma comparação existente,
#' substituindo as cores das categorias por cores dos módulos identificados.
#' @param comparacao Lista criada por comparar_redes(), incluindo modulos.
#' @param cores_modulos Vetor nomeado por M1, M2 etc. para alterar cores, ou NULL
#'   para uma paleta qualitativa HCL (Dark 3). As cores são locais a cada rede.
#' @param mostrar_poligonos Lógico; TRUE (padrão) desenha áreas atrás das arestas e nós.
#' @param alpha_poligonos Número de 0 a 1, padrão 0.10; transparência das áreas.
#' @param margem_poligonos Número positivo, padrão 0.045; expansão dos contornos
#'   na escala do layout. Usado somente quando mostrar_poligonos = TRUE.
#' @return Lista com figura, paineis, metricas, grafos, modulos, poligonos e
#'   especies_modulos (rede, espécie, categoria original, módulo, cor, grau e força).
#'   Pode ser exportada por salvar_comparacao().
#' @details Não refaz o Louvain: usa a mesma partição que gerou o Q observado,
#'   verifica sua concordância com a tabela e mantém todas as interações. A área
#'   dos nós e os pesos das arestas conservam as escalas da primeira figura.
#'   Os textos mantêm a seleção original e recebem a cor de seu módulo.
#'   M1 em uma rede não equivale a M1 em outra. Contornos podem se sobrepor;
#'   cores dos nós, e não a área sombreada, indicam o pertencimento.
#' @examples
#' \dontrun{
#' comparacao_modulos = comparar_modulos(comparacao, alpha_poligonos = 0.10)
#' print(comparacao_modulos$figura)
#' salvar_comparacao(comparacao_modulos, nome = "comparacao_modulos")
#' }
comparar_modulos = function(comparacao, cores_modulos = NULL, mostrar_poligonos = TRUE,
                            alpha_poligonos = 0.10, margem_poligonos = 0.045) {
  if (is.null(comparacao$modulos)) stop("Recrie comparacao com comparar_redes() após carregar o setup atualizado.", call. = FALSE)
  stopifnot(is.logical(mostrar_poligonos), length(mostrar_poligonos) == 1, !is.na(mostrar_poligonos),
            length(alpha_poligonos) == 1, is.finite(alpha_poligonos), alpha_poligonos >= 0, alpha_poligonos <= 1)
  ids = names(comparacao$grafos)
  stopifnot(setequal(ids, names(comparacao$paineis)), setequal(ids, names(comparacao$modulos)))
  maior_modulo = max(unlist(comparacao$modulos))
  nomes_modulos = paste0("M", seq_len(maior_modulo))
  cores = personalizar_cores(setNames(grDevices::hcl.colors(maior_modulo, palette = "Dark 3"), nomes_modulos), cores_modulos)
  paineis = poligonos = especies = setNames(vector("list", length(ids)), ids)
  metricas = comparacao$metricas
  metricas$n_modulos = NA_integer_
  for (id in ids) {
    grafo = comparacao$grafos[[id]]
    modulos = comparacao$modulos[[id]]
    linha = match(id, metricas$rede)
    stopifnot(length(modulos) == igraph::vcount(grafo), !anyNA(modulos),
              all(modulos > 0), all(modulos == as.integer(modulos)),
              abs(igraph::modularity(grafo, modulos, weights = igraph::E(grafo)$weight) - metricas$Q[linha]) < 1e-10)
    figura = comparacao$paineis[[id]]
    nos = figura$data
    nos$modulo = as.integer(modulos[match(nos$id, igraph::V(grafo)$name)])
    especies[[id]] = data.frame(rede = id, id = nos$id, especie = nos$especie,
                                 guilda = nos$guilda, categoria = as.character(nos$grupo),
                                 modulo = nos$modulo, cor_modulo = unname(cores[paste0("M", nos$modulo)]),
                                 grau = nos$grau, soma_interacoes = nos$forca)
    nos$grupo = factor(paste0("M", nos$modulo), levels = nomes_modulos)
    figura$data = nos
    # Clonar as escalas impede que esta versão altere a figura por categorias.
    figura$scales = figura$scales$clone()
    figura$scales$scales = Filter(function(x) !any(x$aesthetics %in% c("fill", "colour")), figura$scales$scales)
    if (mostrar_poligonos) {
      contornos = contornar_modulos(nos, margem = margem_poligonos)
      contornos$grupo = factor(paste0("M", contornos$modulo), levels = nomes_modulos)
      poligonos[[id]] = contornos
      figura = figura +
        geom_polygon(data = contornos, aes(x = x, y = y, group = grupo, fill = grupo, color = grupo),
                     inherit.aes = FALSE, alpha = alpha_poligonos, linewidth = 0.35, show.legend = FALSE)
      figura$layers = figura$layers[c(length(figura$layers), seq_len(length(figura$layers) - 1))]
    }
    n_modulos = length(unique(modulos))
    metricas$n_modulos[linha] = n_modulos
    paineis[[id]] = figura +
      scale_fill_manual(name = "Módulo\n(em cada rede)", values = cores, limits = nomes_modulos, drop = FALSE) +
      scale_color_manual(values = cores, limits = nomes_modulos, guide = "none") +
      labs(subtitle = paste0(figura$labels$subtitle, "\n", n_modulos, " módulos — Louvain ponderado"))
  }
  figura = ggpubr::ggarrange(plotlist = paineis, ncol = length(paineis), nrow = 1,
                            align = "hv", common.legend = TRUE, legend = "right")
  nota = paste("Cores e números identificam módulos dentro de cada rede; não indicam equivalência entre redes.",
               "Contornos são guias visuais e podem se sobrepor. Área dos nós: grau comum; espessura: peso relativo em cada rede.",
               "Métricas observadas mantidas; H2' binário não estima preferência por frequência.", sep = "\n")
  figura = ggpubr::annotate_figure(figura,
                                 top = ggpubr::text_grob("Redes de interações — módulos de Louvain", family = "sans", face = "bold", size = 17),
                                 bottom = ggpubr::text_grob(nota, family = "sans", size = 9))
  especies = do.call(rbind, especies)
  rownames(especies) = NULL
  list(figura = figura, paineis = paineis, metricas = metricas, grafos = comparacao$grafos,
       modulos = comparacao$modulos, poligonos = poligonos, especies_modulos = especies)
}
