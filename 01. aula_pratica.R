# AULA PRÁTICA: REDES DE INTERAÇÕES ECOLÓGICAS ----
# Objetivo: representar uma rede, calcular métricas e compará-la com modelos nulos.
# Execute por seções no RStudio. Faça as perguntas antes de executar os cálculos.
# Fontes conceituais:
# https://mlurgi.github.io/networks_for_r/lesson-2.html
# https://mlurgi.github.io/networks_for_r/lesson-4.html
# https://mlurgi.github.io/networks_for_r/lesson-5.html
# https://fukamilab.github.io/BIO202/09-B-networks.html
# 1. Escolher a rede e carregar o setup ----
source("00. setup.R", encoding = "UTF-8")

# Carregar planilhas
matriz = read.csv("data/polinizacao_matriz.csv", row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
matriz
matriz = as.matrix(matriz)
categorias = read.csv("data/polinizacao_categorias.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
categorias
catalogo = read.csv("data/redes.csv", fileEncoding = "UTF-8", stringsAsFactors = FALSE)
catalogo

redes = unique(catalogo$rede)

pasta_resultados = file.path("results", redes[1])
pasta_figuras = file.path("figures", redes[1])

# 2. Conferir as entradas ----
conferir_dados(matriz, categorias, info$tipo)
head(matriz)
categorias
info
table(categorias$grupo)

# 3. Parâmetros da aula ----
B = 199 # Número de permutações dos testes
# B = 199 agiliza a demonstração; use B = 999 para aumentar a resolução do p.
n_inicios = 10
min_rotulo = 10
n_rotulos = 20
# min_rotulo é a soma mínima de pesos; n_rotulos é o máximo de nomes exibidos.
# Use min_rotulo = 0 para selecionar somente pela quantidade, priorizando maiores somas.
# n_rotulos = Inf inclui todos os elegíveis; n_rotulos = 0 oculta os nomes.
# Esses controles afetam os textos dos gráficos.

## 3.1. Descrever a rede e construir o grafo ----
binaria = 1L * (matriz > 0)
P = nrow(matriz)
A = ncol(matriz)
S = P + A
L = sum(binaria)
total_registros = sum(matriz)
data.frame(plantas = P, animais = A, nos = S, ligacoes = L, soma_da_matriz = total_registros)
# Grafo: nós representam táxons/unidades; arestas representam os pares interagentes.
# Usamos arestas não direcionadas para representar a associação entre as guildas.
grafo = graph_from_biadjacency_matrix(matriz, directed = FALSE, weighted = TRUE)
grafo
# Grau: número de parceiros. Força: soma dos pesos das interações do nó.
grau = degree(grafo)
forca = strength(grafo, weights = E(grafo)$weight)
head(sort(grau, decreasing = TRUE))
head(sort(forca, decreasing = TRUE))
# 4. Visualizar a rede com os nós distribuídos pelo plano ----
# Fruchterman–Reingold: atração entre nós ligados e repulsão entre os nós.
# As duas guildas ficam distribuídas pelo plano, sem separação em duas colunas.
# O desenho continua representando as mesmas interações planta–animal da matriz.
# Cor = categoria; área dos nós = grau; espessura das arestas = peso da interação.
# Na rede binária, todas as arestas têm o mesmo peso e a mesma espessura.
# O limite min_rotulo afeta somente os nomes, com base na soma dos pesos do nó.
# Curvatura = 0 produz linhas retas; 0.12 produz curvas suaves.
# Por padrão, o layout considera os parceiros e usa atração igual nas ligações.
# Use ponderar_layout = TRUE para dar mais atração às interações de maior peso.
# As distâncias são gráficas: não representam distância geográfica ou uma métrica.
# https://r.igraph.org/reference/layout_with_fr.html
# Personalize as cores abaixo usando nomes do R ("darkgreen") ou códigos hexadecimais.
# A paleta vem das categorias importadas, inclusive quando forem dados dos alunos.
# cores = NULL usa as cores originais. Também é possível mudar apenas uma categoria.
paleta = unique(categorias[, c("grupo", "cor")])
cores_grupos = setNames(paleta$cor, paleta$grupo)
# Exemplo de edição: cores_grupos["Plantas"] = "#009E73" (use um grupo dos seus dados).
cor_ligacoes = "#607D8B"
cor_texto = "black"
# Os nomes recebem as cores dos nós e continuam em itálico.
# ggrepel aproxima os nomes e os desloca localmente para reduzir sobreposições.
# As linhas-guia tracejadas identificam os nós; não representam novas interações.
# Aumentar atracao_rotulos aproxima os nomes; aumentar repulsao_rotulos abre espaço.
# posicao_rotulos = "externa" usa as margens, quando isso for desejado.
figura_grafo = desenhar_grafo_ggplot(grafo, categorias, min_rotulo, disposicao = "livre", curvatura = 0.12,
                                   cores = cores_grupos, cor_arestas = cor_ligacoes, cor_borda = "grey20",
                                   cor_texto = cor_texto, cor_guias = "grey55", posicao_rotulos = "proxima",
                                   n_rotulos = n_rotulos, atracao_rotulos = 2, repulsao_rotulos = 50)
print(figura_grafo)
# PERGUNTA: o nó de maior grau também participa das interações de maior peso?

## 4.1. Comparar com a representação em Sankey, se desejado ----
# O objeto abaixo é preparado para a exportação; descomente print para visualizá-lo.
# No Sankey, as larguras dos nós representam a soma dos pesos, não o grau.
figura_rede = desenhar_rede(matriz, categorias, info, min_rotulo, cores = cores_grupos, cor_texto = cor_texto, n_rotulos = n_rotulos,cor_fluxos = "grey70")
print(figura_rede)

# 5. Conectância: fração dos pares possíveis que interagem ----
# Para duas guildas com P plantas e A animais, existem P × A pares possíveis.
# C = L/(P × A), entre 0 e 1. Conte pares, não a soma das frequências.
# Valores maiores indicam maior preenchimento da matriz.
conectancia = L / (P * A)
conectancia

# 6. Aninhamento: parceiros de especialistas como subconjuntos ----
# Uma rede aninhada tende a apresentar espécies com poucos parceiros interagindo
# com um subconjunto dos parceiros das espécies mais generalistas.
# NODF considera sobreposição de parceiros e diferenças nos graus de linhas/colunas.
# Usamos NODF binário: 0 a 100, com valores maiores indicando maior aninhamento.
# Um NODF alto, sozinho, não demonstra diferença em relação ao acaso.
nodf = unname(vegan::nestednodf(binaria, order = TRUE, weighted = FALSE)$statistic["NODF"])
nodf

# 7. Modularidade: compartimentos de interação ----
# Módulos reúnem nós com mais ligações internas do que seria esperado pelo modelo
# de referência do índice Q, com relativamente menos ligações entre módulos.
# Louvain procura uma partição; Q avalia essa partição. O resultado é aproximado.
# Repetimos a busca 10 vezes e mantemos o maior Q, também em cada simulação nula.
# Aqui Q é a modularidade geral ponderada do igraph; não é Q específico bipartido.
modulos = detectar_modulos(grafo, n_inicios, semente = 123)
Q = modularity(grafo, modulos, weights = E(grafo)$weight)
Q
table(modulos)

# 8. Especialização: duas perguntas diferentes ----
## 8.1. Quantos parceiros existem na rede? ----
# L/S = número de ligações dividido pelo número total de nós.
# Valores menores indicam menos ligações por nó nesse índice de rede.
# O grau médio é 2L/S: cada aresta participa do grau de dois nós.
ligacoes_por_especie = L / S
grau_medio = mean(grau)
c(ligacoes_por_especie = ligacoes_por_especie, grau_medio = grau_medio)

## 8.2. Quão específica é a distribuição das interações? ----
# H2' compara as frequências observadas com a expectativa baseada nos totais
# marginais. Varia de 0 (menor especialização) a 1 (maior especialização).
# Não equivale a ter poucos parceiros; também considera como os pesos se distribuem.
# Totais amostrais não são medidas independentes de abundância ou preferência.
# Na herbivoria, H2' é descritivo com pesos unitários; não mede preferência por frequência.
H2 = as.numeric(bipartite::networklevel(matriz, index = "H2", H2_integer = TRUE))
H2
observado = c(conectancia = conectancia, NODF = nodf, Q = Q, ligacoes_por_especie = ligacoes_por_especie, H2 = H2)
data.frame(metrica = unname(nomes_metricas), valor = unname(observado))

# 9. Testar: a estrutura difere do esperado ao acaso? ----
# Redes nulas redistribuem interações, preservando propriedades escolhidas.
# Nas contagens: swap.web mantém os totais por espécie e L; testa NODF, Q e H2'.
# r2dtable mantém os totais por espécie e deixa L variar; testa C e L/S.
# Na matriz binária: quasiswap mantém graus e valores 0/1; testa NODF e Q.
# Nessa condição, C, L/S e H2' não variam e ficam fora da tabela de testes.
# r2dtable não é aplicado à herbivoria: geraria valores de frequência inexistentes.
# Todas as combinações entre guildas são permitidas; não impomos interações proibidas.
# swap.web pode favorecer pesos extremos. A inferência é condicional ao algoritmo.
nulas = gerar_nulas(matriz, info$tipo, B)
names(nulas)
simulacoes = medir_nulas(nulas, n_inicios)
testes = testar_metricas(observado, simulacoes)
testes[, c("modelo", "metrica", "observado", "media_nula", "SES", "p_bilateral", "p_holm", "status")]

# 10. Interpretar tamanho do efeito e valores de P ----
# SES = (observado - média nula)/desvio-padrão nulo.
# SES positivo: acima da referência; negativo: abaixo. Não é um p-valor.
# p_bilateral = 2 × a menor proporção de cauda, limitado a 1.
# Contamos empates e usamos (1 + número de extremos)/(B + 1) em cada cauda.
# Com 999 simulações, o menor p bilateral possível é 0,002; não significa p = 0.
# p_holm corrige o conjunto de testes: compare-o com o nível escolhido, como 0,05.
# C e L/S são proporcionais nesta matriz: representam um único teste no ajuste.
# São quatro testes para contagens e dois para a matriz binária, se houver variação.
# Os quantis 2,5% e 97,5% descrevem a distribuição nula; não são ICs do observado.
# Não escolher uma cauda após ver o resultado nem usar |SES| > 2 como substituto do teste.
cores_nulos = c(swap.web = "#009E73", r2dtable = "#56B4E9", quasiswap = "#0072B2")
graficos_testes = desenhar_testes(simulacoes, testes, cores = cores_nulos,
                                 cor_observado = "#D55E00", cor_borda = "white", cor_texto = cor_texto)
print(graficos_testes)