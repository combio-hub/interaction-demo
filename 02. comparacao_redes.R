# COMPARAR AS TRÊS REDES E SUAS MÉTRICAS OBSERVADAS ----
# Execute a partir da pasta do projeto. Este script lê os CSVs originais da aula.
# Os valores são recalculados; não é necessário executar os modelos nulos antes.
# 1. Carregar as funções e escolher os arquivos ----
source("00. setup.R", encoding = "UTF-8")
arquivos_matrizes = c(polinizacao = "data/polinizacao_matriz.csv",
                      dispersao = "data/dispersao_matriz.csv",
                      herbivoria = "data/herbivoria_matriz.csv")
arquivos_categorias = c(polinizacao = "data/polinizacao_categorias.csv",
                        dispersao = "data/dispersao_categorias.csv",
                        herbivoria = "data/herbivoria_categorias.csv")
tipos = c(polinizacao = "contagem", dispersao = "contagem", herbivoria = "binaria")
titulos = c("Polinização", "Dispersão", "Herbivoria")
# Para dados da turma, edite os caminhos, tipos e títulos; mantenha os mesmos nomes nas listas.
# 2. Importar as matrizes e categorias ----
matrizes = lapply(arquivos_matrizes, read.csv, row.names = 1, check.names = FALSE, fileEncoding = "UTF-8")
matrizes = lapply(matrizes, as.matrix)
categorias = lapply(arquivos_categorias, read.csv, fileEncoding = "UTF-8", stringsAsFactors = FALSE)
# 3. Configurar os painéis ----
n_inicios = 10
semente = 123
min_rotulo = 0
n_rotulos = 8
# n_rotulos vale por rede. Zero oculta os nomes; Inf inclui todos os elegíveis.
# O número de nós e interações permanece completo, independentemente dos nomes.
cores_grupos = NULL
# Exemplo: cores_grupos = c(Plantas = "darkgreen", Aves = "#0072B2", Morcegos = "#CC79A7")
cor_ligacoes = "#607D8B"
cor_guias = "grey55"
atracao_rotulos = 2
repulsao_rotulos = 50
# 4. Calcular e comparar ----
# São as mesmas definições da aula: C, NODF binário, Q ponderado de Louvain, L/S e H2'.
# Os valores aparecem abaixo de cada grafo; as métricas constantes também são descritivas.
# Nenhum teste é realizado aqui: as figuras não mostram p-valores nem SES.
comparacao = comparar_redes(matrizes, categorias, tipos, titulos,
                            n_inicios = n_inicios, semente = semente,
                            min_rotulo = min_rotulo, n_rotulos = n_rotulos, cores = cores_grupos,
                            cor_arestas = cor_ligacoes, cor_guias = cor_guias, curvatura = 0.12,
                            atracao_rotulos = atracao_rotulos, repulsao_rotulos = repulsao_rotulos)
comparacao$metricas
print(comparacao$figura)
# Mesma escala de grau nos três painéis; cores e legenda comuns para as categorias.
# As espessuras variam dentro de cada rede, pois as unidades amostrais são diferentes.
# H2' na herbivoria binária é descritivo; não mede preferência por frequência.
# 5. Salvar figura e tabela ----
arquivos = salvar_comparacao(comparacao, largura = 20, altura = 9)
arquivos
# Saídas: figures/comparacao/comparacao_redes.png, .tiff e .pdf.
# Tabela sem arredondamento: results/comparacao/comparacao_redes.csv (read.csv2).
# PERGUNTAS: redes com conectância semelhante têm os mesmos módulos ou parceiros?
# Por que não podemos interpretar maior Q ou H2' como maior complexidade ecológica?
# 6. Comparar as mesmas redes pelos módulos de Louvain ----
# Usa exatamente as partições que produziram Q na seção 4; não executa outro sorteio.
# Mantém posições, nomes selecionados, tamanhos dos nós, ligações e métricas.
cores_modulos = NULL
# Exemplo: cores_modulos = c(M1 = "#0072B2", M2 = "#D55E00", M3 = "#009E73")
mostrar_poligonos = TRUE
alpha_poligonos = 0.10
margem_poligonos = 0.045
# FALSE remove as áreas. Aumentar alpha_poligonos torna as áreas mais visíveis.
# A margem controla a folga dos contornos ao redor dos nós, na escala do layout.
comparacao_modulos = comparar_modulos(comparacao, cores_modulos = cores_modulos,
                                      mostrar_poligonos = mostrar_poligonos,
                                      alpha_poligonos = alpha_poligonos, margem_poligonos = margem_poligonos)
print(comparacao_modulos$figura)
comparacao_modulos$metricas
# Cores agora representam módulos, que podem conter espécies das duas guildas.
# M1 em uma rede não equivale a M1 em outra; os números são identificadores locais.
# Polígonos podem se sobrepor. A cor do nó identifica seu módulo, não sua posição
# dentro de uma área sombreada. Essas áreas não são limites ecológicos ou geográficos.
# 7. Salvar a versão por módulos e a identificação das espécies ----
arquivos_modulos = salvar_comparacao(comparacao_modulos, nome = "comparacao_modulos", largura = 20, altura = 9)
arquivos_modulos
# Figura: figures/comparacao/comparacao_modulos.png, .tiff e .pdf.
# Tabelas: results/comparacao/comparacao_modulos.csv (read.csv2) e
# comparacao_modulos_especies_modulos.csv (read.csv), com módulo de cada espécie.
# PERGUNTA: os módulos coincidem com as categorias biológicas da primeira figura?
