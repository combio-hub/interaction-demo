# DEMONSTRAÇÃO VISUAL DOS MODELOS NULOS ----
# Material complementar à aula principal. Todas as funções estão no setup.
source("00. setup.R", encoding = "UTF-8")
rede = "polinizacao"
# Também aceita "dispersao" e "herbivoria".
# Nas contagens, compare duas construções:
# 1. Alocação condicional sob margens fixas, com a distribuição de r2dtable.
#    Esta versão didática mostra os sorteios célula a célula.
# 2. Trajetória real de swap.web: parte de uma matriz r2dtable e redistribui
#    pesos em blocos 2 × 2 até igualar o número de ligações ao observado.
# Na herbivoria, a apresentação mostra redes completas de quasiswap.
# As transições entre elas não são apresentadas como passos internos do algoritmo.
arquivo_demonstracao = criar_demonstracao(rede, n_rotulos = 20, min_rotulo = 10)
# Para personalizar as cores: criar_demonstracao(rede, cores = c(Plantas = "#009E73")).
abrir_demonstracao(arquivo_demonstracao)
# No navegador: selecione o modelo, clique em Reproduzir ou avance uma etapa.
# Use a barra para examinar a montagem e alterne entre grafo e matriz.
# PERGUNTAS:
# Quais propriedades permanecem constantes em cada procedimento?
# Por que o número de ligações muda durante swap.web, mas coincide ao final?
# Por que os quadros parciais de uma construção não são réplicas de um teste?
