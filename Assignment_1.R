library('rstudioapi') 
library('igraph')
library('ggplot2')

setwd(dirname(getActiveDocumentContext()$path))
data <- read.csv(file.choose(), header = T, sep = '\t')
ppi_Frame <- data.frame(data$OFFICIAL_SYMBOL_A, data$OFFICIAL_SYMBOL_B)
ppi_Graph <- graph_from_data_frame(ppi_Frame, directed = FALSE)

#1
deg_List <- data.frame('Vertex' = vertex_attr(ppi_Graph,'name'),'Degree' = degree(ppi_Graph, V(ppi_Graph)))
deg_List <- deg_List[order(-deg_List$Degree),]
top_Nodes <- head(deg_List, 10)
View(top_Nodes)

#2
deg_Dist <- degree_distribution(ppi_Graph,v = V(ppi_Graph))
deg_Frame <- data.frame('Degree' = log10(1:length(deg_Dist)), 'Probability' = log10(deg_Dist))
deg_Frame <- deg_Frame[is.finite(deg_Frame$Degree) & is.finite(deg_Frame$Probability),]
ggplot(deg_Frame, aes(x = Degree, y = Probability)) + geom_point() + stat_smooth(method='lm') + labs(title='Degree Distribution', x='log(Degree)', y='log(Probability')
#ggsave('Degree_Distribution.png', dpi=300)
print(lm(deg_Frame$Probability ~ deg_Frame$Degree))

#3
path_Tab <- data.frame('PathLength' = path.length.hist(ppi_Graph)$res)
path_Tab['Index'] <- c(1:dim(path_Tab)[1])
ggplot(path_Tab, aes(x = Index, y = PathLength)) + geom_col() + labs(title = 'Path Length Distribution', x = 'Path Length', y = 'Frequency') + scale_y_continuous(trans = 'log10')
#ggsave('Path_Legnth_Distribution.png', dpi=300)

#4
close_List <- data.frame('Vertex' = vertex_attr(ppi_Graph,'name'), 'Closeness' = closeness(ppi_Graph, v = V(ppi_Graph)))
close_List <- close_List[order(-close_List$Closeness),]
top_Close <- head(close_List, 10)
View(top_Close)

betw_List <- data.frame('Vertex' = vertex_attr(ppi_Graph,'name'), 'Betweenness' = betweenness(ppi_Graph, v = V(ppi_Graph)))
betw_List <- betw_List[order(-betw_List$Betweenness),]
top_Betw <- head(betw_List, 10)
View(top_Betw)
