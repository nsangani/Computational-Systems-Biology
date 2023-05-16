# install network analysis packages
install.packages('igraph')
install.packages("visNetwork")

#Load Library
library(igraph)

# Load data and preview
data = read.csv(file.choose(), header = T, sep = '\t')
View(head(data))
colnames(data)

# Make dataframe of the vertices/nodes
y = data.frame(data$OFFICIAL_SYMBOL_A, data$OFFICIAL_SYMBOL_B)

# Create  network to answer the questions:
net = graph.data.frame(y,directed = F) # Undirected network

# (1) ==== For each node compute the degree and identify top 10 hub nodes ======
V(net) # 26255 nodes
E(net) # 979553 edges

# Network measures
labels = V(net)$label
degrees = V(net)$degree # degree for each node
ecount(net)

hb = data.frame(labels, degrees) #hubs dataframe to search for top 10 hubs
hb[order(hb$degrees, decreasing = T),][1:10,] 

# (2) ====== Draw the graph of degree distribution to see if it is 
#     a scale free network or not? =============================================
# Histogram
hist(V(net)$degree,
     main = 'Degree Distribution',
     col = 'darkmagenta',
     xlab = 'Degree of verticies')

degree_distribution(net, cumulative = FALSE)

# (3) ===== Compute the shortest path for all node pairs and 
# plot the distribution of shortest path length. ===============================

sd = distance_table(
  net, directed = F)


View(sd)

# (4) For each node, compute its betweenness and closeness, and 
# identify top 10 nodes with highest betweenness and closeness respectively. 
btw = betweenness(
  net,
  v = V(net),
  directed = F,
  weights = NULL,
  nobigint = TRUE,
  normalized = FALSE
)

cls = closeness(
  net,
  vids = V(net),
  mode = c("all"),
  weights = NULL,
  normalized = FALSE
)

plot(btw)

View(btw) # Betweenness of all nodes
View(cls) # Closeness of all nodes

tail(sort(btw),10) #Top 10 betweenness
tail(sort(cls),10) #Top 10 closeness






















library("igraph")
library("poweRlaw")
library("ggplot2")

# Just loading my data

G <- labels

# List of degrees
G.degrees = degrees

# Let's count the frequencies of each degree
G.degree.histogram <- as.data.frame(table(G.degrees))

# Need to convert the first column to numbers, otherwise
# the log-log thing will not work (that's fair...)
G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])

# Now, plot it!
ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
  geom_point() +
  scale_x_continuous("Degree\n(nodes with this amount of connections)",
                     breaks = c(1, 3, 10, 30, 100, 300),
                     trans = "log10") +
  scale_y_continuous("Frequency\n(how many of them)",
                     breaks = c(1, 3, 10, 30, 100, 300, 1000),
                     trans = "log10") +
  ggtitle("Degree Distribution (log-log)") +
  theme_bw()

