library(igraph)
library(ggraph)
library(boot)
library(SIN)

# Quick analysis ----------------------------------------------------------

data <- mts # Load Data
hist(data) # Visualize Data

# Variables
n.features <- ncol(data)
n <- nrow(data) # number of observations
data.cor <- cor(data) # matrix of correlations of the sample
feature.names <- colnames(data)

# Initialize edge matrix filled with NA
edge <- matrix(NA, n.features, n.features)
colnames(edge) <- feature.names
rownames(edge) <- feature.names

# Non Parametric Bootstrap ------------------------------------------------

set.seed(123)
B <- 1000
delta.beta <- function(data, indices) {
  d <- data[indices,] # allows boot to select sample
  pearson <- cor(d) # creates a correlation matrix for each boostrap
  data.cor <- cor(data) # matrix of correlations of the sample
  delta <- sqrt(nrow(data))*max(abs(pearson-data.cor))
  return(delta)
}

bootdelta <- boot(data = data, statistic = delta.beta, R = B)

# Now we have the list of delta.b and we can compute the ecdf
delta.b <- bootdelta$t
f.hat <- ecdf(delta.b)

# Plot the ecdf
plot(f.hat, ylab="Fn(x)", verticals = FALSE, lwd = 3)

# let's start with alpha = 0.95
alpha1 <- 0.05
t.alpha = quantile(f.hat, probs = (1-alpha1), type =1) # inverse of the ECDF function
ci.lower = data.cor-t.alpha/sqrt(n) # lower bound of the CI
ci.upper = data.cor+t.alpha/sqrt(n) # upper bound of the CI

alpha2 <- 0.01
t.alpha2 = quantile(f.hat, probs = (1-alpha2), type =1) # inverse of the ECDF function
ci.lower2 = data.cor-t.alpha2/sqrt(n) # lower bound of the CI
ci.upper2 = data.cor+t.alpha2/sqrt(n) # upper bound of the CI


# Edge matrix -------------------------------------------------------------

# 1. EPSILON = 0

eps <- 0 # we initially set epsilon as 0

for (i in 1:n.features){
  for(j in 1:n.features){
    
    # If 0 is not included in the CI it means there is an edge between the two nodes 
    if ((ci.lower[i,j] < eps & eps < ci.upper[i, j]) == FALSE) edge[i,j] = 1 # we set one in the adjacency matrix
    else edge[i,j] = 0 # if 0 is included in CI we set 0 in the adjacency matrix
    if (i == j) edge[i,j] = 0 # self edge is excluded
  }
}

graph.eps1 = graph.adjacency(edge, mode = "undirected")
# tkplot(graph.eps1)

# Density 
density.eps1 <- edge_density(graph.eps1, loops=FALSE)
# the density of this graph is very high, it's not a complete graph but 

# Connected components
components(graph.eps1) # all nodes are connected

# Number of edges 
sum(edge)

# 2. EPSILON = 0.5

eps2 <- 0.5 # we initially set epsilon as 0
edge2 <- matrix(NA, n.features, n.features)

for (i in 1:n.features){
  for(j in 1:n.features){
    # If 0 is not included in the CI it means there is an edge between the two nodes
    if (eps2 < ci.lower[i, j] | -eps2 > ci.upper[i, j]) edge2[i,j] = 1
    else edge2[i,j] = 0
    if (i == j) edge2[i,j] = 0 # self edge is excluded
  }
}
graph.eps2 = graph.adjacency(edge2, mode = "undirected")

# Components
components(graph.eps2)

# Number of edges
sum(edge2) # number of edges


# 3. EPSILON = 0.65

eps3 <- 0.65 # we initially set epsilon as 0
edge3 <- matrix(NA, n.features, n.features)

for (i in 1:n.features){
  for(j in 1:n.features){
    # If 0 is not included in the CI it means there is an edge between the two nodes
    if (eps3 < ci.lower[i, j] | -eps3 > ci.upper[i, j]) edge3[i,j] = 1
    else edge3[i,j] = 0
    if (i == j) edge3[i,j] = 0 # self edge is excluded
  }
}
graph.eps3 = graph.adjacency(edge3, mode = "undirected")

# Components
components(graph.eps3)

# Number of edges
sum(edge3) # number of edges


# PLOTS
# Change colours
V(graph.eps1)$color[1:39] <- "red"
V(graph.eps1)$color[40:81] <- "green"

# Change colours
V(graph.eps2)$color[1:39] <- "red"
V(graph.eps2)$color[40:81] <- "green"

# Plot only nodes with more than two edges
clean <- function(graph){
  Isolated = which(degree(graph)<2)
  G2 = delete.vertices(graph, Isolated)
  return (G2)
}

par(mfrow=c(3,1))
# Epsilon = 0
plot.igraph(clean(graph.eps1), ylab = expression(epsilon == 0), vertex.color=V(graph.eps1)$colorx, vertex.label = NA)
# Epsilon = 0.5
plot.igraph(clean(graph.eps2), ylab = expression(epsilon == 0.5), vertex.color=V(graph.eps2)$color, vertex.label = NA)
# Epsilon = 0.65
plot.igraph(clean(graph.eps3), ylab = expression(epsilon == 0.65), vertex.color=V(graph.eps2)$color, vertex.label = NA)
legend("topright",
       c("right side","left side"),
       fill=c("red","green"))


# Summary  ----------------------------------------------------------------

recap <- matrix(c(sum(edge), sum(edge2), sum(edge3)), nrow = 1, ncol = 3)
colnames(recap) <- c("Eps = 0","Eps = 0.5","Eps = 0.65")
rownames(recap) <- c("Number of nodes")
recap

# Relationship between epsilon and edges ----------------------------------

# we repeat the same function for different values of epsilon
n.edges <- function(epsilon, n.features, l, u){
  e <- matrix(NA, n.features, n.features)
  for (i in 1:n.features){
    for(j in 1:n.features){
      # If 0 is not included in the CI it means there is an edge between the two nodes
      if (epsilon < l[i, j] | -epsilon > u[i, j]) e[i,j] = 1
      else e[i,j] = 0
      if (i == j) e[i,j] = 0 # self edge is excluded
    }
  }
  return (sum(e))
}

a <- seq(0,1, 0.01) # epsilon from 0 to 1

# alpha = 0.05
n.edges.vec <- list()
for (eps in a){
  n.edges1 <- n.edges(eps, n.features, ci.lower, ci.upper)
  n.edges.vec <- append(n.edges.vec, n.edges1)
}

# alpha = 0.01
n.edges.vec2 <- list()
for (eps in a){
  n.edges2 <- n.edges(eps, n.features, ci.lower2, ci.upper2)
  n.edges.vec2 <- append(n.edges.vec2, n.edges2)
}


# PLOT
par(mfrow=c(1,1))
plot(a, lwd = 4,  main="Different alphas", n.edges.vec, type = 'l',xlab = 'epsilon', ylab = 'n. edges', col = 'blue')
lines(a, n.edges.vec2, col="red", lwd = 4)
legend("topright",
       c("alpha = 0.05","alpha = 0.01"),
       fill=c("blue","red"))

# Analysis  ---------------------------------------------------------------

# We are analyzing epsilon = 0 and alpha = 0.05
edgeR <- edge [1:39, 1:39] # we get all the right-right edges 
edgeL <- edge [40:81, 40:81] # we get all the left-left edges
edgeLR <- edge [40:81, 1:39] # we get all the left-right edges


sum(edgeR) # number of nodes for right-right
sum(edgeL) # number of nodes for left-left
sum(edgeLR) # number of nodes for left-right

# We are analyzing epsilon = 0.5 and alpha = 0.05
edge2R <- edge2 [1:39, 1:39] # we get all the right-right edges 
edge2L <- edge2 [40:81, 40:81] # we get all the left-left edges
edge2LR <- edge2 [40:81, 1:39] # we get all the left-right edges


sum(edge2R) # number of nodes for right-right
sum(edge2L) # number of nodes for left-left
sum(edge2LR) # number of nodes for left-right

relation <- matrix(c(sum(edgeR), sum(edgeL), sum(edgeLR), sum(edge2R), sum(edge2L), sum(edge2LR) ), nrow = 2, ncol = 3)
colnames(relation) <- c("R-R","L-L","R-L")
rownames(relation) <- c("N. nodes, eps = 0", "N. nodes, eps = 0.5")
relation


# Partial Correlation -----------------------------------------------------

# Package
library(help = SIN) # documentation 

# Run & take a look
out <- sinUG(data.cor, n)
round(out,3) # estimated partial correlation matrix
plotUGpvalues(out) # take a look at the p-values, they are almost all at 1


# the edges grow when alpha increases
par(mfrow = c(1,1))
alpha = 0.1
E.SIN = getgraph(out, alpha)
G.SIN = graph.adjacency(E.SIN, mode = "undirected")

# Change colours
V(G.SIN)$color[1:39] <- "red"
V(G.SIN)$color[40:81] <- "green"

plot.igraph(clean(G.SIN), ylab = expression(alpha == 0.1), vertex.color=V(G.SIN)$color, vertex.label = NA)
legend("topright",
       c("right side","left side"),
       fill=c("red","green"))

# As we can see in this summary the number of edges grow when alpha grows
components(G.SIN, mode = c("weak", "strong"))
