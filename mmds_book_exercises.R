#### == MMDS - Book Exercises ==

### Exercise 3.1.1
# Compute the Jaccard similarities of each pair of the following three sets:
# {1, 2, 3, 4}, {2, 3, 5, 7}, and {2, 4, 6}.

a <- c(1, 2, 4, 5)
b <- c(2, 3, 5, 7)
c <- c(2, 4, 6)

jaccardSimPair1 <- length(intersect(a, b)) / length(union(a, b))
jaccardSimPair2 <- length(intersect(a, c)) / length(union(a, c))
jaccardSimPair3 <- length(intersect(b, c)) / length(union(b, c))

### Exercise 3.1.2
# Compute the Jaccard bag similarity of each pair of the following three bags: 
# {1, 1, 1, 2}, {1, 1, 2, 2, 3}, and {1, 2, 3, 4}.

a <- c(1, 1, 1, 2)
b <- c(1, 1, 2, 2, 3)
c <- c(1, 2, 3, 4)

# Example to calculate bagIntersect1

bagIntersect <- function(x, y) {
  intersect <- intersect(x, y)
  bagIntersect <- c()
  for (i in intersect) {
    a <- x[x == i]  
    b <- y[y == i]
    value <- min(length(a), length(b))
    bagIntersect <- c(bagIntersect, rep(i, value))
  }
  return(bagIntersect)
}

bagUnion1 <- c(a, b)
bagUnion2 <- c(a, c)
bagUnion3 <- c(c, b)

bagIntersect1 <- bagIntersect(a, b)
bagIntersect2 <- bagIntersect(a, c)
bagIntersect3 <- bagIntersect(c, b)

jaccardBagSim1 <- length(bagIntersect1) / length(bagUnion1) # 0.33
jaccardBagSim2 <- length(bagIntersect2) / length(bagUnion2) # 0.25
jaccardBagSim3 <- length(bagIntersect3) / length(bagUnion3) # 0.33

### Exercise 3.1.3
# Suppose we have a universal set U of n elements, and we choose two subsets S and T 
# at random, each with m of the n elements. What is the expected value of the Jaccard similarity of S and T ?

# 1 - Probability of being in S: m / n
# 2 - Probability of being in T: m / n
# 3 - Probability of being in S & T (commonElements): (m / n) * (m / n) = m^2 / n^2
# 4 - commonElements / (length(S) + length(T) - commonElements)
# 5 - This equals: intersect / union
# 6 - You need to subtract commonElements once, else you count them twice in the union!

### Exercise 3.2.1
# What are the first ten 3-shingles in the first sentence of Section 3.2?

sentence <- "The most effective way to represent documents as sets, for the purpose of identifying lexically similar documents is to construct from the document the set of short strings that appear within it."
shingles <- c()

for (i in 1:(nchar(sentence) - 2)) {
  shingles[i] <- substring(sentence, i, (i+2))
}

shingles[1:10]

### Exercise 3.2.2
# If we use the stop-word-based shingles of Section 3.2.4, and
# we take the stop words to be all the words of three or fewer letters, then what
# are the shingles in the first sentence of Section 3.2?

### Exercise 3.2.3
# What is the largest number of k-shingles a document of n
# bytes can have? You may assume that the size of the alphabet is large enough
# that the number of possible strings of length k is at least as n.

# Solution: n - k + 1
# Example:

n <- nchar(sentence)
k <- 3

length(shingles) # 193
n - k + 1 # 193

### Exercise 5.1.1
# Compute the PageRank of each page in Fig. 5.7, assuming no taxation.

# Number of nodes 
n = 3
# Matrix of vectors 
M <- matrix(c(1/3, 1/3, 1/3, 1/2, 0, 1/2, 1/3, 1/3, 1/3), nrow = n)
# rank vector (i.e., initialize page rank) 
r <- matrix(c(1/n, 1/n, 1/n), ncol = 1)

for (i in 1:10) {
  r <- M %*% r 
}

### Exercise 5.1.2
# Compute the PageRank of each page in Fig. 5.7, assuming β = 0.8.
b <- 0.8

# Re-initialise r:
r <- matrix(c(1/n, 1/n, 1/n), ncol = 1)

# Jumping to a random page/teleportation (1 - b):
t <- matrix(nrow = n, ncol = 1, rep(1/n, n))

# PageRank: Power Iteration (iterate until r does not change):
for (i in 1:10) {
  r <- b * M %*% r + (1 - b) * t
}

### Exercise 5.1.3
# Suppose the Web consists of a clique (set of nodes with all
# possible arcs from one to another) of n nodes and a single additional node that
# is the successor of each of the n nodes in the clique. Figure 5.8 shows this graph
# for the case n = 4. Determine the PageRank of each page, as a function of n
# and β.

b <- 0.8
n <- 5
clique <- matrix(c(0, 1/4, 1/4, 1/4, 1/4, 
                   1/4, 0, 1/4, 1/4, 1/4,
                   1/4, 1/4, 0, 1/4, 1/4,
                   1/4, 1/4, 1/4, 0, 1/4,
                   0, 0, 0, 0, 0), ncol = n)

r <- matrix(nrow = n, ncol = 1, rep(1/n, n))
t <- matrix(nrow = n, ncol = 1, rep(1/n, n))

for (i in 1:20) {
  r <- b * clique %*% r + (1 - b) * t
}

for (i in 1:20) {
  r <- b * clique %*% r
  s <- sum(r) # see what potentially "leaked out"
  r <- r + (1-s) * t # feed this back to rank-vector evenly
}

### 6.1.1

# Suppose there are 100 items, numbered 1 to 100, and also 100
# baskets, also numbered 1 to 100. Item i is in basket b if and only if i divides b
# with no remainder. Thus, item 1 is in all the baskets, item 2 is in all fifty of the
# even-numbered baskets, and so on. Basket 12 consists zof items {1, 2, 3, 4, 6, 12},
# since these are all the integers that divide 12. Answer the following questions:

### a) If the support threshold is 5, which items are frequent?

items <- 1:100
numBaskets <- 1:100
baskets <- list()
for (i in numBaskets) { baskets[[i]] <- c(0) }

for (i in numBaskets) {
  x <- 1
  for (j in items) {
    if (i %% j == 0) {
      baskets[[i]][x] <- j 
      x <- x + 1
    }
  }
}

masterSet <- c()
for (i in numBaskets) { masterSet <- c(masterSet, baskets[[i]]) }

tableMasterSet <- table(masterSet)

### b) If the support threshold is 5, which pairs of items are frequent?

masterItemSet <- c()
for (i in numBaskets) { 
  items <- baskets[[i]]
  if (length(items) > 1) {
    sets <- t(combn(items, 2))
    sets <- paste(sets[,1], sets[,2], sep = "-")
    masterItemSet <- c(masterItemSet, sets)
  }
}

tableMasterItemSet <- table(masterItemSet)
tableMasterItemSet[tableMasterItemSet >= 5] # pairs and their counts
sum(tableMasterItemSet[tableMasterItemSet >= 5]) # sum of all pairs
length(names(tableMasterItemSet[tableMasterItemSet >= 5])) # distinct count of all pairs

### c) What is the sum of the sizes of all the baskets?

listLength <- lapply(baskets, length)
sumSizes <- sum(unlist(listLength))
sumSizes

### Exercise 6.1.2 - For the item-basket data of Exercise 6.1.1, which basket is the largest?

listLength <- lapply(baskets, length)
maxLength <- max(unlist(listLength))
index <- unlist(listLength) == maxLength
numBaskets[index]

### 6.1.3 

# Suppose there are 100 items, numbered 1 to 100, and also 100
# baskets, also numbered 1 to 100. Item i is in basket b if and only if b divides i
# with no remainder. For example, basket 12 consists of items:
# 
# {12, 24, 36, 48, 60, 72, 84, 96}
# 
# Repeat Exercise 6.1.1 for this data.

items2 <- 1:100
numBaskets2 <- 1:100
baskets2 <- list()
for (i in numBaskets2) { baskets2[[i]] <- c(0) }

for (i in numBaskets2) {
  x <- 1
  for (j in items2) {
    if (j %% i == 0) { # swap j & i
      baskets2[[i]][x] <- j 
      x <- x + 1
    }
  }
}

# Continue as 6.1.1

### Exercise 6.1.4

# This question involves data from which nothing interesting
# can be learned about frequent itemsets, because there are no sets of items that
# are correlated. Suppose the items are numbered 1 to 10, and each basket is
# constructed by including item i with probability 1/i, each decision being made
# independently of all other decisions. That is, all the baskets contain item 1,
# half contain item 2, a third contain item 3, and so on. Assume the number of
# baskets is sufficiently large that the baskets collectively behave as one would
# expect statistically. Let the support threshold be 1% of the baskets. Find the
# frequent itemsets.

# Answer:
# - Say you have 100 baskets
# - 1 will be in 100% of them
# - 2 will be in 50% of them
# - 3 will be in 33%
# - ...
# -- 100 will be in 1%.
# -- 101 will be in 0.99%

### Exercise 6.1.5

# For the data of Exercise 6.1.1, what is the confidence of the following association rules?
# 6.1.5 a) {5, 7} --> 2.

a <- 0 
b <- 0

for (i in numBaskets) {
  if (sum(c(5, 7) %in% baskets[[i]]) == 2) {
    a <- a + 1
  }
  if (sum(c(2, 5, 7) %in% baskets[[i]]) == 3) {
    b <- b + 1
  }
}

confidence <- b/a
confidence # 0.5

# 6.1.5 b) {2, 3, 4} --> 5.

a <- 0 
b <- 0

for (i in numBaskets) {
  if (sum(c(2, 3, 4) %in% baskets[[i]]) == 3) {
    a <- a + 1
  }
  if (sum(c(2, 3, 4, 5) %in% baskets[[i]]) == 4) {
    b <- b + 1
  }
}

confidence <- b/a
confidence # 0.125

### Exercise 6.1.6
# For the data of Exercise 6.1.3, what is the confidence of the following association rules?

# a) {24, 60} --> 8

a <- 0 
b <- 0

for (i in numBaskets2) {
  if (sum(c(24, 60) %in% baskets2[[i]]) == 2) {
    a <- a + 1
  }
  if (sum(c(24, 60, 8) %in% baskets2[[i]]) == 3) {
    b <- b + 1
  }
}

confidence <- b/a
confidence # 0.5

# b) {2, 3, 4} --> 5

a <- 0 
b <- 0

for (i in numBaskets2) {
  if (sum(c(2, 3, 4) %in% baskets2[[i]]) == 3) {
    a <- a + 1
  }
  if (sum(c(2, 3, 4, 5) %in% baskets2[[i]]) == 4) {
    b <- b + 1
  }
}

confidence <- b/a
confidence # 1.0

## Exercise 6.1.7 : Describe all the association rules that have 100% confidence for the market-basket data of:
# a) Exercise 6.1.1.

# {2} --> {1}
# {4} --> {2, 1}
# {8} --> {4, 2, 1}
# {16} --> {8, 4, 2, 1}
# {32} --> {16, 8, 4, 2, 1} etc.

# {x} --> y 
# if you multiply x by 2, it will always inherit its previous members!
# E.g. any basket which contains 16 will also always have 8, 4, 2 & 1

# b) Exercise 6.1.3.

# Every basket contains all multiples of that specific basket #
# Hence whenever you have for instance 12 in a basket, you will always also have 24, 36, 48, etc. up to 96.

## Exercise 6.1.8
# Prove that in the data of Exercise 6.1.4 there are no interesting association rules; 
# i.e., the interest of every association rule is 0.

# Say we do {1} --> 2.
# 1 will be in 100% (1/1) of the baskets
# 2 will be in 50% (1/2) of the baskets
# {1, 2} will therefore only be in 50% of the baskets
# {2} will be in 50% of the baskets
# Interest = 1/2 - 1/2 = 0

# Say we do {2} --> 3
# {2, 3} will appear in at most 1/3 of the baskets
# {3} will appear in 1/3 of the baskets
# 1/3 - 1/3 = 0 ???
# Or {2, 3} = 1/2 * 1/3 = 1/6
# 1/6 - 1/3 = -1/6 = 0?

#### 6.2.7 Exercises for Section 6.2

### Exercise 6.2.1
# If we use a triangular matrix to count pairs, and n, the number of items, is 20, what pair’s count is in a[100]?

n <- 1:20
pairs <- as.data.frame(t(combn(n, 2)))
colnames(pairs) <- c('i', 'j')
pairs$k <- (pairs$i - 1) * (length(n) - (pairs$i/2)) + pairs$j - pairs$i
pairs[(pairs$k == 100),] # i == 7, j == 8

### ! Exercise 6.2.2
# In our description of the triangular-matrix method in Section 6.2.2, the formula for k involves dividing an arbitrary 
# integer i by 2. Yet we need to have k always be an integer. Prove that k will, in fact, be an integer.

# Answer:
# When i is an even number, i/2 will always return an integer. 
# When i is uneven, the multiplication of (i/2) with (i - 1) will make it always even again.

### ! Exercise 6.2.3
# Let there be I items in a market-basket data set of B baskets. Suppose that every basket contains exactly K items. 
# As a function of I, B, and K:

# Example:

I <- 100000
B <- 10000000
K <- 10

# a) How much space does the triangular-matrix method take to store the counts of all pairs of items, 
# assuming four bytes per array element?

space <- 4*((I^2)/2) # bytes
spaceMb <- space / 1024 / 1024 # bytes -> kilobytes -> megabytes

# b) What is the largest possible number of pairs with a nonzero count?

numPairs <- B * ((K * (K-1))/2)

# c) Under what circumstances can we be certain that the triples method will use less space than the triangular array?

# Answer:
# If the real number of frequent pairs is likely to be much smaller then the total number of possible pairs.

### !! Exercise 6.2.4
# How would you count all itemsets of size 3 by a generalization of the triangular-matrix method? 
# That is, arrange that in a one-dimensional array there is exactly one element for each set of three items.

NA

### ! Exercise 6.2.5
# Suppose the support threshold is 5. Find the maximal frequent itemsets for the data of:
# a) Exercise 6.1.1.

n <- length(unique(unlist(baskets)))
possiblePairs <- (n * (n - 1)) / 2
array <- c(rep(0, possiblePairs))

for (i in numBaskets) { 
  items <- baskets[[i]]
  if (length(items) > 1) {
    sets <- t(combn(items, 2))
    for (i in 1:dim(sets)[1]) {
      k <- (sets[i,1] - 1) * (n - (sets[i,1]/2)) + sets[i,2] - sets[i,1]
      array[k] <- array[k] + 1
    }
  }
}

lookupTable <- as.data.frame(t(combn(1:n, 2)))
colnames(lookupTable) <- c('i', 'j')
lookupTable$k <- (lookupTable$i - 1) * (n - (lookupTable$i/2)) + lookupTable$j - lookupTable$i

supportIndex <- array >= 5
frequentPairs <- lookupTable[supportIndex, ]

# b) Exercise 6.1.3.

# Same as a), swap baskets with baskets2

### Exercise 6.2.6
# Apply the A-Priori Algorithm with support threshold 5 to the data of:
# a) Exercise 6.1.1.

## Support Functions:

library(plyr)

countGenerator <- function(sourceTable, candidateSets, relevantInts, subsets) {
  for (i in 1:length(sourceTable)) {
    basket <- sourceTable[[i]]
    basket <- basket[basket %in% relevantInts]
    if (length(basket) > subsets) {
      combinations <- t(combn(basket, subsets))
      for (j in 1:dim(combinations)[1]) {
        candidate <- as.data.frame(t(combinations[j,]))
        index <- rownames(candidateSets) == suppressMessages(rownames(match_df(candidateSets, candidate)))
        candidateSets[index, 'count'] <- candidateSets[index, 'count'] + 1
      }
    }
  }
  return(candidateSets)
}

candidateGenerator <- function(frequentTable) {
  possibleSets <- as.data.frame(t(combn(unique(unlist(frequentTable[,1:dim(frequentTable)[2]-1])), dim(frequentTable)[2])))
  possibleSets$candidate <- 0
  for (i in 1:dim(possibleSets)[1]) {
    subSets <- t(combn(possibleSets[i,(1:dim(possibleSets)[2]-1)], (dim(frequentTable)[2]-1)))
    subSets <- matrix(as.matrix(sapply(subSets, as.numeric)), ncol = (dim(frequentTable)[2]-1))
    x <- 0
    for (j in 1:dim(subSets)[1]) {
      y <- as.data.frame(t(subSets[j,]))
      index <- rownames(frequentTable) == suppressMessages(rownames(match_df(frequentTable, y)))
      if (sum(index) > 0) {
        x <- x + max(index * 1)
      }
    }
    if (x == dim(subSets)[1]) { possibleSets[i, "candidate"] <- 1 }
  }
  return(possibleSets)
}

# 1) A-Priori - Pass 1
# Generating frequent integers

integers <- c()
for (i in numBaskets) { integers <- c(integers, baskets[[i]]) }

intArray <- c(rep(0, 100))
for (i in integers) { intArray[i] <- intArray[i] + 1}

# In-Between Passes
x <- 1:100
index <- intArray >= 5
freqItemTable <- c(rep(0, 100))
freqItemTable[index] <- x[intArray >= 5]

# 2) A-Priori - Pass 2
# Generating frequent pairs

masterSet <- matrix(nrow = 1, ncol = 2)

for (i in numBaskets) {
  set <- baskets[[i]]
  freqSingletons <- set[set %in% freqItemTable]
  if (length(freqSingletons) > 1) {
    sets <- t(combn(freqSingletons, 2))
    masterSet <- rbind(masterSet, sets)
  }
}

df <- as.data.frame(masterSet)
df <- count(df, vars=c('V1','V2'))
df <- df[(complete.cases(df) == TRUE),]
dfFrequent <- df[(df$freq >= 5),]

# 3.1) Checking which of possible triples are candidates (only those were all sub-sets are frequent):

candidateTriples <- candidateGenerator(dfFrequent)
candidateTriples <- candidateTriples[(candidateTriples$candidate == 1), 1:3]

# 3.2) Making another pass through the data, collecting counts of candidate triples:

candidateTriples$count <- 0
tripleInts <- unique(unlist(candidateTriples[,1:3]))
candidateTriples <- countGenerator(baskets, candidateTriples, tripleInts, 3)
frequentTriples <- candidateTriples[(candidateTriples$count >= 5),]

# 4.1) Checking which of possible quadruples are candidates (only those were all sub-sets are frequent):

candidateQuadruples <- candidateGenerator(frequentTriples)
candidateQuadruples <- candidateQuadruples[(candidateQuadruples$candidate == 1), 1:4]

# 4.2) Making another pass through the data, collecting counts of candidate quadruples:

candidateQuadruples$count <- 0
quadrupleInts <- unique(unlist(candidateQuadruples[,1:4]))
candidateQuadruples <- countGenerator(baskets, candidateQuadruples, quadrupleInts, 4)
frequentQuadruples <- candidateQuadruples[(candidateQuadruples$count >= 5),]

# 5.1) Generating all possible quintuples:

candidateQuintuples <- candidateGenerator(frequentQuadruples)
candidateQuintuples <- candidateQuintuples[(candidateQuintuples$candidate == 1), 1:5]

# 5.2) Making another pass through the data, collecting counts of candidate quintuples:

candidateQuintuples$count <- 0
quintupleInts <- unique(unlist(candidateQuintuples[,1:5]))
candidateQuintuples <- countGenerator(baskets, candidateQuintuples, quintupleInts, 5)
frequentQuintuples <- candidateQuintuples[(candidateQuintuples$count >= 5),]

# 6.1) Generating all possible sextuples:

candidateSextuples <- candidateGenerator(frequentQuintuples)
candidateSextuples <- candidateSextuples[(candidateSextuples$candidate == 1), 1:6]

# 6.2) Making another pass through the data, collecting counts of candidate sextuples:

candidateSextuples$count <- 0
quintupleInts <- unique(unlist(candidateSextuples[,1:6]))
candidateSextuples <- countGenerator(baskets, candidateSextuples, quintupleInts, 6)
frequentSextuples <- candidateSextuples[(candidateSextuples$count >= 5),]

# None of the candidate sextuples are frequent!

# b) Exercise 6.1.3.

# Same as 6.1.2 but swap baskets with baskets2.