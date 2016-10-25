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

## 6.2.7
# Suppose we have market baskets that satisfy the following assumptions:
#
# 1. The support threshold is 10,000.
# 2. There are one million items, represented by the integers 0, 1, . . . , 999999.
# 3. There are N frequent items, that is, items that occur 10,000 times or more.
# 4. There are one million pairs that occur 10,000 times or more.
# 5. There are 2M pairs that occur exactly once. Of these pairs, M consist of
#    two frequent items; the other M each have at least one nonfrequent item.
# 6. No other pairs occur at all.
# 7. Integers are always represented by 4 bytes.
#
# Suppose we run the A-Priori Algorithm and can choose on the second pass
# between the triangular-matrix method for counting candidate pairs and a hash
# table of item-item-count triples. Neglect in the first case the space needed to
# translate between original item numbers and numbers for the frequent items,
# and in the second case neglect the space needed for the hash table. As a function
# of N and M, what is the minimum number of bytes of main memory needed to
# execute the A-Priori Algorithm on this data?

# Case 1 - Triangular Matrix
# Neglect: Item names to integers

triangularBytes <- function(N) {
  freqItemsBytes <- N * 4 # not counts, only items (!= * 2)
  triangularArraySize <- (N * (N - 1)) / 2
  triangularArrayBytes <- triangularArraySize * 4
  bytes <- freqItemsBytes + triangularArrayBytes
  return(bytes)
}

# Case 2 - Triples Method
# Neglect: hash table (Bitmap)

tripleBytes <- function(M) {
  freqItemsBytes <- N * 4 # not counts, only items (!= * 2)
  frequentPairs <- 1000000 # pairs whose items are both frequent and are frequent itself
  nonFrequentPairs <- M # pairs whose items are both frequent, yet are not frequent itself
  totalCounts <- frequentPairs + nonFrequentPairs
  bytes <- (totalCounts * 3 * 4) + freqItemsBytes # we need to store 3 integers of each 4 bytes
  return(bytes)
}

# Note: I am not sure if we do not also need to neglect "Item names to integers" for Case 2 (which I do).

### Exercise 6.3.1
# Here is a collection of twelve baskets. Each contains three of the six items 1 through 6.

b1 <- c(1,2,3); b2 <- c(2,3,4); b3 <- c(3,4,5); b4 <- c(4,5,6)
b5 <- c(1,3,5); b6 <- c(2,4,6); b7 <- c(1,3,4); b8 <- c(2,4,5)
b9 <- c(3,5,6); b10 <- c(1,2,4); b11 <- c(2,3,5); b12 <- c(3,4,6)

baskets <- list(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12)
rm(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12)

# Suppose the support threshold is 4. On the first pass of the PCY Algorithm
# we use a hash table with 11 buckets, and the set {i, j} is hashed to bucket i * j mod 11.

# a) By any method, compute the support for each item and each pair of items.

s <- 4 # support threshold
items <- 1:6
itemsCount <- vector(length = 6)
buckets <- vector(length = 11)

# PCY - Pass 1

for (basket in baskets) {
  
  # This loop counts all the items
  for (item in basket) {
    itemsCount[item] <- itemsCount[item] + 1
  }
  
  pairs <- t(combn(basket, 2))
  
  # This loop hashes the pairs into buckets
  for (x in 1:dim(pairs)[1]) {
    i <- pairs[x,1]
    j <- pairs[x,2]
    bucket <- (i * j) %% 11
    buckets[bucket] <- buckets[bucket] + 1
  }
}

# PCY - Between Passes

frequentItems <- items[itemsCount >= s] # all of them are frequent
bitmap <- (buckets >= s) * 1

# PCY - Pass 2

countedPairs <- data.frame(x = 0, y = 0, count = 0)

for (basket in baskets) {
  
  # Only looking at items that are frequent
  basket <- basket[basket %in% frequentItems]
  
  pairs <- t(combn(basket, 2))
  
  for (z in 1:dim(pairs)[1]) {
    i <- pairs[z,1]
    j <- pairs[z,2]
    bucket <- (i * j) %% 11
    
    # Only keeping counts if the pair hashes to a frequent bucket
    if (bitmap[bucket] == 1) {
      
      # Using the triples method to store counts  
      if (dim(countedPairs[(countedPairs$x == i & countedPairs$y == j),])[1] == 0) {
        countedPairs <- rbind(countedPairs, c(i,j,1))
      } else {
        countedPairs[(countedPairs$x == i & countedPairs$y == j), 3] <- 
          countedPairs[(countedPairs$x == i & countedPairs$y == j), 3] + 1
      }
    } 
  }
}

countedPairs <- countedPairs[2:dim(countedPairs)[1],]

# Final Results

frequentPairs <- countedPairs[(countedPairs$count >= s),] # 2-4, 3-4, 3-5

# b) Which pairs hash to which buckets?

bucketList <- list()

for (i in 1:11) {
  bucketList[[i]] <- c(i)
}

for (basket in baskets) {
  
  pairs <- t(combn(basket, 2))
  
  # This loop hashes the pairs into buckets
  for (x in 1:dim(pairs)[1]) {
    i <- pairs[x,1]
    j <- pairs[x,2]
    bucket <- (i * j) %% 11
    bucketList[[bucket]] <- c(bucketList[[bucket]], c(i,j))
  }
}

bucketPairs <- function(vec) {
  bucketNum <- vec[1]
  pairs <- vec[-1]
  pairsStr <- vector(length = length(pairs)/2)
  count <- 1
  
  if (length(pairs) > 1) {
    for (x in seq(1, length(pairs), 2)) {
      i <- pairs[x]
      j <- pairs[x+1]
      pairsStr[count] <- paste0(as.character(i), '-', as.character(j))
      count <- count + 1
    }
    writeLines(paste0("Bucket: ", as.character(bucketNum)))
    writeLines(unique(pairsStr))
  }
}

for (i in bucketList) {
  bucketPairs(i)
}

# c) Which buckets are frequent?

numBuckets <- 1:11
numBuckets[buckets > s] # 1, 2, 4, 8

# d) Which pairs are counted on the second pass of the PCY Algorithm?

countedPairs # 1-2, 2-4, 3-4, 3-5, 4-6, 5-6, 2-6, 1-4

## Exercise 6.3.2
# Suppose we run the Multistage Algorithm on the data of Exercise 6.3.1, with the same support 
# threshold of 4. The first pass is the same as in that exercise, and for the second pass, 
# we hash pairs to nine buckets, using the hash function that hashes {i, j} to bucket i + j mod 9. 
# Determine the counts of the buckets on the second pass. Does the second pass reduce the
# set of candidate pairs? Note that all items are frequent, so the only reason a
# pair would not be hashed on the second pass is if it hashed to an infrequent bucket on the first pass.

s <- 4 # support threshold
items <- 1:6
itemsCount <- vector(length = 6)
buckets <- vector(length = 11)

# Multistage - Pass 1

for (basket in baskets) {
  
  # This loop counts all the items
  for (item in basket) {
    itemsCount[item] <- itemsCount[item] + 1
  }
  
  pairs <- t(combn(basket, 2))
  
  # This loop hashes the pairs into buckets
  for (x in 1:dim(pairs)[1]) {
    i <- pairs[x,1]
    j <- pairs[x,2]
    bucket <- (i * j) %% 11
    buckets[bucket] <- buckets[bucket] + 1
  }
}

# PCY - Between Passes

frequentItems <- items[itemsCount >= s] # all of them are frequent
bitmap <- (buckets >= s) * 1

# Multistage - Pass 2

buckets2 <- vector(length = 9)

for (basket in baskets) {
  
  # Only looking at items that are frequent
  basket <- basket[basket %in% frequentItems]
  
  pairs <- t(combn(basket, 2))
  
  for (z in 1:dim(pairs)[1]) {
    i <- pairs[z,1]
    j <- pairs[z,2]
    bucket <- (i * j) %% 11
    
    # Only re-hashing pairs that hashed to a frequent bucket in Pass 1
    if (bitmap[bucket] == 1) {
      bucket2 <- (i * j) %% 9
      buckets2[bucket2] <- buckets2[bucket2] + 1
    }
  }
}

bitmap2 <- (buckets2 >= s) * 1

sum(buckets) # 36
sum(buckets2) # 22 

# A: Yes it does, by 14.

## Exercise 6.3.3
# Suppose we run the Multihash Algorithm on the data of Exercise 6.3.1. 
# We shall use two hash tables with five buckets each. For one, the set {i, j}, is hashed to bucket 2i+3j+4 mod 5, 
# and for the other, the set is hashed to i+4j mod 5. Since these hash functions are not symmetric in i and j, 
# order the items so that i < j when evaluating each hash function. Determine the counts of each of the 10 buckets. 
# How large does the support threshold have to be for the Multistage Algorithm to eliminate more pairs than the PCY
# Algorithm would, using the hash table and function described in Exercise 6.3.1?
# Personal Note: I assume they mean Multihash in l827

# Multihash - Pass 1

# s <- 4 # support threshold
items <- 1:6
itemsCount <- vector(length = 6)
bucketsA <- vector(length = 5)
bucketsB <- vector(length = 5)

for (basket in baskets) {
  
  # This loop counts all the items
  for (item in basket) {
    itemsCount[item] <- itemsCount[item] + 1
  }
  
  # Pairs are automatically generated so that i < j
  pairs <- t(combn(basket, 2))
  
  for (z in 1:dim(pairs)[1]) {
    i <- pairs[z,1]
    j <- pairs[z,2]
    bucketA <- (2*i + 3*j + 4) %% 5
    bucketB <- (i + 4*j) %% 5
    bucketsA[bucketA] <- bucketsA[bucketA] + 1
    bucketsB[bucketB] <- bucketsB[bucketB] + 1
  }
}

bitmapA <- (bucketsA >= s) * 1
bitmapB <- (bucketsB >= s) * 1

# Counts of the buckets:
bucketsA # 2 14 6 0 0 | sum == 22
bucketsB # 2 6 14 14 0 | sum == 36

# PCY had 22 pairs hashed into buckets with a support threshold of 4.
# If we want to eliminate more pairs, we should remain with fewer pairs. 
# If we set s <- 3, we will for sure have fewer pairs:

sum(bucketsA[bucketsA >= 3]) # 20 <-- upper limit
sum(bucketsB[bucketsB >= 3]) # 34

# bucketsA will let max 20 pairs through, which can only be when all pairs that hash into a 
# frequent bucket in bucketsA also always hashes into a frequent bucket in bucketsB. This is however unlikely: 
# in several cases a pair that hashes to a frequent bucket in bucketsA will hash into an infrequent bucket in bucketsB.
# Therefore it is likely we end up with less than 20. We could also do the exercise and simply check how many candidate pairs
# both algorithms ultimately generate, for varying levels of s.

## ! 6.3.4
# Suppose we perform the PCY algorithm to find frequent pairs, 
# with market-basket data meeting the following specifications:
#
# - s, the support threshold, is 10,000.
# - There are one million items, which are represented by the integers 0,1,...,999999.
# - There are 250,000 frequent items, that is, items that occur 10,000 times or more.
# - There are one million pairs that occur 10,000 times or more.
# - There are P pairs that occur exactly once and consist of 2 frequent items.
# - No other pairs occur at all.
# - Integers are always represented by 4 bytes.
#
# When we hash pairs, they distribute among buckets randomly, but as evenly as possible; 
# i.e., you may assume that each bucket gets exactly its fair share of the P pairs that occur once.
# Suppose there are S bytes of main memory. In order to run the PCY algorithm successfully, 
# the number of buckets must be sufficiently large that most buckets are not large. 
# In addition, on the second pass, there must be enough room to count all the candidate pairs. 
# As a function of S, what is the largest value of P for which we can successfully run the PCY algorithm on this data? 

NA

## ! Exercise 6.3.5
# Under the assumptions given in Exercise 6.3.4, will the Multihash
# Algorithm reduce the main-memory requirements for the second pass?
# As a function of S and P, what is the optimum number of hash tables to use on the first pass?

NA

## ! Exercise 6.3.6
# Suppose we perform the 3-pass Multistage Algorithm to find frequent pairs, with market-basket 
# data meeting the following specifications:
# 
# 1. The support threshold is 10,000.
# 2. There are one million items, represented by the integers 0, 1, . . . , 999999.
#    All items are frequent; that is, they occur at least 10,000 times.
# 3. There are one million pairs that occur 10,000 times or more.
# 4. There are P pairs that occur exactly once.
# 5. No other pairs occur at all.
# 6. Integers are always represented by 4 bytes.
# 7. When we hash pairs, they distribute among buckets randomly, but as
# evenly as possible; i.e., you may assume that each bucket gets exactly its
# fair share of the P pairs that occur once.
# 8. The hash functions used on the first two passes are completely independent.
# 
# Suppose there are S bytes of main memory. As a function of S and P, what
# is the expected number of candidate pairs on the third pass of the Multistage Algorithm?

NA

## Exercise 6.4.1
# Suppose there are eight items, A, B, .., H, and the following
# are the maximal frequent itemsets: {A,B}, {B,C}, {A,C}, {A,D}, {E}, and
# {F}. Find the negative border.

freqItemSets <- list(c('A','B'), c('B','C'), c('A','C'), c('A','D'), c('E'))
itemsChr <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')

# 1 - Items that are not present in any of the maximal frequent itemsets are in the negative border (because their subsets [the empty set {}] are frequent).

freqItems <- unique(unlist(freqItemSets))
negBorderItems <- itemsChr[!itemsChr %in% freqItems]

# 2 - Get all possible pairs of frequent items, if they are not in the maximal frequent itemsets, they are in the negative border.

pairs <- t(combn(freqItems, 2))
pairs <- cbind(pairs, 1)

for (i in 1:dim(pairs)[1]) {
  vector <- c(pairs[i,1], pairs[i,2])
  for (j in freqItemSets) {
    logical <- vector == j
    if (sum(logical) == 2) {
      pairs[i,3] <- '0'
    }
  }
}

negBorderPairs <- pairs[(pairs[,3] == 1), 1:2]

# 3 - Get all possible triples of items that are in a frequent pair, if all sub-pairs are frequent they are in the negative border.
#     - Because there are no triples among the maximal frequent itemsets, if we find any they must be part of the negative border.

# 3.1) Generating candidate triples (all subsets must be frequent, hence the items use to expand pairs to triples must be too)

# A B -> C | D | E
# B C -> D | E 
# A C -> D | E 
# A D -> E

candidateTriples <- matrix(ncol = 3, nrow = 0)

for (i in freqItemSets) {
  if (length(i) == 2) {
    logical <- (freqItems %in% i) * 1
    logical[which.max(logical)] <- 0
    appendItems <- tail(freqItems, -which.max(logical))
    for (j in appendItems) {
      candidateTriples <- rbind(candidateTriples, c(i, j))
    }
  }
}

candidateTriples <- cbind(candidateTriples, 0)

# 3.2) Final triples: checking if all subpairs are frequent, if so, the triple is in the negative border.

for (i in 1:dim(candidateTriples)[1]) {
  pairs <- t(combn(candidateTriples[i,1:3],2))
  count <- 0
  for (j in 1:dim(pairs)[1]) {
    vec <- pairs[j,]
    for (k in freqItemSets) {
      if (length(k) == 2) {
        # or: all(vec %in% k)
        if (sum(vec %in% k) == 2) {
          count <- count + 1
        }
      }
    }
  }
  if (count == 3) {
    candidateTriples[i,4] <- 1
  }
}

negBorderTriples <- candidateTriples[(candidateTriples[,4] == 1), 1:3]

# A:
negBorderItems # F, G, H
negBorderPairs # A-E, B-D, B-E, C-D, C-E, D-E
negBorderTriples # A-B-C

## Exercise 6.4.2
# Apply Toivonen’s Algorithm to the data of Exercise 6.3.1,
# with a support threshold of 4. Take as the sample the first row of baskets:
# {1, 2, 3}, {2, 3, 4}, {3, 4, 5}, and {4, 5, 6}, i.e., one-third of the file.

sample <- list(c(1,2,3), c(2,3,4), c(3,4,5), c(4,5,6))

# Our scaled-down support theshold will be 1. 
# Personal note: Original s was 4. You scale it down more than 1/3 of 4 (1.33) to reduce the risk of False-Negatives, 
# by taking a value that is slightly lower: 1.
 
# a) What are the itemsets frequent in the sample?

# A: Since our support threshold is 1, basically all basket triples are also the maximum frequent itemsets!

# b) What is the negative border?

# b1) - Get all frequent items, compute possible pairs -> which have occurred? The ones that have not
# occurred are in the negative border.

freqItems <- unique(unlist(sample)) # 1 2 3 4 5 6

possiblePairs <- cbind(t(combn(freqItems, 2)), 0)

for (i in sample) {
  if (length(i) == 3) {
    pairs <- t(combn(i, 2))
    for (j in 1:dim(pairs)[1]) {
      a <- pairs[j,1]
      b <- pairs[j,2]
      possiblePairs[(possiblePairs[,1] == a & possiblePairs[,2] == b), 3] <- 1
    }
  }
}

negBorderPairs <- possiblePairs[(possiblePairs[,3] == 0), 1:2] # 1-4, 1-5, 1-6, 2-5, 2-6, 3-6

# b2) - Take the frequent pairs, see which triples they can generate, those that did not occur are in the negative border.

freqPairs <- possiblePairs[(possiblePairs[,3] == 1), 1:2] # 1-2, 1-3, 2-3, 2-4, 3-4, 3-5, 4-5, 4-6, 5-6

candidateGen <- function(subSets, freqItems) {
  candidates <- matrix(ncol = (dim(subSets)[2] + 1), nrow = 0)
  for (i in 1:dim(subSets)[1]) {
    vec <- c(subSets[i,])
    logical <- freqItems %in% vec * 1
    x <- 0
    repeat {
      logical[which.max(logical)] <- 0
      x <- x + 1
      if (x == (dim(subSets)[2] - 1)) { break }
    }
    appendItems <- tail(freqItems, -which.max(logical))
    for (j in appendItems) {
      candidates <- rbind(candidates, c(vec, j))
    }
  }
  return(cbind(candidates, 0))
}  

candidateTriples <- candidateGen(freqPairs, freqItems)

for (i in sample) {
  candidateTriples[(candidateTriples[,1] == i[1] & 
                      candidateTriples[,2] == i[2] & 
                      candidateTriples[,3] == i[3]), 4] <- 1
}

freqTriples <- candidateTriples[(candidateTriples[,4] == 1), 1:3] # 1-2-3, 2-3-4, 3-4-5, 4-5-6
negBorderTriples <- candidateTriples[(candidateTriples[,4] == 0), 1:3] # 1-2-4, 1-2-5, 1-2-6, ...

# b3) - Check if there any possible quadruples in the negative border.

freqTriples <- matrix(ncol = 3, nrow = 0)
for (i in sample) {
  freqTriples <- rbind(freqTriples, i)
}

candidateQuadruples <- candidateGen(freqTriples, freqItems)

freqItemSets <- sample

for (i in 1:dim(candidateQuadruples)[1]) {
  triples <- t(combn(candidateQuadruples[i,1:4],3))
  count <- 0
  for (j in 1:dim(triples)[1]) {
    vec <- triples[j,]
    for (k in freqItemSets) {
      if (length(k) == 3) {
        if (sum(vec %in% k) == 3) {
          count <- count + 1
        }
      }
    }
  }
  if (count == 4) {
    candidateQuadruples[i,5] <- 1
  }
}

# There are no quadruples in the negative border.
negBorderQuadruples <- candidateQuadruples[(candidateQuadruples[,5] == 1), 1:4] # None

# c) What is the outcome of the pass through the full dataset? Are any of the
#    itemsets in the negative border frequent in the whole?

negBorderTriples <- cbind(0, negBorderTriples)
negBorderPairs <- cbind(0, negBorderPairs)

for (basket in baskets) {
  
  for (row in 1:dim(negBorderTriples)[1]) {
    if (all(negBorderTriples[row, -1] %in% basket)) {
      negBorderTriples[row, 1] <- negBorderTriples[row, 1] + 1
    }
  }
  
  for (row in 1:dim(negBorderPairs)[1]) {
    if (all(negBorderPairs[row, -1] %in% basket)) {
      negBorderPairs[row, 1] <- negBorderPairs[row, 1] + 1
    }
  }
}

s <- 4
# No items in the negative border were frequent in the whole
negBorderTriples[(negBorderTriples[,1] >= s), -1]
negBorderPairs[(negBorderPairs[,1] >= s), -1]

# Are the sample frequent items also frequent in the whole?
itemCounts <- table(unlist(baskets))
freqItems <- itemCounts[itemCounts >= s] # all items are frequent

freqPairs <- cbind(0, freqPairs)
freqTriples <- cbind(0, freqTriples)

for (basket in baskets) {
  
  for (row in 1:dim(freqPairs)[1]) {
    if (all(freqPairs[row, -1] %in% basket)) {
      freqPairs[row, 1] <- freqPairs[row, 1] + 1
    }
  }
  
  for (row in 1:dim(freqTriples)[1]) {
    if (all(freqTriples[row, -1] %in% basket)) {
      freqTriples[row, 1] <- freqTriples[row, 1] + 1
    }
  }
}

freqPairsWhole <- freqPairs[(freqPairs[,1] >= s), 2:3] # 2-4, 3-4, 3-5
freqTriplesWhole <- freqTriples[(freqTriples[,1] >= s), 2:4] # None

## !! Exercise 6.4.3
# Suppose item i appears exactly s times in a file of n baskets,
# where s is the support threshold. If we take a sample of n/100 baskets, and
# lower the support threshold for the sample to s/100, what is the probability
# that i will be found to be frequent? You may assume that both s and n are divisible by 100.

NA

### Chapter 12.4 - Nearest Neighbours

## Example 12.13
data <- matrix(c(1,1,2,2,4,3,8,4,16,5,32,6), ncol = 2, byrow = TRUE)
q <- 3.5

# Weight decays as the square of the distance
weight <- function(x, q) {
  value <- 1 / (x - q)^2
  return(value)
}

data <- cbind(data, weight(data[,1], q))
data <- cbind(data, data[,2] * data[,3])
sum(data[,4]) / sum(data[,3]) # 5.514039

# A popular choice is to use a normal distribution
weightBellCurve <- function(x, q, s) {
  # s == stdev of the distribution
  # q == acts as the mean (how far from this 'mean' are the neighbors)
  value <- exp(1)^-((x-q)^2)/s^2
  return(value)
}

x <- data[,1]
s <- sd(x)
q <- 3

weightBellCurve(x, q, s)[2:4]

## Exercise 12.4.2
# Suppose we have the following training set:

trainingSet <- matrix(c(1,2,1,3,4,-1,2,1,-1,4,3,1), ncol = 3, byrow = TRUE)

# which is the training set used in Example 12.9. If we use nearest-neighbor
# learning with the single nearest neighbor as the estimate of the label of a query
# point, which query points are labeled +1?

# http://flowingdata.com/2016/04/12/voronoi-diagram-and-delaunay-triangulation-in-r/
# install.packages("deldir", dependencies=TRUE)
library(deldir)

x <- trainingSet[,1]
y <- trainingSet[,2] 
vtess <- deldir(x, y)

plot(x, y, type="n", asp=1)
points(x, y, pch=20, col="red", cex=0.5)
plot(vtess, wlines="tess", wpoints="none", number=FALSE, add=TRUE, lty=1)
text(trainingSet[,1], trainingSet[,2], trainingSet[,3], cex=0.6, pos=1, col="red")

# A: View plot for decision boundaries.

## Exercise 12.4.3
# Consider the one-dimensional training set

data <- matrix(c(1,1,2,2,4,3,8,4,16,5,32,6), ncol = 2, byrow = TRUE)

# Describe the function f(q), the label that is returned in response to the query
# q, when the interpolation used is:

# Rule 1 - Nearest Neighbor.
# a) The label of the nearest neighbor:

nearestNeighbor1 <- function(query, data) {
  vec <- vector(length = dim(data)[1])
  for (i in 1:dim(data)[1]) {
    distance <- abs(data[i,1] - query)
    vec[i] <- distance
  }
  index <- which.min(vec)
  return(data[index,2])
}

nearestNeighbor1(3, data)

# Rule 2 - Average of Two Nearest Neighbors.
# b) The average of the labels of the two nearest neighbors:

nearestNeighbor2 <- function(query, data, k) {
  vec <- vector(length = dim(data)[1])
  for (i in 1:dim(data)[1]) {
    distance <- abs(data[i,1] - query)
    vec[i] <- distance
  }
  indices <- vector(length = k)
  for (j in 1:k) {
    indices[j] <- which.min(vec)
    vec <- vec[-indices[j]]
  }
  labelsWeighted <- (data[indices,2] * 1/k) # weighted by 1/k
  return(sum(labelsWeighted))
}

nearestNeighbor2(3, data, 2)

# Rule 3 - Weighted Average of the Two Nearest Neighbors.
# c) The average, weighted by distance, of the two nearest neighbors:

nearestNeighbor3 <- function(query, data, k) {
  vec <- vector(length = dim(data)[1])
  for (i in 1:dim(data)[1]) {
    distance <- abs(data[i,1] - query)
    vec[i] <- distance
  }
  indices <- vector(length = k)
  weights <- vector(length = k)
  for (j in 1:k) {
    indices[j] <- which.min(vec)
    weights[j] <- 1 / vec[indices[j]]
    vec <- vec[-indices[j]]
  }
  if (sum(weights == Inf) > 0) {
    label <- data[indices[weights == Inf], 2]
    return(label)
  } else {
    labelsWeighted <- (data[indices,2] * weights)
    return(sum(labelsWeighted) / sum(weights))
  }
}

nearestNeighbor3(32, data, 2)

# Note: 
# When following the formula, when there is an exact match, distance becomes 1 / 0 --> Inf, which results in NaN; 
# theoretically its weight is infinite, hence it dominates and should reflect this label, which is what you want.

# Rule 4 - Average of Three Nearest Neighbors.
# d) The (unweighted) average of the three nearest neighbors:

nearestNeighbor2(3, data, 3)

## ! Exercise 12.4.4
# Apply the kernel function of Example 12.13 to the data of
# Exercise 12.4.3. For queries q in the range 2 < q < 4, what is the label of q?

q <- seq(2, 4, 0.1) # Technically 2 < q < 4 means q == 3
for (i in q) {
  index <- which.max(weight(data[,1], i))
  print(paste0('q: ', as.character(i), ' | Label: ', as.character(data[index,2])))
}

# A: From 2.1-2.9 it is 2, at 3 it is tied between 2 & 3 and from 3.1-3.9 it is 3.

## Exercise 12.4.5
# What is the function that estimates the label of query points
# using the data of Example 12.12 and the average of the four nearest neighbors?

data2 <- matrix(c(1,1,2,2,3,4,4,8,5,4,6,2,7,1), ncol = 2, byrow = TRUE)
q <- 3
k <- 4
nearestNeighbor2(q, data2, 4) # 2.25

# !! Exercise 12.4.6
# Simple weighting functions such as those in Example 12.12
# need not define a continuous function. We can see that the constructed functions
# in Fig. 12.22 and Fig. 12.23(b) are not continuous, but Fig. 12.23(a) is. Does the
# weighted average of two nearest neighbors always give a continuous function?

NA

### Chapter 10

## Video Lecture Example - From AGM to BigCLAM

Fu <- c(0,1.2,0,0.2)
Fv <- c(0.5,0,0,0.8)
Fw <- c(0,1.8,1,0)

# %*% Generates the dot-product between 2 vectors
Puv <- 1 - exp(-Fu %*% Fv) # 0.14
Puw <- 1 - exp(-Fu %*% Fw) # 0.88
Pvw <- 1 - exp(-Fv %*% Fw) # 0

## 10.5.2) MLE - Maximum Likelihood Estimation

# Consider a graph with 15 nodes and 23 edges. 

n <- 15
e <- 23
pairs <- (n * (n - 1)) / 2

# Each edge is independently chosen with probability p.
# Then, the probability of generating exactly this graph is: 
prob <- function(p, e, pairs) {
  return(p^e*(1-p)^(pairs - e))
} 
# Or: p^23*(1 - p)^82
# But we do not know what p is.

# However, we do know that this function has a maximum, a value for p that returns the maximum value.
# We can determine this by taking its derivate, and setting this to 0.

# Original:
# p^23 * (1 - p)^82
# Derivate:
# 23*p^22*(1-p)^82 - 82p^23(1-p)81 = 0
# Derivative re-written:
# p^22(1 - p)^81 * (23(1 - p) - 82p) = 0
# The only way the right side = 0 if if p is 0 or 1, or if (23(1 - p) - 82p) is 0.

# When p is 0 or 1, the value of the likelihood function p^23*(1 − p)^82 is
# minimized, not maximized: 
test <- function(p) {
  return(p^22*(1-p)^81 * (23*(1-p)-82*p))
}
test(1) # 0
test(0) # 0

# So it must be the last factor that is 0: | Note: Why not p^22*(1-p)^81 = 0?
# 23 − 23p − 82p = 0

p <- 23/105 # should come as no surprise

# This is the value of p that maximizes the occurence of the graph we 
# found, with e edges and n nodes.
prob(p, e, pairs) # 1.065614e-24
# Any other value for p will return a smaller probability.
q <- 24/105
r <- 22/105
prob(p, e, pairs) > prob(q, e, pairs) # TRUE
prob(p, e, pairs) > prob(r, e, pairs) # TRUE

### Example - AGM
## Model B(V, C, M, {Pc})
# V - Nodes 
V <- 1:20

# C - Communities 
C <- list()

# Pc - probability of connect per community
Pc <- c(0.2, 0.4, 0.6)

x <- 0
for (i in Pc) {
  x <- x + 1
  size <- i * length(V)
  C[[x]] <- sample(V, size, replace = FALSE)
}

# M - Memberships
M <- cbind(0, t(combn(V, 2)))

# e - Default probability (no shared communities)
e <- 0.001

# I - C as List Implementation:

for (i in 1:dim(M)[1]) { 
  pair <- M[i,-1]
  sharedC <- rep(0, length(C))
  count <- 1
  for (j in C) {
    if (all(pair %in% j) == TRUE) {
      sharedC[count] <- 1
      count <- count + 1
    } else {
      count <- count + 1
    }
  }
  if (sum(sharedC) == 0) {
    M[i,1] <- e
  } else {
    M[i,1] <- 1 - prod((1 - Pc)[sharedC == 1])
  }
}

# II - C as Matrix Implementation:

CM <- matrix(ncol = 0, nrow = 20)
for (i in Pc) {
  CM <- cbind(CM, rbinom(20,1, prob = i))
}
CM <- cbind(V, CM)

edges <- rep(0, length = ((length(V) * (length(V)-1)) / 2))

for (i in 1:dim(CM)[1]) {
  memb <- CM[i,2:4]
  if (i == max(V)) { break }
  for (j in (i+1):dim(CM)[1]) {
    index <- (i - 1) * (length(V) - (i/2)) + j - i # Triangular Matrix Approach
    logical <- CM[i,2:4] > 0 & CM[j,2:4] > 0
    if (sum(logical) == 0) {
      edges[index] <- e
    } else {
      edges[index] <- 1 - prod((1 - Pc)[logical])
    }
    # print(paste0(as.character(i),'-',as.character(j)))
  }
}

## From AGM to BigCLAM (relaxed; continuous membership values instead of binary)

# Generate the Community Membership Strength Matrix
CMR <- matrix(ncol = 0, nrow = 20)
for (i in Pc) { 
  vec <- round(runif(length(V), min = 0, max = 3), 2)
  vec[rbinom(length(V), 1, prob = 0.5) == 1] <- 0
  CMR <- cbind(CMR, vec)
}
CMR <- cbind(V, CMR)

edges2 <- rep(0, length = ((length(V) * (length(V)-1)) / 2))

# Generate Probability of Connections
for (i in 1:dim(CMR)[1]) {
  comStrI <- CMR[i,2:4]
  if (i == length(V)) { break }
  for (j in (i+1):dim(CMR)[1]) {
    index <- (i - 1) * (length(V) - (i/2)) + j - i # Triangular Matrix Approach
    comStrJ <- CMR[j,2:4]
    probConn <- 1 - exp(-comStrI %*% comStrJ) # PA(u,v) = 1 - exp(-FuA * FvA)
    if (probConn == 0) {
      edges2[index] <- e
    } else {
      edges2[index] <- probConn
    }
    # print(paste0(as.character(i),'-',as.character(j)))
  }
}

# BigCLAM: How to find F 
# Task: given a network G(V,E), estimate F. 
# Find F that maximizes the likelihood that F generated G. 

V <- 1:20
E <- matrix(rbinom(n = 400, size = 1, prob = 0.3), ncol = 20) + diag(20)
E[E > 1] <- 1

# Random Initialisation of F
set.seed(1)
MF <- matrix(round(runif(60, max = 1, min = 0), 2), ncol = 3) 
MF <- cbind(V, MF)

eta <- 0.001

# Applying Gradient Ascent to Find the Likeliest Matrix F that Generated E
repeat {
  MFO <- MF
  MFI <- MF
  MFColSums <- colSums(MF[,-1])
  for (i in 1:dim(MF)[1]) {
    a <- MF[i, 2:4]
    logical <- E[i,]
    logical[i] <- 0
    neighborNodes <- matrix(MF[(logical == 1), 2:4], ncol = 3)
    neighborGradient <- vector(length = 3)
    for (j in 1:dim(neighborNodes)[1]) {
      b <- neighborNodes[j,]
      gradient <- b * (exp(-a %*% b) / (1 - exp(-a %*% b)))
      neighborGradient <- neighborGradient + gradient
    }
    nonNeighborGradient <- MFColSums - a - colSums(neighborNodes)
    finalGradient <- neighborGradient - nonNeighborGradient
    vec <- MF[i, 2:4] + (eta * finalGradient)
    vec[vec < 0] <- 0.0
    MFI[i, 2:4] <- vec
  }
  MF <- MFI
  if (sum(MFO[,2:4] - MF[,2:4]) < 0.01) { 
    writeLines(paste0("Convergence! \n Iterations needed: ", as.character(x), "\n F Delta's of the final two iterations: "))
    print(cbind(V,round(MF[,2:4] - MFO[,2:4], 3)))
    break 
  }
}

edges3 <- rep(0, length = ((length(V) * (length(V)-1)) / 2))

# Generate Probability of Connections
for (i in 1:dim(MF)[1]) {
  comStrI <- MF[i,2:4]
  if (i == length(V)) { break }
  for (j in (i+1):dim(MF)[1]) {
    index <- (i - 1) * (length(V) - (i/2)) + j - i # Triangular Matrix Approach
    comStrJ <- MF[j,2:4]
    probConn <- 1 - exp(-comStrI %*% comStrJ) # PA(u,v) = 1 - exp(-FuA * FvA)
    if (probConn == 0) {
      edges3[index] <- e
    } else {
      edges3[index] <- probConn
    }
    # print(paste0(as.character(i),'-',as.character(j)))
  }
}

# Transforming Triangular Array into Node Matrix
nodeMatrixResult <- matrix(rep(0, 400), ncol = 20, nrow = 20)
for (row in 1:dim(nodeMatrixResult)[1]) {
  for (column in 1:dim(nodeMatrixResult)[2]) {
    if (row == column) {
      nodeMatrixResult[row, column] <- 1.000
    } else {
      vec <- sort(c(row, column))
      i <- vec[1]
      j <- vec[2]
      index <- (i - 1) * (length(V) - (i/2)) + j - i
      value <- round(edges3[index], 3)
      nodeMatrixResult[row, column] <- value
    }
  }
}

# Compare: 
E[1:10,1:10]
nodeMatrixResult[1:10,1:10]
