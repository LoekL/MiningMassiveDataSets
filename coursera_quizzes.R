# Mining Massive Datasets -- Coursera Quizzes

### Week 1 - Quiz 1

## Q1
# Suppose we compute PageRank with a β of 0.7, and we introduce the additional constraint
# that the sum of the PageRanks of the three pages must be 3, to handle the problem that
# otherwise any multiple of a solution will also be a solution. Compute the PageRanks a, b,
# and c of the three pages A, B, and C, respectively. Then, identify from the list below, the
# true statement.

# probablity of following the link
b <- 0.7

# number of nodes 
n = 3

# matrix of vectors 
M <- matrix(c(0, 0.5, 0.5, 0, 0, 1, 0, 0, 1), nrow = 3)

# rank vector (i.e., initialize page rank) 
r <- matrix(c(1/3, 1/3, 1/3), ncol = 1)

# jumping to a random page/teleportation (i.e., 1 - B)
t <- matrix(nrow = n, ncol = 1, rep(1/n,n))

# PageRank: Power Iteration (i.e., iterate until r does not change)
r <- b * M %*% r + (1-b)*t
r

(r[1] + r[2]) * 3
(r[1] + r[2]) * 3
(r[2] + r[3]) * 3
(r[1] + r[3]) * 3 # Correct

## Q2
# Suppose we compute PageRank with β=0.85. Write the equations for the PageRanks a, b,
# and c of the three pages A, B, and C, respectively. Then, identify in the list below, one of the
# equations.

beta <- 0.85
n <- 3
M <- matrix(c(0, 0.5, 0.5, 0, 0, 1, 1, 0, 0), nrow = 3)
r <- matrix(nrow = n, ncol = 1, rep(1/n, n))
t <- matrix(nrow = n, ncol = 1, rep(1/n, n))

for (i in 1:10) {
  r <- b * M %*% r 
  s <- sum(r) # see what potentially "leaked out"
  r <- r + (1-s) * t # feed this back to rank-vector evenly (t = 1/3 * 3)
}

r[2] == (0.575 * r[1]) + (0.15 * r[3])
r[2] == (0.475 * r[1]) + (0.05 * r[3])
r[1] == r[3] + (0.15*r[2])
0.95*r[2] == (0.475*r[1]) + (0.05*r[3])

# -- No answer found..

## Q3
# Assuming no "taxation," compute the PageRanks a, b, and c of the three pages A, B, and # --> Taxation is everything around β, the evenly back distributing of leaked rank across nodes
# C, using iteration, starting with the "0th" iteration where all three pages have rank a = b = c
# = 1. Compute as far as the 5th iteration, and also determine

n <- 3
M <- matrix(c(0, 0.5, 0.5, 0, 0, 1, 1, 0, 0), nrow = 3)
r <- matrix(nrow = n, ncol = 1, rep(1, n))

for (i in 1:5) {
  r <- M %*% r 
}

a <- r[1]
b <- r[2]
c <- r[3]

# After 5th iteration:
b == 5/8 # True

## Q4
# Suppose our input data to a map-reduce operation consists of integer values (the keys are not important). 
# The map function takes an integer i and produces the list of pairs (p,i) such that p is a prime divisor of i. For example, map(12) = [(2,12), (3,12)].
# The reduce function is addition. That is, reduce(p, [i , i , ...,i ]) is (p,i +i +...+i ).
# Compute the output, if the input is the set of integers 15, 21, 24, 30, 49. Then, identify, in the list below, one of the pairs in the output.

t1 <- list(2, c(24, 30))
t2 <- list(3, c(15, 21, 24, 30))

sum(t1[[2]])

# (2, 54)

### Week 2 - Quiz 1

## Q1
# The edit distance is the minimum number of character insertions and character deletions required
# to turn one string into another. Compute the edit distance between each pair of the strings he,
# she, his, and hers. Then, identify which of the following is a true statement about the number of
# pairs at a certain edit distance.

words <- c('he', 'she', 'his', 'hers')
pairs <- t(combn(words, 2))
pairs <- as.data.frame(pairs)

for (i in 1:dim(pairs)[1]) {
  pairs$dist[i] <- adist(pairs$V1[i], pairs$V2[i])[1,1]
}

# A: There are 2 pairs at distance 3.

## Q2
# Consider the following Matrix:

R1 <- matrix(c(0,1,1,0), ncol = 4)
R2 <- matrix(c(1,0,1,1), ncol = 4)
R3 <- matrix(c(0,1,0,1), ncol = 4)
R4 <- matrix(c(0,0,1,0), ncol = 4)
R5 <- matrix(c(1,0,1,0), ncol = 4)
R6 <- matrix(c(0,1,0,0), ncol = 4)

Matrix <- rbind(R1, R2, R3, R4, R5, R6)

# Perform a minhashing of the data, with the order of rows: R4, R6, R1, R3, R5, R2. 
# Which of the following is the correct minhash value of the stated column?
# Note: we give the minhash value in terms of the original name of the row, rather than the order of the row in the
# permutation. These two schemes are equivalent, since we only care whether hash values
# for two columns are equal, not what their actual values are.

MatrixPerm <- rbind(R4, R6, R1, R3, R5, R2)

# A: The minhash value for C4 is R3
# The minhash value of MatrixPerm[,4] is in the 4th row, which corresponds to R3.

## Q3
# Here is a matrix representing the signatures of seven columns, C1 through C7.

Matrix <- matrix(c(
  1,2,1,1,2,5,4,
  2,3,4,2,3,2,2,
  3,1,2,3,1,3,2,
  4,1,3,1,2,4,4,
  5,2,5,1,1,5,1,
  6,1,6,4,1,1,4
  ),
  ncol = 7,
  byrow = TRUE)

# Suppose we use locality-sensitive hashing with three bands of two rows each. Assume
# there are enough buckets available that the hash function for each band can be the identity
# function (i.e., columns hash to the same bucket if and only if they are identical in the band).
# Find all the candidate pairs, and then identify one of them in the list below.

# Hash Function: Cantor's Pairing Function: (a + b) * (a + b + 1) / 2 + a

bucket1 <- vector(length = dim(Matrix)[2])
bucket2 <- vector(length = dim(Matrix)[2])
bucket3 <- vector(length = dim(Matrix)[2])

for (i in 1:dim(Matrix)[2]) {
  group1 <- (Matrix[1,i] + Matrix[2,i]) * (Matrix[1,i] + Matrix[2,i] + 1) / 2 + Matrix[1,i]
  bucket1[i] <- group1
  group2 <- (Matrix[3,i] + Matrix[4,i]) * (Matrix[3,i] + Matrix[4,i] + 1) / 2 + Matrix[3,i]
  bucket2[i] <- group2
  group3 <- (Matrix[5,i] + Matrix[6,i]) * (Matrix[5,i] + Matrix[6,i] + 1) / 2 + Matrix[5,i]
  bucket3[i] <- group3
}

candidateGenerator <- function(bucket) {
  index <- 1:7
  count <- 0
  for (i in bucket[duplicated(bucket)]) {
    sub <- i == bucket
    if (count == 0) {
      matrix <- as.matrix(index[sub], ncol = 2)
      count <- count + 1
    } else {
      result <- as.matrix(index[sub], ncol = 2)
      matrix <- cbind(matrix, result)
    }
  }
  return(matrix)
}

candidatesBucket1 <- candidateGenerator(bucket1) # 1 & 4 | 2 & 5
candidatesBucket2 <- candidateGenerator(bucket2) # 1 & 6
candidatesBucket3 <- candidateGenerator(bucket3) # 1 & 3 | 4 & 7

# A: C2 and C5

## Q4
# Find the set of 2-shingles for the "document":

doc1 <- 'ABRACADABRA'

# and also for the "document":

doc2 <- 'BRICABRAC'

# Answer the following questions:
# 1. How many 2-shingles does ABRACADABRA have?

shingles1 <- c()

for (i in 1:(nchar(doc1) - 1)) {
  shingles1[i] <- substring(doc1, i, (i+1))
}

length(unique(shingles1)) # 7

# 2. How many 2-shingles does BRICABRAC have?

shingles2 <- c()

for (i in 1:(nchar(doc2) - 1)) {
  shingles2[i] <- substring(doc2, i, (i+1))
}

length(unique(shingles2)) # 7

# 3. How many 2-shingles do they have in common?

uniqueShingles1 <- unique(shingles1)
uniqueShingles2 <- unique(shingles2)

sum(uniqueShingles1 %in% uniqueShingles2) # 5
# or:
length(intersect(uniqueShingles1, uniqueShingles2)) # 5

# 4. What is the Jaccard similarity between the two documents"?

jaccardSim <- length(intersect(uniqueShingles1, uniqueShingles2)) / 
  length(union(uniqueShingles1, uniqueShingles2)) # 0.556

# A: The Jaccard similarity is 5/9.

##Q5
# Suppose we want to assign points to whichever of the points (0,0) or (100,40) is nearer.
# Depending on whether we use the L or L norm, a point (x,y) could be clustered with a different
# one of these two points. For this problem, you should work out the conditions under which a point
# will be assigned to (0,0) when the L norm is used, but assigned to (100,40) when the L norm is
# used. Identify one of those points from the list below.

LNorm <- function(power, a, b) {
  distances <- abs(a - b)
  sum <- 0
  for (i in 1:length(distances)) {
    sum <- sum + distances[i]^power
  }
  result <- sum^(1/power)
  return(result)
}

point1 <- c(0, 0)
point2 <- c(100, 40)

option1 <- c(61, 8)
option2 <- c(53, 18)
option3 <- c(58, 13)
option4 <- c(50, 18)

# option1 closer to point1 using L1
LNorm(1, point1, option1)
LNorm(1, point2, option1)

# Whereas option1 closer to point2 using L2
LNorm(2, point1, option1)
LNorm(2, point2, option1)

# A: (61, 8)

### Week 2 - Quiz 2

## Q1
# Suppose we have transactions that satisfy the following assumptions:
# - s, the support threshold, is 10,000.
# - There are one million items, which are represented by the integers 0,1,...,999999.
# - There are N frequent items, that is, items that occur 10,000 times or more.
# - There are one million pairs that occur 10,000 times or more.
# - There are 2M pairs that occur exactly once. M of these pairs consist of two frequent items, the other M each have at least one nonfrequent item.
# - No other pairs occur at all.
# - Integers are always represented by 4 bytes.
#
# Suppose we run the a-priori algorithm to find frequent pairs and can choose on the second pass between the triangular-matrix method 
# for counting candidate pairs (a triangular array count[i][j] that holds an integer count for each pair of items (i, j) where i < j) 
# and a hash table of item-item-count triples. Neglect in the first case the space needed to translate between original item numbers 
# and numbers for the frequent items, and in the second case neglect the space needed for the hash table. Assume that item numbers and counts are always 4-byte integers.
# As a function of N and M, what is the minimum number of bytes of main memory needed to execute the a-priori algorithm on this data? 
# Demonstrate that you have the correct formula by selecting, from the choices below, the triple consisting of values for N, M, and the 
# (approximate, i.e., to within 10%) minumum number of bytes of main memory, S, needed for the a-priori algorithm to execute with this data.

s <- 10000
items <- 0:999999

# Triangular Matrix Approach

N1 <- 10000
N2 <- 20000
N3 <- 30000

triangularBytes <- function(N) {
  possiblePairs <- (N * (N - 1)) / 2
  bytes <- possiblePairs * 4
  return(bytes)
}

R1 <- triangularBytes(N1) # 199980000 ~ 199M bytes
R2 <- triangularBytes(N2) # 799960000 ~ 799M bytes
R3 <- triangularBytes(N3) # 1799940000 ~ 1.8B bytes

# Triples Approach

M1 <- 50000000
M2 <- 60000000
M3 <- 100000000
M4 <- 200000000

# 1 milion pairs that are frequent
# 2M pairs that occur exactly once, 1M of these pairs have 2 freq items. other M at least one nonfreq.

tripleBytes <- function(M) {
  frequentPairs <- 1000000 # pairs whose items are both frequent and are frequent itself
  nonFrequentPairs <- M # pairs whose items are both frequent, yet are not frequent itself
  totalCounts <- frequentPairs + nonFrequentPairs
  bytes <- totalCounts * 3 * 4 # we need to store 3 integers of each 4 bytes
  return(bytes)
}

T1 <- tripleBytes(M1) # ~ 612M bytes
T2 <- tripleBytes(M2) # ~ 732M bytes
T3 <- tripleBytes(M3) # ~ 1.2B bytes
T4 <- tripleBytes(M4) # ~ 2.4B bytes

# A: N = 30,000; M = 200,000,000; S = 1,800,000,000

## Q2
# Imagine there are 100 baskets, numbered 1,2,...,100, and 100 items, similarly numbered. 
# Item i is in basket j if and only if i divides j evenly. 
# For example, basket 24 is the set of items {1,2,3,4,6,8,12,24}. 
# Describe all the association rules that have 100% confidence. Which of the following rules has 100% confidence?

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

# 100% Confidence is achieved when the superset is as frequent as its subset.
# Example: {1, 2, 3} --> 4
# - If there is a 100% confidence, {1, 2, 3} goes together with {4}.
# - Hence {1, 2, 3} will be exactly as frequent as {1, 2, 3, 4}.

s1 <- c(4, 10, 12)
s2 <- c(2, 4)
s3 <- c(1)
s4 <- c(1, 4, 7)

a1 <- c(80)
a2 <- c(8)
a3 <- c(2)
a4 <- c(14)

frequency <- function(x, baskets) {
  size <- length(x)
  count <- 0
  for (i in 1:length(baskets)) {
    if (sum(x %in% baskets[[i]]) == size) {
      count <- count + 1
    }
  }
  return(count)
}

# Option 1
freqS1 <- frequency(s1, baskets) # 1
freqS1A1 <- frequency(c(s1, a1), baskets) # 0
freqS1 == freqS1A1 # FALSE

# Option 2
freqS2 <- frequency(s2, baskets) # 25
freqS2A2 <- frequency(c(s2, a2), baskets) # 12
freqS2 == freqS2A2 # FALSE

# Option 3
freqS3 <- frequency(s3, baskets) # 100
freqS3A3 <- frequency(c(s3, a3), baskets) # 50
freqS3 == freqS3A3 # FALSE

# Option 4
freqS4 <- frequency(s4, baskets) # 3
freqS4A4 <- frequency(c(s4, a4), baskets) # 12
freqS4 == freqS4A4 # TRUE

# A: Option 4

## Q3
# Suppose ABC is a frequent itemset and BCDE is NOT a frequent itemset. 
# Given this information, we can be sure that certain other itemsets are frequent 
# and sure that certain itemsets are NOT frequent. Other itemsets may be either frequent 
# or not. Which of the following is a correct classification of an itemset?

# A:
# 1 - CDE is not frequent.
#     + FALSE: Even though BCDE is not frequent, that does not mean a subset cannot be.
# 2 - BCD can be either frequent or not frequent. 
#     + TRUE: Even though BCDE is not frequent, it could be in the negative border, meaning BCD could still be frequent.
#             It could also simply not be frequent.
# 3 - ABCDEF can be either frequent or not frequent.
#     + FALSE: BCDE is a subset of ABCDEF and is NOT frequent. Therefore ABCDEF cannot be frequent.  
# 4 - AB can be either frequent or not frequent.
#     + FALSE: Since ABC is frequent, all its subsets must also be frequent.

### Week 2 - Quiz 3

## Q1
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
# Demonstrate that you have the correct formula by indicating which of the following is a value for 
# S and a value for P that is approximately (i.e., to within 10%) the largest possible value of P for that S.

# Pass 1

S <- S4
P <- P4

delta <- function(s, p) {
  
  numFreqItems <- 250000
  numFreqPairs <- 1000000
  numBuckets <- (s - (numFreqItems * 4 * 2)) / 4 # 4-bytes per bucket | largest possible number of buckets we can maintain
  
  # Pass 1
  SP1 <- (numFreqItems * (4 * 2)) + # Item count
    (numBuckets * 4) # Hash table with max number of buckets
  
  # Pass 2
  SP2 <- (numFreqItems * 4) + # Frequent items (no longer counts, hence not * 2)
    (numBuckets / 8) + # bitmap
    (((numFreqPairs * p / numBuckets) + 1) * (4 * 3))
  
  print(paste0("Memory Usage Pass 1: ", as.character(round(SP1 / 1000000, digits = 0)), "M Bytes"))
  print(paste0("Memory Usage Pass 2: ", as.character(round(SP2 / 1000000, digits = 0)), "M Bytes"))
  
  delta1 <- s - SP1 # This one should always be near 0
  delta2 <- s - SP2 # This one can vary
  
  # Large negative numbers mean OoM!
  print(paste0("Memory Left during Pass 1: ", as.character(round(delta1 / 1000000, digits = 0)), "M Bytes"))
  print(paste0("Memory Left during Pass 2: ", as.character(round(delta2 / 1000000, digits = 0)), "M Bytes"))
  
}

# A:
option1 <- delta(100000000, 120000000)
option2 <- delta(300000000, 3500000000)
option3 <- delta(1000000000, 10000000000)
option4 <- delta(200000000, 800000000) # Closest

## Q2
# During a run of Toivonen's Algorithm with set of items {A,B,C,D,E,F,G,H} a sample is found to 
# have the following maximal frequent itemsets: {A,B}, {A,C}, {A,D}, {B,C}, {E}, {F}. 
# Compute the negative border. Then, identify in the list below the set that is NOT in the negative border.

# {H} - {H} is not frequent, yet it's subset {} (empty set) is, hence {H} is in the negative border | FALSE
# {A,B,D} - {A,B,D} is not frequent, two of its subsets are frequent: {A,B} & {A,D} but not {B,D} hence it is NOT in the negative border | TRUE 
# {B,E} - {B,E} is not frequent, yet {E} is frequent as is {B} ({B,E} is frequent, monotonicity dicatates that singleton {B} must then be frequent as well)
#         therefore {B,E} must be in the negative border | FALSE
# {D,F} - Similar as above, singleton {F} is frequent, {D} as well because of {A,D} being frequent, hence {D,F} is in the negative border | FALSE

# A: {A,B,D} is not in the negative border.
