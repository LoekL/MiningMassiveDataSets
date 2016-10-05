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
