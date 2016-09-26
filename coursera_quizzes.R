# Mining Massive Datasets -- Coursera Quizzes

## Quiz 1

### Q1

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

