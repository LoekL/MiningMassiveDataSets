== Mining Massive Datasets ==

-- Week 1 --

## 2 - 1 - Distributed File Systems (15-50).mp4

- Distributed File System
- Computational Model
- Scheduling and Data FLow
- Refinements

[Single Node Architecture]

- CPU + Memory: Machine Learning, Statistics # All is stored in memory, no disk required
- CPU + Memory + Disk: 'Classical' Data Mining # However sometimes this is still not enough!

Motivation: Google Example
- 10 billion web pages
- Average size of webpage = 20kb
- 10 billion * 20kb = 200TB
- Disk read bandwidth = 50 MB/sec
- Time to read = 4 milion seconds = 46+ days # ...
- Even longer to do something useful with the data # Hence we need an alternative!

The main idea is to cut the work into several pieces and run them in parallel

[Cluster Architecture]

- Each rack contains 16-64 commodity linux nodes
- 1 Gbps between any pair of nodes in a rack
- 2-10 Gbps backbone between racks

# In 2011 it was guestimated that Google had 1M machines - to get a sense of the scale.

Cluster Computing Challenges
1 - Node failure
    + A single server can stay up for 3 years (1000 days)
    + 1000 servers in cluster => 1 failure/day
    + 1M servers in cluster => 1000 failures/day
  - How to store data persistently and keep it available if nodes can fail? # persistently means you can always retrieve it again
  - How to deal with node failures during a long-running computation?
    + That node had data pertaining to the calculation.
    + Worst case you need to start over, but then another node may fail.
      * You need an infrastructure that can deal with these kind of failures.
2 - Network bottleneck
    + Network bandwidth = 1 Gbps # Gigabit per second
    + Moving 10 TB (terabytes) takes approximately 1 day
3 - Distributed programming is hard!
    + Need a simple model that hides most of the complexity.

[Map-Reduce]

- Map-Reduce addresses the challenges of cluster computing
  1 - Store data redundantly on multiple nodes for persistence and availability.
  2 - Move computation close to data to minimize data movement.
  3 - Simple programmatic model to hide the complexity of all this magic.

[Redundant Storage Inrastructure]

1 - Distributed File System # provides a Redundant Storage Infrastructure
    + Provides global file namespace, redundancy, and availability
    + E.g., Google GFS, Hadoop HDFS
    + Typical usage pattern
      * Huge files (100s of GB to TB)
      * Data is rarely updated in place
      * Reads and appends are common
    + Data kept in 'chunks' spread across machines
    + Each chunk replicated on different machines
      * Ensures persistence and availability

Chunk server 1 -- Chunk Server 2 --- Chunk Server 3 ... Chunk Server n
c0 c1 c5 c2       c5 c1 c3 d0        d1 c2 c4 d0         d1 c0 c4 c3

- Each chunk is saved twice, but never on the same machine (Chunk Server)!
- Chunk servers also serve as compute servers, you move the computation to where the data is # and avoid moving data around

[Components of a Distributed File System]

1 - Chunk Servers
    + File is split into contiguous chunks (16-64MB) # next or together in sequence
    + Each chunk is replicated (usually 2x or 3x) # 3x most common
    + Try to keep replicas in different racks # not only on different servers, also different racks
      * In case switch fails
2 - Master node # called as such in GFS / Google File System
    + A.k.a. Name Node in Hadoops HDFS
    + Stores metadata about where files are stored
    + Might be replicated itself
3 - Client library for file access
    + Talks to master to find chunk servers
    + Connects directly to chunk servers to access data

## 2 - 2 - The MapReduce Computational Model (22-04)

[Programming Model: MapReduce]

Warm-up task:
- We have a huge text document.
- Count the number of times each distinct word appears in the file.
- Same applications
  + Analyze web server logs to find popular URLs
  + Term statistics for search

Task: Word Count

- Case 1
  + File too large for memory, but all <word, count> pairs fit in memory.
    * Create a hash table with the words and respective counts.
    * Every time you come across the word you either initialise it to 1 or add 1. # .get method in Python
- Case 2
  + Even the <word, count> pairs do not fit in memory.
    * Use a unix command line way of doing this: "words(doc.txt) | sort | uniq -c"
    * Where 'words' takes a file and outputs the words in it, one per line.
    * 'sort' sorts the words, 'unique -c' runs through the words and collects a count per word
    * The output will be word/count pairs
  + Case 2 captures the essence of MapReduce
    * Great thing is that it is naturally parallelizable

[MapReduce: Overview]

1 - Map # words(doc.txt)
    + Scan input file record-at-a-time
    + Extract something you care about from each records (keys)
2 - Group by key # sort
    + Sort and Shuffle
3 - Reduce # uniq -c
    + Aggregate, summarize, filter or transform
    + Write the results

Outline stays the same, Map and Reduce change to fit the problem.

MapReduce
1 - the Map Step
    + Input key-value pairs --> Intermediate key-value pairs
2 - the Reduce Step
    + Intermediate key-value pairs --> group by key --> key-value groups
      # for example in several chunks of text the word 'apple' was found,
      # in varying counts, we need to group these together.
    + key-value groups --> reduce --> Output key-value-paris
      # Here you sum (or do something more complicated) up the values to get to a single value output key-value pair.

# Full cycle:
Input key-value pairs --> Intermediate key-value pairs --> group by key --> key-value groups --> reduce --> Output key-value-paris

More formally..

- Input: a set of key-value pairs
- Programmer specifies two methods:
  + Map(k, v) --> <k1, v1>*
    * Takes a key-value pair and outputs a set of (intermediate) key-value pairs
    * There is one Map call for every (k,v) pair
  + Reduce(k1, <v1>*) --> <k1, v11>* # all v1 values are added to the same key in some way (e.g. sum)
    * All values v1 with same key k1 are reduced together
    * There is one Reduce function call per unique key k

-- MapReduce: Word Counting

1 - MAP: Read input and produces a set of key-value pairs
# Provided by the programmer.

Big document --> (key,value)

The crew of the space                   (The, 1)
shuttle Endeavor recently               (crew, 1)
returned to Earth as           -->      (of, 1)
ambassadors, harbringers of             (the, 1)
a new era of...                         (space, 1)
                                          ...

Group by key: collect all pairs with the same key

(crew, 1)
(crew, 1)
(space, 1)
(the, 1)
(the, 1)
(the, 1)
  ...

2 - Reduce: Collect all values belonging to the key and output
# Provided by the programmer.

(crew, 2)
(space, 1)
(the, 3)
(shuttle, 1)
  ...

- Within Map-Reduce this all does not run on a single node.
- The source data (document) is cut up into several chunks, which reside on differing nodes.
- All MAP tasks are run in parralel on each node.
- The MAP output of all nodes is then copied onto a single node.
  + Once this is done the data is sorted, and can then progress into the final reduce stop, done by a single node.
- You can opt to use more reduce nodes.
  + The system will be smart enough that for instances all instances of 'the' will end up
    in the same reduce node.
  + The system does this via a hash function, which maps a specific key to a specific reduce node.
    * It hashes each map,key and determines a single reduce node to ship that tuple to.
    * This ensures that all tuples with the same key end up at the same reduce node.
    * Then they can get sorted, and reduced in parallel, which is fine because keys are grouped together.

The operative attribute however is that the device can access any required record immediately on demand.
The opposite is sequential access, where a remote element takes longer time to access.

- Random Access # https://en.wikipedia.org/wiki/Random_access
- Sequential Access # Map-Reduce - https://en.wikipedia.org/wiki/Sequential_access

A typical illustration of this distinction is to compare an ancient scroll (sequential;
all material prior to the data needed must be unrolled) and the book (direct: can be immediately flipped open to any arbitrary page).
A more modern example is a cassette tape (sequential — one must fast forward through earlier songs to get to later ones) and a CD
(direct access — one can skip to the track wanted, knowing that it would be the one retrieved).

Sequential reads are much faster, which is why Map-Reduce is built around these kind of reads only.

# JavaScript

map(key, value):
// key: document name; value: text of the document
  for each word w in value:
    emit(w,1) # produce a tuple (w, 1)

reduce(key, values):
// key: a word; value; an iterator over counts
    results = 0
    for each count v in values:
        result += v
    emit(key, result)

-- Example: Host Size

- Suppose we have a large web corpus with a metadata file formatted as follows:
  + Each record of the form: (URL, size, date, ...)
- For each host, find the total number of bytes # there can be many URLs with the same host name

1 - Map
    + For each record, output (hostname(URL), size)
      * Look at the URL of the record, output the host name
2 - Reduce
    + Sum the sizes for each host

-- Example: Language Model

- Count number of times each 5-word sequence occurs in a large corpus of documents.

1 - Map
    + Extract (5-word sequence, count) from document
2 - Reduce
    + Combine the counts (add them up)

## 2 - 3 - Scheduling and Data Flow (12-43).mp4

Map-Reduce environment takes care of:
- Partitioning the input data
- Scheduling the programs execution across a set of machines
- Performing the group by key step
- Handling node failures
- Managing required inter-machine communication

[Data Flow]

- Input and final output are stored on the distributed file system (DFS):
  + Scheduler tries to schedule map tasks 'close' to physical storage location
    of input data.
- Intermediate results are stored on local FS (LFS) of Map and Reduce workers.
- Output is often input to another Map-Reduce task

[Coordination: Master]

- Master node takes care of coordination:
  + Task status: (idle, in-progress, completed)
  + Idle tasks get scheduled as workers become available
  + When a map task completes, it sends the master the location and sizes of its
    *R* intermediate files, one for each reducer (R = number of required reducers).
  + Master pushes this info to reducers
- Master pings workers periodically to detect failures.

[Dealing with Failures]

- Map worker failure
  + Map tasks completed or in-progress at workers are reset to idle on that failed node.
  + Idle tasks eventually rescheduled on other worker(s).
- Reduce worker failure
  + Only in-progress tasks are reset to idle.
    * The output of the reducer is the final output, which is written to the DFS and not LFS
      meaning the results are not lost!
  + Idle Reduce tasks restarted on other worker(s).
- Master failure
  + MapReduce task is aborted and client is notified.
    * The one scenario where the task will have to be restarted from scratch.
    * User can reboot node.

[How many Map and Reduce jobs?]

- M map tasks, R reduce tasks
- Rule of thumb:
  + Make M much larger than the number of nodes in the cluster.
  + One DFS chunk per map is common.
  + Improves dynamic load balancing and speeds up recover from worker failures.
- Usually R is smaller than M
  + Because output is spread across R files

## 2 - 4 - Combiners and Partition Functions (12-17) [Advanced].mp4

[Refinement: Combiners]

- Often a Map task will produce many pairs of the form (k,v1), (k,v2), ... for the same key k.
  + E.g. popular words in the word count example
- We can save network time by pre-aggregating values in the mapper:
  + combine(k, list(v1)) --> v2 # output of the combiner function is key/value pair where value is scalar. Input is a key + list of values.
  + Combiner is usually the same as the reduce function

- Back to our word counting example:
  + Combiner combines the values of all keys of a single mapper (single node)
  + Much less data needs to be copied and shuffled!

(A,B)          (A,1)           (A,3)
(A,C)          (B,1)           (B,4)
(A,D)  mapper  (A,1)  combiner (C,2)
(B,E)   -->    (C,1)     -->   (D,1)
(B,D)          (A,1)           (E,1)
(C,B)          (D,1)

This output is then sent to reducer, where it is grouped by key & eventually again combined into single key/value tuples.

- Combiner trick works only if reduce function is commutative and associative.
  + Sum # both commutative and associative
    * Commutative --> a + b = b + a
    * Associative --> (a + b) + c = a + (b + c)
  + The idea is that you can break the data in two buckets, sum the values in each buckets & afterwards sum the sub-totals to get
    to the final result. This needs to be accurate in any form or shape.
  + Average # not commutative and associative -- if samples are not same size?
    * A trick would be to let the combiner output: (sum, count) --> the reducer will calculate by dividing sum(sums) with the sum(counts)
  + Median --> this will never work within a combiner, can only be done in reducer stage.

[Refinement: Partition Function]

- Want to control how keys get partitioned
  + The set of keys that go to a single reduce worker
- System uses a default partition function:
  + hash(key) % R # modulo remainder
    * This produces a number from 0 to R-1 --> to send the key to one of reducers R
- Sometimes useful to override the hash function:
  + E.g. hash(hostname(URL)) mod R ensures URLs from a host end up in the same output file

[Implementations]

- Google MapReduce
  + Uses Google File System (GFS) for stable storage
  + Not available outside Google
- Hadoop
  + Open-source implementation in Java
  + Uses HDFS for stable storage
  + Download: http://lucene.apache.org/hadoop/
- Hive, Pig
  + Provide SQL-like abstractions on top of Hadoop Map-Reduce layer

[Cloud Computing]

- Ability to rent computing by the hour
  + Additional services e.g., persistent storage
- E.g., Amazons 'Elastic Comput Cloud' (EC2)
  + S3 (stable storage)
  + Elastic Map Reduce (EMR)

## 2 - 5 - Link Analysis and PageRank (9-39).mp4

 Graph Data, 3 modules:
 1 - PageRank (+ SimRank) # Link analysis methods
 2 - Community Detection # Finding clusters of nodes in a network
 3 - Spam Detection # Identify nodes that are spam

Graph Data: Social Networks # Graph as a set of nodes with connections between them

- Web as a Graph
  + Nodes: Webpages
  + Edges: Hyperlinks
- Broad Question: How to organize the Web?
  + First try: Human curated
    * Web directories
      - Yahoo, DMOZ, LookSmart
  + Second try: Web Search
    + Information Retrieval investigates:
      * Find relevant docs in a small and trusted set
        - Newspaper articles, Patents, etc.
      * But: Web is huge, full of untrusted documents,
        random things, web spam, etc. # How do we retrieve only good, trustworthy webpages in this huge webgraph?

- Web Search - 2 Challenges:
  1 - Web contains many sources of information - Who to 'trust'?
      + Trick: trustworthy pages may point to each other!
  2 - What is the 'best' answer to query 'newspaper'?
      + No single right answer
      + Trick: pages that actually know about newspapers might
        all be pointing to many newspapers

We address these challenges by ranking nodes on the graph
- All web pages are not equally important
  + www.joe-schmoe.com vs. www.stanford.edu
  + There is large diversity in the web-graph node connectivity.
    Lets rank the pages by the link structure!

[Link Analysis Algorithm]

- We will cover the following Link Analysis approaches for computing importances
  of nodes in a graph:
  + Page Rank
  + Hubs and Authorities (HITS) # webpages can be called hubs or authorities
  + Topic-Specific (Personalized) Page Rank
  + Web Spam Detection Algorithms # how spammers try to manipulate the web graph to make websites seem important whereas they're not

#@ 2 - 6 - PageRank- The Flow Formulation (9-16).mp4

[Links as Votes]

- Idea: Links as votes
  + Page is more important if it has more links
    * In-coming links? Out-going links? # in-links are better, harder to fabricate
- Think of in-links as votes:
  + www.stanford.edu has 23,400 in-links
  + www.joe-schmoe.com has 1 in-link
- Are all in-links equal?
  + Links from important pages count more
  + Recursive question!

[Simple Recursive Formulation]

- Each links vote is proportional to the importance of its source page.
- If page j with importance rj has n out-links, each links gets rj/n votes.
- Page j s own importance is the sum of the votes on its in-links.

[PageRank: The 'Flow' Model]
- A 'vote' from an important page is worth more.
- A page is important if it is pointed to by other important pages.
- Define a 'rank' rj for page j.

rj =   Σ ri / di
    i -> j
di ... out-degree of node i

'Flow' equations:
ry = ry / 2 + ra / 2 # y had one self-link, one in-link from a, which had 2 outlinks
ra = ry / 2 + rm # rm only had one outlink
rm = ra / 2 # rm only had one in-link, from a which had 2 outlinks, etc.

[Solving the Flow Equations]

- 3 equations, 3 unknowns, no constants
  + No unique solution
  + All solutions equivalent modulo the scale factor
- Additional constraint forces uniqueness:
  + ry + ra + rm = 1
  + Solution: ry = 2/5, ra = 2/5, rm = 1/5
- Gaussian elimination method works for small examples, but we need a
  better method for large web-size graphs
  + https://en.wikipedia.org/wiki/Gaussian_elimination
- We need a new formulation!

## 2 - 7 - PageRank- The Matrix Formulation (8-02).mp4

- Stochastic adjacency matrix M
  + Let page i has di out-links
  + if i --> j, then Mji = 1/di else Mji =  0
  # if i points to j, then that value in matrix M (Mji) is non-zero,
  # being one / out-degree, or 1 / num_out_links_of_i
    * M is a column stochastic matrix
      - 'Column Stochastic' --> Means that the columns sum up to 1

- Rank vector r: vector with an entry per page
  + ri is the importance score of page i
  + Σi ri = 1

[Eigenvector Formulation]

- The flow equations can be written:

r = M * r

Note: x is an eigenvector with the corresponding eigenvalue λ if: A*x = λ*x # A = matrix, x = vector, λ = scalar
# https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors

- So the rank vector r is an eigenvector of the stochastic web matrix M
  + In fact, its first or principal eigenvector, with correpsonding eigenvalue 1
    - Largest eigenvalue of M is 1 since M is column stochastic
      * Why? We know r is a stochastic vector and each column of M sums to one, so Mr <= 1
    - [1 * r = M * r] == [A * x = λ * x]

- We can now efficiently solve for r!
  + The method is called Power Iteration
  + We can use it to compute the eigenvector of matrix M that corresponds to an eigenvalue of 1

## 2 - 8 - PageRank- Power Iteration (10-34).mp4

[Power Iteration Method]

- Given a web graph with N nodes, where the nodes are pages and edges are hyperlinks
- Power iteration: a simple iterative scheme
  + Suppose there are N web pages
  + Initialize: r^(0) = [1/N, ..., 1/N]
  + Iterate: r^(t+1) = M * r^(t)
  + Stop when |r^(t+1) - r^(t)| < epsilon --> using L1 norm # https://rorasa.wordpress.com/2012/05/13/l0-norm-l1-norm-l2-norm-l-infinity-norm/

[PageRank: How to Solve?]

- Power Iteration
  + Set rj = 1/N
  + 1: r''j = Σ ri / di
            i -> j
  + 2: r = r''
  + If not converged: goto 1

 [Random Walk Interpretation]

 - Imagine a random web surfer
   + At any time t, surfer is on some page i
   + At time t + 1, the surfer follows an out-link from i uniformly at random
   + Ends up on some page j linked from i
   + Process repeats indefinitely
- Let:
  + p(t) - vector whose ith coordinate is the probability that the surfer is at page i at time t
  + So, p(t) is a probability distribution over pages

- Where is the furfer at time t+1?
  + Follows a linke uniformly at random:
    p(t + 1) = M * p(t)
    # It's basically the sum of the probability of being on any of the nodes that point towards
    # the current node, where for each node you divide the probability by the number of out-links
    # since they had to choose the path that pointed to the current node
- Suppose the random walk reaches a steady  state
  p(t + 1) = M * p(t) = p(t)
  then p(t) is a stationary distribution of a random walk
- Our original rank vector r satisfies r = M * r
  + So, r is a stationary distribution for the random walk

[Existence and Uniqueness]

- A central result from the theory of random walks (a.k.a. Markov processes):

For graphs that satisfy certain conditions, the stationary distribution is unique,
and eventually will be reached no matter what the initial probability distribution is
at time t = 0.

# The power iteration will always converge to the same vector regardless of how we initialise it!

--> Importance r can be viewed as probability of being on a node <--

- After a long enough walk all nodes and connections have been mapped out and
  it will reach the stationary state, which is the eigenvector.

## 2 - 9 - PageRank- The Google Formulation (12-08).mp4

[PageRank: Three Questions]

rj^(t+1) = Σ (ri^(t)) / di # or equivalently r = M * r
         i -> j

1 - Does this converge?
2 - Does it converge to what we want?
3 - Are results reasonable?

[Does This Converge?]

A - The "Spider Trap" problem: a <--> b
- Example:

ra =  1  0  1  0
rb =  0  1  0  1

--> There is no convergence, ranks just constantly swap. No convergence!

B - The "Dead End" problem: a --> b

ra = 1  0  0  0
rb = 0  1  0  0

# a was important at first (manually set to 1), hence b gained in importance.
# However no inbound links to a makes a unimportant in next iteration, which in
# turn leads to b becoming unimportant as well

--> Again no convergence!

[PageRank: Problems]

- 2 problems:
  1 - Some pages are dead ends (have no out-links)
      + Such pages cause importance to 'leak out'
        * It receives its rank, but there is no way for the page to
          pass this importance on to other nodes
  2 - Spider traps: all out-links are within the group.
      + Eventually spider traps absorb all importance.
      + This part of the graph will get very high weight, whereas other parts of
        the graph very low weight.

## 2 - 10 - Why Teleports Solve the Problem (12-26).mp4

[Why Teleports Solve the Problem?]

r^(t+1) = M * r^(t)

Markov chains
- Set of states X
- Transition matrix P where Pij = P(Xt = i | Xt-1 = j) # Given I was in j step before, how likely will I transition to i
- π specifying the stationary probability of being at each state x ∈ X
- Goal is to find π such tat π = P *  π # r = M * r

[Why is This Analogy Useful?]

- Theory of Markov chains
- Fact: for any start vector, the power method applied to a Markov transition
  matrix P will converge to a unique positive stationary vector as long as P is
  stochastic, irreducible and aperiodic.

[Make M Stochastic]

- Stochastic: Every column sums up to 1.
- If a website has no out-links, it sums up to 1 and violates this requirement.
- Solution: Add green links - a fake link to all other nodes, including a node to itself (self-link)
  + By linking to all others you do not increase/decrease the importance of any specific node

A = M + a^(T) * (1/n * e)
# This process basically fills in fake links for all nodes with out-degree 0
# It fills them up with 1/n
# If out-degree is not-0, a defaults to 0 and nothing happens

- ai = 1 if node i has out degree 0, else it is 0
- e is a vector of all 1s

[Make M Aperiodic]

- A chain is periodic (!= aperiodic) if there exists k > 1 such that the interval between two
  visits to some state s is always a multiple of k.
  + The path is deterministic; alwasy the same, e.g.: a --> b --> c --> c --> ... # this is periodic, an infinite loop
- Solution: Add green links # we break this (periodicness) by introducing new paths (teleports)

[Make M Irreducible]

- From any state, there is a non-zero probability of going from any one state to
  any another # i.e. never get 'stuck' at a dead end, with no out-links
- Solution: add green links

[Solution: Random Jumps]

- Googles solution that does it all:
  + Makes M stochastic, aperiodic, irreducible
- At each step, random surfer has two options:
  + With probability β, follow a link at random
  + With probability 1-β, jump to some random page
- PageRank equation [Brin-Page, 98]

rj =  (Σ β * ri/di) + (1 - β) * 1/n
     i -> j

- rj = importance (r) of node j
- ri = probability of being on node i # or its importance
- d = out-degree

- The first part of the formula just sums up the importance of all the incoming links -- Σ β * ri/di
  + This happens with probability β --> the walker can decide to either follow a link, or teleport
- The random-walker could have also arrived through a teleport, with probability 1 - β -- (1 - β) * 1/n
  + Where n is the number of nodes

# The above formulation assumes that M has no dead ends. We can either preprocess
# M (bad!) or explicitly follow random teleport links with probability 1.0 from
# dead-ends.

[The Google Matrix]

- The Google Matrix A:

A = β * M + (1 - β) * 1/n * e * e^(T) # where e is a vector of all 1s

- A is stochastic, aperiodic and irreducible, so:

r^(t+1) = A * r^(t)

- What is β?
  + In practice β = 0.8, 0.9 (0.8 --> make 5 steps and jump)

## 2 - 11 - How we Really Compute PageRank (13-49).mp4

[Computing Page Rank]

- Key step is matrix-vector multiplication
  + r^(new) = A * r^(old)
- Easy if we have enough main memory to hold A, r^(old), r^(new)

A = β * M + ((1 - β) * (1/n)nxn) # where (1/n)nxn is an n*n matrix, populated with 1/n for all values # where n is number of nodes

- Even if M starts out with zeros (which requires minimal memory), after a few steps all
  values will be non-zero [because of the random teleportation stuff: ((1 - β) * (1/n)nxn)]
  + A is now what is called a dense-matrix
- Say N = 1 billion pages
  + We need 4 bytes for each entry (say)
  + 2 billion entries for vectors, approx 8GB
  + Matrix A has N^2 entries
    * 10^18 is a large number

Matrix Formulation
- Suppose there are N pages
- Consider page j, with dj out-links
- We have Mij = 1/dj when j --> i
     and Mij = 0 otherwise
- The random teleport is equivalent to:
  + Adding a teleport link from j to every other page
    and setting transition probability to (1 - β)/N # == (1 - β) * 1/N
  + Reducing the probability of following each out-link from 1/dj
    to β/dj # β * 1 / dj
  + Equivalent: Tax each page a fraction (1-β) of its score and
    redistribute evenly

[Rearranging the Equation]

r = A * r, where Aij = β * Mij + ((1-β)/N)

     N
ri = Σ Aij * rj
    j=1
==
     N
ri = Σ [β * Mij + ((1-β)/N)] * rj
    j=1
==
     N                  N
ri = Σ [β * Mij * rj] + Σ [(1-β)/N * rj]
    j=1                j=1
==
     N
ri = Σ β * Mij * rj + (1-β)/N # since Σ rj = 1
    j=1

So we get: r = β M * r + [1-β/N]N # [X]N = a vector of length N with all entries X

[Spare Matrix Formulation]

- We just rearranged the PageRank equation

r = (βM * r) + [(1-β)/N]N # βM * r returns a vector

# Where [(1-β)/N]N is a vector with all N entries (1-β)/N

- M is a sparse matrix (with no dead-ends)!
  + 10 links per node, approx 10N entries per website
  + therefore the 1bn * 1bn will contain almost all 0 values, making it non-dense --> sparse!
    * now we need much less memory to be able to work with it
- So in each iteration, we need to:
  + Compute r^(new) = βM * r^(old)
  + Add a constant value (1-β)/N to each entry in r^(new)
    * Note if M contains dead-ends then Σi ri^new < 1 and we also
      have to renormalize r^(new) so that it sums to 1

The symbol "∀" means "for all" or "for every"
"∀i∈{0,1,2}" means "for all i in the set {0,1,2}"

[PageRank: The Complete Algorithm]

- Input: Graph G and parameter β # "teleportation parameter beta"
  + Directed graph G with spider traps and dead ends # the algorithm is robust against spider traps and dead ends
  + Parameter β
- Output: PageRank vector r
  - Set rj^(0) = 1/N, t = 1
  - Do:
    * ∀j: rj`^(t) =  Σ β * (ri^(t-1)/di) # rj`^(t) == intermediate value of rj^(t), still needs more work to get to final value
                   i -> j

          rj`^(t) = 0 if in-degree of j is 0

    * Now re-insert the leaked PageRank: # if we have dead-ends the page-rank will leak out
    # so we determine how much has leaked out, and re-insert this into our page-ranks via S:
      ∀j: rj^(t) = rj`^(t) + (1-S)/N # where S = Σj rj`^(t) -- rj`^(t) == page-ranks BEFORE adding leaked page-rank
      # --> 1 - S is what is missing --> we divvy that equally and give back to all nodes
      # Now everything will sum up to one again

    * t = t + 1
  - While:

    Σj |rj^(t) - rj^(t-1)| > ε # --> using L1 norm -- https://rorasa.wordpress.com/2012/05/13/l0-norm-l1-norm-l2-norm-l-infinity-norm/

-- Week 2 --

## 3 - 1 - Finding Similar Sets (13-37).mp4

[Finding Similar Sets]

- Applications
- Shingling
- Minhashing
- Locality-Sensitive Hashing

[Applications of Set-Similarity]

- Many data-mining problems can be espressed as finding 'similar' sets:
  1 - Pages with similar words, e.g., for classification by topic.
  2 - NetFlix users with similar tastes in movies, for recommendation systems.
  # Think of a matrix where users are rows, columns movies and their ratings: 4/5 is high, blank if not given
  # We can now check for similarities between users and recommend other movies to similar users that they haven't seen yet.
  3 - Dual: movies with similar sets of fans, probably belong to a certain genre.
  4 - Entity resolution. # Determining the set of records that belong to the same individual.
  # That's why you're often asked for your phone number - but you can also input land/mobile number inconsistently.

[Similar Documents]

- Given a body of documents, e.g., the Web, find pairs of documents with a lot of text in common such as:
  + Mirror sites, or approximate mirrors.
    * Application: do not want to show both in a search.
- Plagiarism, including large quotations.
- Similar news articles at many news sites.
  + Application: cluster articles "by story".

[Three Essential Techniques for Similar Documents]

1 - Shingling: convert documents, emails, etc., to sets.
2 - Minhashing: convert large sets to short signatures, while preserving similarity.
3 - Locality-sensitive hashing: focus on pairs of signatures likely to be similar.

[The Big Picture]

Document --> Shingling --->  Minhashing ---> Locality-sensitive Hashing ---> Candidate pairs: those pairs of signatures
                        ^                ^                                   that we need to test for similarity.
                        |                |
                        |    Signatures: short integer vectors that
                        |    represent the sets, and reflect their similarity.

            The set of strings of length k
            that appear in the document

[Shingles]

- A k-shingle (or k-gram) for a document is a sequence of k (consecutive) characters that appears in the document.
  + Blanks are normally considered characters.
  + HTML tags may be considered characters, or be ignored.
  + Example:
    * k = 2
    * doc = 'abcab'
      - Set of 2-shingles = {'ab', 'bc', 'ca'} # no dupes (ab), since we are constructing a set
      - A k in the range of 5 to 10 is generally used
      - Represent a doc by its set of k-shingles

[Shingles and Similarity]

- Documents that are intuitively similar will have many shingles in common.
- Changing a word only affects k-shingles within distance k from the word.
# only the shingles left and right, as well as within, the word, can be affected
- Reordering paragraphs only affects the shingles that cross paragraph boundaries.
- Example:
  + k = 3
  + Old = "The dog which chased the cat"
  + New = "The dog that chased the cat"
    * Only 3-shingles that are replaced are: {'g_w', '_wh', 'whi', 'hic', 'ich', 'ch_', 'h_c'} # _ indicates 'blank' or ' '

[Shingles: Compression Option]

- To compress long shingles, we can hash them to (say) 4 bytes.
  + Called tokens.
- Represent a doc by its tokens, that is, the set of hash values of its k-shingles.
- Two documents could (rarely) appear to have shingles in common, when in fact only the hash-values
  were shared.

## 3 - 2 - Minhashing (25-18).mp4

[Minhashing]

- Jaccard Similarity Measure
- Constructing Signatures

[Jaccard Similarity Measure]

- The Jaccard similarity of two sets is the size of their intersection divided by the size of their union.
  + Sim(C1, C2) = |C1 ∩ C2| / |C1 ∪ C2| # ratio of the intersect / union

[From Sets to Boolean Matrices]

- Rows = elements of the universal set.
  + Example: the set of all k-shingles. # of all documents
- Columns = sets # all inidividual documents
- 1 in row e and column S if and only if e is a member of S. # otherwise that entry is 0
- Column similarity is the Jaccard similarity of the sets of their rows with 1.
  + Where they both are 1 == intersect.
  + Where either one of them is 1 == union.
- Typical matrix is sparse --> many more 0s than 1s. # which also allows for processing in memory

C1  C2
0    1
1    0
1    1
0    0
1    1
0    1

Sim(C1, C2) = 2 / 5 = 0.4

[Four Types of Rows]

- Given columns C1 and C2, rows may be classified as:

    C1   C2
a :  1    1 # TRUE TRUE
b :  1    0 # TRUE FALSE
c :  0    1 # FALSE TRUE
d :  0    0 # FALSE FALSE

- Also, a = number of rows of type a, etc.
- Note Sim(C1, C2) = a / (a + b + c)

[Minhashing]

- Imagine the rows permuted randomly. # shuffled randomly, entire rows, vertically
- Define the minhash function h(C) = the number of the first (in permuted order) row
  in which column C has 1.
- Use several (e.g. 100) independent hash functions to create a signature for each column.
- The signature can be displayed in another matrix - the signature matrix - whose columns
  represent the sets and the rows represent the minhash values, in order for that column.
  # The row contains the result of the chosen minhash function for that row, to each of the columns

[Surprising Property]

- The probability (over all permutations of the rows) that h(C1) = h(C2) is the same as Sim(C1, C2).
  # The probability (over all permuatations of the rows) that the minhash functions return the same hash value is the same as the Jaccard probability
  + Both are a / (a + b + c)
  + Why?
    * Look down the permuted columns C1 and C2 until we see a 1.
    * If it is a type-a row, then h(C1) = h(C2). If a type-b or type-c row, then not.
    # If not the same, the value of the other column will have to wait until it moves down lower,
    # and since the row # increases, its hash value will definitely be higher and won't be the same --> this does not make sense when using modulo!

Thus the probability that two columns will have the same minhash-value equals the probability that the first
row that is not of type d (0 0 | FALSE FALSE) is a type a row (1 1 | TRUE TRUE). That probability is the number
of type a rows divided by the number of rows of any of the types a, b or c (a / (a + b + c))

[Similarity for Signatures]

- The similarity of signatures is the fraction of the minhash functions in which they agree.
  + Thinking of signatures as columns of integers, the similarity of signatures is the fraction
    of rows in which they agree (have the same value).
- Thus the expected similarity of two signatures equals the Jaccard similarity of the columns or
  sets that the signatures represent.
  + And the longer the signatures, the smaller will be the expected error.

Similarity can mean two things:
1 - For columns or sets it is the Jaccard Similarity
2 - For signatures it is the fraction of components in which the signatures agree.

    Row           Input
Permutations      Matrix

1 | 4 | 3   -   1  0  1  0
3 | 2 | 4   -   1  0  0  1
7 | 1 | 7   -   0  1  0  1
6 | 3 | 6   -   0  1  0  1
2 | 6 | 1   -   0  1  0  1
5 | 7 | 2   -   1  0  1  0
4 | 5 | 5   -   1  0  1  0
[1] [2] [3]

The above yields the following 3-signature matrix M (rows), for 4 documents (columns):

2  1  2  1  #  [3]
2  1  4  1  #  [2]
1  2  1  2  #  [1]

Columns:  1-3   2-4   1-2
Col/Col   0.75  0.75  0.0
Sig/Sig   0.67  1.00  0.0

[Implementation of Minhashing]

- We have implemented Minhashing as if we repeatedly permutated the order of the rows randomly.
  + This is not feasible in practice:
    * Suppose 1 billion rows.
    * Hard to pick a random permutation of 1 billion.
    * Representing a random permutation requires 1 billion entries.
    * Accessing rows in permuted order leads to thrashing. # thrashing occurs when a computer's virtual memory subsystem is in a constant state of paging, rapidly exchanging data in memory for data on disk, to the exclusion of most application-level processing

- Here is how we simulate permutation without actually doing so:
  + For each minhash function, we pick a normal sort of hash function which hashes integers to some number of buckets.
  + A good approximation to permuting rows: pick, say, 100 hash functions. # for each minhash function we want to simulate
  + For each column c and each hash function hi, keep a 'slot' M(i, c) # in the signature matrix
  + Intent: M(i, c) will become the smallest value of hi(r) for which column c has 1 in row r.
  # !!!
  # Note that it is possible that hi(r) maps two or more rows to the same value!
  # But if we make the number of buckets hi hashes very large, larger than the number of rows,
  # then the probability of a collision at the smallest value (minhash) is very small and we can then ignore the probability of a collision
  # !!!
    * I.e. hi(r) gives order of rows for ith permutation.

for each row r do begin
  for each hash function hi do
     compute hi(r);
  for each column c
     if c has 1 in row r
       for each hash function hi do
         if hi(r) is smaller than M(i, c) then
           M(i, c) := hi(r);
end;

- Often, data is given by column, not row.
  + Example: columns = documents, rows = shingles.
- If so, sort matrix once so it is by row.

Row  C1  C2
 1   1   0
 2   0   1
 3   1   1
 4   1   0
 5   0   1

 Hash Functions # since we have 2 our signature length will be 2

 h(x) = x % 5
 g(x) = (2x + 1) % 5

- Note that these hash functions work (% 5) because we only have 5 rows!
- If more than 5 h(1) & h(6) would return the same (minimum) value of 1, which is not correct!

           Sig1    Sig2
h(1) = 1    1       ∞
g(1) = 3    3       ∞

h(2) = 1    1       2
g(2) = 3    3       0

h(3) = 1    1       2
g(3) = 3    2       0

h(4) = 1    1       2
g(4) = 3    2       0

h(5) = 1    1       0
g(5) = 3    2       0

So ultimately our signatures do not agree at all (1 != 0 & 2 != 0), whereas the actual
Jaccard similarity is 0.2 (see the input matrix). Adding more hash functions (increasing the
length of the signature) will likely make this error less large.

My interpretation:
- Uses differing hash functions is the same as randomly permuting the rows.
- If you hit a type a row both columns will get the same signature value for that hash function.
- If you do not hit such a row, either type b or c, the hash function will output two differing values.
- Its then left to chance, by having enough hash functions / permutations, you will ultimately end up with about
  the same as the jaccard similarity when comparing signatures (where you look at the % where the min hash values are equal).

# Say one function hashes at minimum at row 100000.
# Thats possible and its hi value will be X.
# If and only if the other set also had a 1 in this row will the other column have the same value (type a row).
# It could have had 1s in the rows above or below, and the min value may be in one of those instead (type b or c row).

# Always looking for the minimum — is that what technically permutes the rows? Yes!

# Thats why set a and b have X as value and can be similar, and why set c and d can have Y.
# Their minimum can lie in various rows.

# Every minhash function basically checks for row type a b c between two columns: if same then a, if not then b or c.

# C1 C2

# Row 1 - 0 1
# Row 2 - 1 1
# Row 3 - 1 0

# Now one minhash function may have the minimum at row 1, then 2, etc.
# Then this minhash function will output a 'b' or 'c' row type.
# Another minhash function may return the minimum value at row 2.
# This function will output an a row type, etc. Having many functions will therefore reduce error by chance.
# By always looking for the minimum, you are permuting the rows.

# 3 - 3 - Locality-Sensitive Hashing (19-24).mp4

[Locality-Sensitive Hashing]

- Focusing on Similar Minhash Signatures
- Other Applications Will Follow

- General idea: Generate from the collection of all elements (signatures in our example)
  a small list of candidate pairs: pairs of elements whose similarity must be evaluated.
  # There will be false negatives: pairs that are similar but won't be checked for similarity
- For signature matrices: Hash columns to many buckets, and make elements of the same bucket pairs. # by ordinary hash functions, not minhash functions
  # A pair becomes a candidate pair if any one or more of the ordinary hash functions puts both signatures in the same bucket
  # We need to tune the number of hash functions and the number of buckets per hash function so that the buckets have relatively few signatures in them
  # That way there are not too many candidate pairs generated - but we cannot use too many buckets or else pairs that are truly similar will not wind up in
  # the same bucket for even one of the hash functions we use.

[Candidate Generation From Minhash Signatures]

- Pick a similarity threshold t, a fraction < 1. # the minimum value of jaccard similarity that we regard as two sets being similar
- We want a pair of columns c and d of the signature matrix M to be a candidate pair if
  and only if their signatures agree in at least fraction t of the rows.
  + I.e., M(i,c) = M(i,d) for at least fraction t values of i.
- Big idea: hash columns of signature matrix M several times.
- Arrange that (only) similar columns are likely to hash to the same bucket.
- Candidate pairs are those that hash at least once to the same bucket

[Partition into Bands]

- We divide the signature matrix M in to a b number of bands
- Each band has a certain number of rows r per band
  - r * b == number signature elements == number of minhash functions
- We will create one hash function from each band!

- Divide matrix M into b bands of r rows.
- For each band, hash its portion of each column to a hash table with k buckets.
  + Make k as large as possible. # so that the hash function almost resembles its identity function (band of the signature has to fully match in order to hash to the same bucket)
- Candidate column pairs are those that hash to the same bucket for >= 1 band.
- Tune b and r to catch most similar pairs, but few nonsimilar pairs.
  # If we make b large and r small then there are lots of bands and therefore lots of opportunities for a pair to wind
  # up in the same bucket. Since r, the width of the band, is small it's not hard for a pair to hash to the same bucket for
  # one of the bands. Thus making b large is good if the similarity threshold is relatively low. --> Easy to be labeled as similar
  # Conversely if we make b small and r large, then it wil be very hard for two signatures to hash to the same bucket for a given
  # band and there are few bands to give them the opportunity to do so. Thus a small number of bands is best if we have a high threshold
  # of similarity. --> Hard to be labeled as similar

[Example - Bands]

- Suppose 100,000 columns
- Signatures of 100 integers
- Therefore, signatures take 40 Mb
- Want all 80%-similar pairs of documents # t
- 5,000,000,000 pairs of signatures can take a while to compare # (100,000 * 99,999)/2 ~ approximated by (100,000^2)/2
- Choose 20 bands of 5 integers/band

[Suppose C1, C2 are 80% Similar]

- Probability C1, C2 identical in one particular band: (0.8)^5 = 0.328 # not that high, but we have 20 chances to hit this at least once!
- Probability C1, C2 are not similar in any of the 20 bands: (1-0.328)^20 = 0.00035
  # not similar being (0.2^5) + (0.2^4 * 0.8) + (0.2^3 * 0.8^2) etc, which == (1-(0.8^5))
- 1 - 0.0035 = the chance that it will be a candidate pair = 0.99965
  + I.e. about 1/3000th (0.99965) of the 80%-similar underlying sets are false negatives. # shown as non-candiates being not similar, yet are similar

[Suppose C1, C2 Only 40% Similar]

- Probability C1, C2 identical in any one particular band: (0.4)^5 = 0.01.
- Probability C1, C2 identical in >= 1 of 20 bands:
  20 * 0.01 = 0.2 # this is not great, 20% false-positives
- But false positives much lower for similarities << 40%.
  + For 20% similarity we only get 1% FPs
  + It is difficult to assess the % of FPs beforehand, given this is dependent on the
    underlying distribution of similarities between all pairs (the massive number of pairs: num_pairs^2/2)
    # For instance if the underlying similarity of all underlying sets is 79%, almost all would be FPs

[Analysis of LSH - What We Want]

- Y-axis: Probability of sharing a bucket
- X-axis: Jaccard Similarity s of two sets
- t: threshold of s at which we deem the sets similar or not
  + We want an S-shaped curve, a vertical line splitting what is before and after t

1 - No chance to become a candidate pair if s < t
2 - Probability should be 100% if s > t to become a candidate pair

[What One Band of One Row Gives You]

- Remember: probability of equal minhash values = Jaccard Similarity
  + The more similar the sets are, the likelier they get to hash the same minimum value!
- One band of one row will give you a straight, positive linear line.
  + At least the probability goes in the right direction.
  + It is linear because even if they are very dissimilar, there is still a chance
    that at the shingles they are similar, for both sets the minhash value is found.
    * This is however unlikely and therefore the probability is low.
    * This probability however becomes linearly larger the more similar the sets become.
    * Hence we get this linear relationship.

[What b Bands of r Rows Gives You]

- When you combine many minhash functions into b bands of r rows each, you begin to
  get the S-curve that you desire (as this resembles a vertical cut-off line on the X-axis)
  + This has greatly reduced FP & FN regions.

1 - If two columns share similarity s (say 60%), the probability of a band of 5 rows being fully equal is: s^r # 0.6^5
2 - Probability that some row r of a band is unequal: 1 - s^r # 1 - 0.6^5 | Remember: S-similar, D-dissimilar: SSSSD + SSSDD + SSDDD etc.
3 - Probability that there are no bands (say there are 10 bands) identical: (1 - s^r)^b # (1 - 0.6^5)^10
4 - At least one band identical: 1 - (1 - s^r)^b # 1 - (1 - 0.6^5)^10

# As b and r get large, this function increasingly resembles a step function, where you take a vertical step on the X-axis.

- The threshold at which the rise occurs, is approximately: t ~ (1/b)^1/r

Example: b = 20; r = 5

 s   -   1-(1-s^r)^b
0.2       0.006
0.3       0.047
0.4       0.186
0.5       0.470
0.6       0.802
0.7       0.975
0.8       0.9996

# Not exactly a step function, but it does get rather steep in the middle & is much
# better than the function we got from a single row.

[LSH Summary] # What we need to do to find sets with a given threshold t of Jaccard Similarity

- Decide on your values of b and r: tune to get almost all pairs with similar signatures,
  but eliminate most pairs that do not have similar signatures.
  + You can use t ~ (1/b)^1/r and see which combinations yields the desired t.
  + Note that the longer your signatures are, the larger b and r can be and the less FPs you will have.
    * However, the longer the signatures, the more space they will take and the more work it will be to
      perform all the minhashing.
- Define a hash function for each band!
- Check that candidate pairs really do have similar signatures.
  + Check if enough elements of the signature agree to match the corresponding minimum threshold t.
- Optional: In another pass through data, check that the remaining candidate pairs really represent similar sets.
  + By looking at the acutal sets (!= signature similarity) and computing the Jaccard Similarity # a / a + b + C
    * It could be that certain candidate pairs using the signature strategy were over-estimated in similarity.
    + This eliminates the final FPs. Unfortunately this does not take care of FNs.

## 3 - 4 - Applications of LSH (11-40).mp4

[Applications of LSH]

1 - Entity Resolution
2 - Fingerprints
3 - Similar News Articles

[1 - Entity Resolution]

- The entity-resolution problem is to examine a collection of records and determine
  which refer to the same entity.
  + Entities could be people, events, etc.
- Typically, we want to merge records if their values in corresponding fields are similar.

[Matching Customer Records]

- I once took a consulting job solving the following problem:
  + Company A agreed to solicit customers for Company B, for a fee.
  + They then argued over how many customers. # A supplied customers to B, hence wanted fees for those referrals.
  + Neither recorded exactly which customers were involved.
    * Figure out how many customers appeared in both databases.

- Each company had about 1 million records describing customers that might have
  been sent from A to B. # 1,000,000 * 1,000,000 = 1 Trilion pairs --> to much to compare.
- Records had name, address, and phone, but for various reasons, they could be
  different for the same person.

- Step 1: Design a measure ('score') of how similar records are:
  # Get 100 points if name/address/phone matches 100% --> 300 is the maximum score.
  # Interestingly only 7000 records received this score, whereas 180.000 candidate pairs were found.
  + E.g. deduct points for small misspellings ('Jeffrey' vs. 'Jeffrey') or same phone
    with different area code.
- Step 2: Score all pairs of records that the LSH scheme identified as candidates; report high scores
  as matches.

- Problem: (1 million)^2 is too man pairs of records to score.
- Solution: A simple LSH. # Not shingling or minhashing
  + Three hash functions: exact values of name, address, phone.
    # First had a bucket for each name
    # Second had a bucket for each address
    # Third had a bucket for each phone number
    * Compare if records are identical in at least one. # Or: compare all columns that share a bucket.
  + Misses similar records with a small difference in all three fields.

- You can use the same hash function to hash each of them to buckets in case of signatures.
  + The values in the bands will always be different so a 100% similar signature of two sets will hash the two signatures/sets b times to the same bucket.
- In this example, simply get all records with each bucket, in SQL conduct a cross-join on their unique ID to get all pairs.
- Do the same for all buckets.
- Union everything into one master table and run a distinct on rows to remove dupes --> you are left with your candidate pairs!

[Aside: Hashing Names, etc.]

- How do we hash strings such as names so there is one bucket for each string?
  + Answer: sort the strings instead. # whenever there is a name change --> next bucket, e.g. RANK() function in SQL
  + If you want to look for similarities with a string (large text for instance), then you move to shingling & minhasing -- but that is now unneccessary.
- You can do this for name, phone number (can be saved as string ;) & address.
- Another option was to use a few milion buckets, and deal with buckets that contain
  several different strings. # You do have the chance that certain strings hash to the same buckets -- meaning they are not similar.
  # But if you make the number of buckets much larger than than the number of records, this probability should be low.

[Aside: Validation of Results]

- We were able to tell what values of the scoring function were reliable in an interesting way.
- Identical records (the ~ 7,000k) had a creation date difference of 10 days. # record created date @ A ten days earlier than B
- We only looked for records created within 90 days of each other, so bogus matches had a 45-day average.
- By looking at the pool of matches with a fixed score, we could compute the average time-difference, say x,
  and deduce that fraction (45-x)/35 of them were valid matches.
  # Say you take a pool of records with similarity score 200
  # Some of them will be accurate matches with an average time difference of 10 days --> given by the 7k right matches
  # Some of them will be false matches with an average time difference of 45 days
  # Within this pool, the average time difference is x. The fraction of values that are valid are therefore: (45-x)/35
  # The formula here is based upon the average difference for bogus matches and equating to 100% when x == what you found in the 7k exact matches
  # For example, if x = 10 then (45-35)/35 = 1 = 100%, as 10 was the accurate difference in days.
  # If x would have been 20, then 5/7 of the matches were valid, meaning 5/7 will have an average difference of 10, and 2/7 an average of 45 days delay!
  # Thus the weighted average of the averages is 20.
- Alas, the lawyers did not think the jury would understand.

[Validation - Generalized]

- Any field not used in the LSH could have been used to validate, provided corresponding values
  were closer for true matches than false.
- Example: if records had a height field, we would expect true matches to be close, false matches to
  have the average difference for random people. # !!!
  # If you have some gold standard of records (exact matches) you can determine what the average measurement
  # error or height difference is between height values of the same entity. You can also determine the average
  # difference in height between two random records. You can then use the above reasoning again to see within
  # your set of similar records which fraction is likely to be correct based on what the observed average height difference is,
  # compared to what it should be according to the golden, exact sample.

## 3 - 5 - Fingerprint Matching (7-07).mp4

[2 - Fingerprint Matching]

- Minutiae
- A New Way of Bucketing

# Usually fingerprint matching is not many:many, rather one:many -- matching one fingerprint with a database of fingerprints
# Still it can help reduce the number of pairs to check

[Fingerprint Comparison]

- Represent a fingerprint by the set of positions of minutiae # set of coordinates on a 2-dimensional space
  + These are features of a fingerprint, e.g., points where two ridges
    come together or a ridge ends.

[LSH for Fingerprints]

- Place a grid on a fingerprint.
  + Normalize so identical prints will overlap.
- Set of grid squares where minutiae are located represents the fingerprint.
- Possibly, treat minutiae near a grid boundary as if also present in adjacent
  grid points.
  # Problem is matrix is not sparse - grid cannot be too fine or it will be unclear where the minutiae belong.
  # Matrix rows: grid squares, columns fingerprint sets. This matrix will therefore not be sparse.
  # That means minhashing will not work very well, each minhash will have relatively few different values so we don't get
  # a good distribution into a large number of buckets when we do the LSH.
  # We need to change things to get LSH to work.

[Applying LSH to Fingerprints]

- Fingerprint = set of grid squares. # those with minutiae
- No need to minhash, since the number of grid squares is not too large.
- Represent each fingerprint by a bit-vector with one position for each square.
  + 1 in only those positions whose squares have minutiae # 0 otherwise --> a bit = 1 or 0
- Example: pick 1024 sets of 3 grid squares (components of the bit vectors), randomly.
- For each set of three squares, we look for all the sets of 3 that have 1 in all 3 elements.
  + Two fingerprints that each have 1 for all three squares in the same 3-set are candidate pairs.
- Funny sort of 'bucketization':
  + Each set of three squares creates one bucket.
  + Prints can be placed in many buckets. # It will be common now that a single fingerprint lands in several different buckets.

[Example: LSH/Fingerprints]

- Suppose typical fingerprints have minutiae in 20% of the grid squares.
- Suppose fingerprints from the same finger agree in at least 80% of the minutiae.
  + It helps that we include nearby grid fields to show minutiae helps with this assumption.
- Probability two random fingerprints each have minutiae in all three squares = (0.2)^6 = 0.000064.
  + For a single finger it is 0.2^3 --> for two fingers it is (0.2^3)^2 == 0.2^6
- Probability two fingerprints from the same finger each have 1s in three given squares:
  + ((0.2)(0.8))^3 = 0.004096 # still tiny, but 64X larger than 0.000064
    * 0.2 = probability there is minutiae
    * 0.8 = probability that there will also then be minutiae in the second fingerprint
    * ^3 = has to happen for three consecutive fingerprints
- Probability for at least one of 1024 sets of three points = 1 - (1 - 0.004096)^1024 = 0.985 # only 1.5% False Negatives
  # What is the probability that at least one of 1024 sets match when coming from the same finger?
  # 1 - 0.004096 == chance they do not match but come from the same finger
  # this has to happen 1024 times --> (1 - 0.004096)^1024
  # then we are left with the odds that none of them matched == 1.5% of False Negatives (since the fingers acutally matched!)
  # Probability that therefore at least one set DID match = 1 - (1 - 0.004096)^1024 == 98.5%
  # Remember that probability there is minutiae is already taken into consideration (20%)!
- But for random fingerprints: 1-(1-0.00064)^1024 = 0.063 # 6.3% FPs

-- Why sometimes shingles + minhasing + lsh and sometimes straight to LSH --

- When comparing documents for similarity you cut them up in shingles. Two documents can match on several shingles:

Take: 'abcdude'

r1 abc  1
r2 bcd  1
r3 cde  0
r5 def  0

In contrast, when matching whole strings:

r1 abcdude 1
r2 bcddude 0
r3 cdedude 0

You can only match in one. Therefore you might as well just sort and rank them and skip minhashing. E.g.

abcdude 1
abcdude 1
abcdude 1
bcddude 2
bcddude 2 # etc.

You can consider the rank the bucket.

Now, why minhashing -- only when you have a lot of shingles, and there is overlap among records, is this neccessary.

''' p.81 from the book:
Sets of shingles are large. Even if we hash them to four bytes each, the space
needed to store a set is still roughly four times the space taken by the document.
If we have millions of documents, it may well not be possible to store all the
shingle-sets in main memory.3
Our goal in this section is to replace large sets by much smaller representations
called “signatures.” The important property we need for signatures is
that we can compare the signatures of two sets and estimate the Jaccard similarity
of the underlying sets from the signatures alone.
'''

## 3 - 6 - Finding Duplicate News Articles (6-08).mp4

[3 - Similar News Articles]

- A New Way of Shingling
- Bucketing by Length

[Application: Same News Articles]

- The Political-Science Department at Stanford asked a team from CS to help them with the
  problem of identifying duplicate, on-line news articles.
- Problem: the same article, say from the Associated Press, appears on the Web site of many
  newspapers, but looks quite different.

[News Articles]

- Each newspaper srrounds the text of the article with:
  + Its own logo and text.
  + Ads.
  + Perhaps links to other articles.
- A newspaper may also 'crop' the article (delete parts).
- The team came up with its own solution, that included shingling, but not
  minhashing or LSH.
  + A special way of shingling that appears quite good for this application.
  + LSH substitute: candidates are articles of similar length.

[Enter LSH]

- I told them the story of minhashing + LSH.
- They implemented it and found it faster for similarities below 80%.
  + Aside: That is no surprise. When similarity is high, there are better
    methods, as we shall see.
- Their first attempt at minhasing was very inefficient.
- They were unaware of the importance of doing the minhashing row-by-row.
- Since their data was column-by-column, they needed to sort once before minhashing.

[Specialized Shingling Technique]

- The team observed that news articles have a lot of stop words while ads do not.
  + "Buy Sudzo" (Ad) vs. "I recommend that you buy Sudzo for your laundry." (Article)
- They defined a shingle to be a stop word and the next two following words.

[Why it Works]

- By requiring each shingle to have a stop word, they biased the mapping from documents to shingles
  so it picked more shingles from the article than from the ads.
- Pages with the same article, but different ads, have higher Jaccard similarity than those with the
  same ads, yet different articles.

## 3 - 7 - Distance Measures (22-39).mp4

[Theory of LSH]

- Distance Measures
- LS Families of Hash Functions
- S-Curves
- LSH for Different Distance Measures

[Distance Measures]

- Generalized LSH is based on some kind of 'distance' between points.
  + Similar points are 'close'.
- Jaccard similarity is not a distance; 1 minus Jaccard similarity is.
- Two major classes of distance measures:
  1 - Euclidian
  2 - Non-Euclidian

[Euclidian vs. Non-Euclidian]

- A Euclidian space has some number of real-valued dimensions and 'dense' points.
  + The typical ones are 1 & 2-dimensional, but it can have any number of dimensions.
  + There is a notion of 'average' of two points.
  + A Euclidian distance is based on the locations of points in such a space.
- Any other space is Non-Euclidian.
  + Distance measures for non-Euclidian spaces are based on properties of points, but
    not their 'location' in a space.

[Axioms of a Distance Measure]

- d is a distance measure if it is a function from pairs of points to real numbers such that:
  1 - d(x,y) >= 0 # it never has a negative value
  2 - d(x,y) = 0 if x = y # distance only 0 when the points are the same
  3 - d(x,y) = d(y,x) # distance is symmetric
  4 - d(x,y) <= d(x,z) + d(z,y) [triangle inequality] # direct distance between x,y is always smaller than when adding an arbitrary point z to go through

[Some Euclidian Distances]

- L2 norm: d(x,y) = square root of the sum of the squares of the differences between x and y in each dimension.
  + The most common notion of 'distance'.
  + [3, 3, 3] & [1, 1, 1] --> (2^2 + 2^2 + 2^2)^0.5 =~ 3.46
- L1 norm: sum of the differences in each dimension.
  + Manhattan distance = distance if you had to travel along coordinates only.
  + [3, 3, 3] & [1, 1, 1] --> (2 + 2 + 2) == 6
- L∞ norm: d(x,y) = the maximum of the differences between x and y in any dimension.
  + [3, 4, 6] & [10, 9, 8] --> [7, 5, 2] == 7
  + Note: the maximum is the limit as r goes to ∞ of the Lr norm: what you get by taking the
    r^th power of the differences, summing and taking the r^th root.

The general formula is: sum(distances^r)^1/r --> both L1 (r = 1) & L2 (r = 2) adhere to this one!
- Raising r larger and larger causes the largest difference to dominate the sum, all the other ones become negligible.
- When you take the r^th root of that sum you essentially are taking the r^th root of the r^th power of the largest difference,
  which gives you back essentially just the largest of the differences.

[Non-Euclidian Distances]

- Jaccard distance for sets == 1 - Jaccard similarity.
- Cosine distance for vectors = angle between the vectors.
- Edit distance for strings = minimum number of inserts and deletes to change one string into another.
  + Sometimes we allow a mutation as one edit, where a mutation changes one character to another: abc --> abd.
  + Without mutations this would have costed 2 edits: abc --> ab --> abd.
- Hamming Distance for bit vectors (of the same lenght) = the number of positions in which they differ.

[Example: Jaccard Distance]

- Consider x = {1, 2, 3, 4} and y = {1, 3, 5}
- Size of intersection = 2; size of union = 5, Jaccard similarity (not distance) = 2 / 5.
- Jaccard Distance : 1 - 2/5 == 3/5

[Why J.D. Is a Distance Measure]

# ∩ - intersect
# ∪ - union

- d(x,y) >= 0 because |x ∩ y| <= |x ∪ y| # the intersect is always smaller or equal to the union
- d(x,x) = 0 because |x ∩ x| = |x ∪ x| # union == intersect e.g. 1 - 2/2 = 0
  + And if x ≠ y, then the size of x∩y is strictly less than the size of x∪y
- d(x,y) = d(y,x) because union and intersection are symmetric
- d(x,y) <= d(x,z) + d(z,y) [triangle inequality] --> tricker:

#  Jaccard Distance xz         Jaccard Distance yz          Jaccard Distance xy
[1 - (|x ∩ z| / |x ∪ z|)] + [1 - (|y ∩ z| / |y ∪ z|)] >= [1 - (|x ∩ y | / |x ∪ y|)]

- Remember: |a ∩ b| / |a ∪ b| = probability that minhash(a) = minhash(b).
- Thus, 1 - |a ∩ b| / |a ∪ b| = probability that minhash(a) ≠ minhash(b).
- Claim: prob[minhash(x) ≠ minhash(y)] <= prob[minhash(x) ≠ minhash(z)] + prob[minhash(z) ≠ minhash(y)]
- Proof: whenever minhash(x) ≠ minhash(y), at least one of minhash(x) ≠ minhash(z) and minhash(z) ≠ minhash(y) must be true.
  + A ≠ B, B ≠ C, C ≠ A
  + IF C ≠ A --> then either A ≠ B OR B ≠ C (or both) must be TRUE, because if both are FALSE THEN A = B and B = C, which must mean that A = C. # by "transitivity of equals"

[Cosine Distance]

- Think of a point as a vector from the origin ([0, 0, ..., 0]) to its location.
- Two points vectors make an angle, whose cosine is the normalized dot-product of the vectors: [p1 * p2] / [|p2| * |p1|]
  + The dot product is the sum of the product of the corresponding components: p1 * p2
  + Which you divide by the lengths of the two vectors: |p2||p1| # L1 norm of the point of the head of the vector to the origin
  # That is it is the square root of the sum of the squared components of the vector.
  + Example: p1 = [0,0,1,1,1]; p2 = [1,0,0,1,1]
    * p1 * p2 = 2; |p1| = |p2| = 3^0.5 # there are 3 1's, hence (1^2 * 3)^0.5 = sqrt(3)
    * 3^0.5 * 3^0.5 = 3
    * cos(θ) = 2/3; hence θ is about 48 degrees.

[Why C.D. Is a Distance Measure]

- We only think about direction of a vector, not size.
  + A vector in its negation (the opposite direction) can be thought of as the same vector.
- d(x,x) = 0 because arccos(1) = 0 # arccos == angle with itself is 0
- d(x,y) >= 0 because any two intersecting vectors make an angle in the rage 0 to 180 degrees.
- d(x,y) = d(y,x) by symmetry.
- Triangle inequality: physical reasoning. If I rotate an angle from x to z and then from z to y,
  I cannot rotate less than from x to y.

[Edit Distance]

- The edit distance of two strings is the number of inserts and deletes of characters needed to
  turn one into the other.
- An equivalent definition: d(x,y) = |x| + |y| - 2|LCS(x,y)|
  + |x| = length of string x; |y| = length of string y; |LCS(x,y)| = length of LCS.
  + LCS = longest common subsequence = any longest string obtained both by deleting from x and deleting from y.
  + Example: x = 'abcde'; y = 'bcduve';
    * Turn x into y by deleting 'a', then inserting 'u' and 'v' after 'd'.
      - Edit distance: 3
    * Or, LCS(x,y) = 'bcde' # they do not have to be contiguous!
      - |x| + |y| - 2|LCS(x,y)| == 5 + 6 - 2*4 == 3 == edit distance

[Why Edit Distance is a Distance Measure]

- d(x,x) = 0 because 0 edits suffice.
- d(x,y) >= 0: no notion of negative edits.
- d(x,y) = d(y,x) because insert/delete are inverses of each other.
- Triangle inequality: changing x to z and then to y is one way to change x to y.
  + The minimum number of edits required to make that transformation is the sum of the edit distances from
    x to z and from z to y. This sequence of edits is one of the possible ways to transform x to y. So the
    total number of edits is at least the edit distance from x to y.

[Hamming Distance]

- Hamming distance: is the number of positions in which bit-vectors differ.
- Example: p1 = 10101; p2 = 10011
  + d(p1, p2) = 2 because the bit vectors differ in the 3rd and 4th positions.

[Why Hamming Distance Is a Distance Measure]

- d(x,x) = 0 since no positions differ.
- d(x,y) >= 0 since strings cannot differ in a negative number of positions.
- d(x,y) = d(y,x) by symmetry of "different from".
- Triangle inequality: changing x to z and then to y is one way to change x to y.
  + The sum of these two numbers cannot be less than the number of bits you have to flip to turn x to y directly.

## 3 - 8 - Nearest Neighbor Learning (11-39).mp4

[Large Scale Machine Learning]

- Nearest Neighbours

[Supervised Learning]

- Would like to do prediciton:
  estimate a function f(x) so that y = f(x)
- Where y can be:
  + Real number: Regression
  + Categorical: Classification
  + Complex object:
    * Ranking of items, Parse tree, etc.
- Data is labeled:
  + Have many pairs {(x,y)}
    * x .. vector of binary, categorical, real valued features
    * y .. class ({+1, -1}, or a real number)

[Large Scale Machine Learning]

- We will talk about the following methods:
  + k-Nearest Neighbor (Instance based learning)
  + Support Vector Machines
  + Decision trees
- Main question:
  + How to efficiently train (build a model/find model parameters)?

[Instance Based Learning]

- Instance based learning
- Example: Nearest neighbour
  + Keep the whole training dataset: {(x, y)}
  + A query example (vector) q comes
  + Find closes example(s) x*
  + Predict y*
- Works both for regression and classification
  + Collaborative filtering is an example of k-NN classifier
    * Find k most similar people to user x that have rated movie y
    * Predict yx of x as an average of yk

[1-Nearest Neighbor]

- To make Nearest Neighbor work we need 4 things:
  + Distance metric:
    * Euclidian
  + How many neighbors to look at?
    * One
  + Weighting function (optional):
    * Unused (only one neighbor)
  + How to fit with the local points?
    * Just predict the same output as the nearest neighbor (just use one NN in this example)

[k-Nearest Neighbor]

- Distance metric:
  + Euclidian
- How many neighbors to look at?
  + k
- Weighting function (optional):
  + Unused (if you want to weigh certain neighbors higher)
- How to fit with the local points?
  + Just predict the average output among k nearest neighbors

[How to find nearest neighbors?]

- Given: a set P of n points in R^d
- Goal: Given a query point q
  + NN: Find the nearest neighbor p of q in P
  + Range search: Find one/all points in P within distance r from q

# Can use LSH for large scale data!

## 3 - 9 - Frequent Itemsets (29-50).mp4

[Frequent Itemsets]

- The Market-Basket Model
- Association Rules
- A-Priori Algorithm

[The Market-Basket Model] # Essentially a representation for a many:many relationship between two concepts, which we call items & baskets

- A large set of items, e.g., things sold in a supermarket.
- A large set of baskets, each of which is a small set of the items, e.g. the
  things one customer buys on one day.

[Support]

- Simplest question: find sets of items that appear 'frequently' in the baskets.
- Support for itemset I = the number of baskets containing all items in I.
  + Sometimes given as a percentage.
- Given a support threshold s, sets of items that appear in at least s baskets
  are called frequent itemsets.

[Example: Frequent Itemsets]

- Items = {milk, coke, pepsi, beer, juice}
- Support = 3 baskets

B1 = {m, c, b}    B2 = {m, p, j}
B3 = {m, b}       B4 = {c, j}
B5 = {m, p, b}    B6 = {m, c, b, j}
B7 = {c, b, j}    B8 = {b, c}

- Frequent itemsets:
  + Singletons: {m}, {c}, {b}, {j}
  + Doubletons: {m, b}, {b, c}, {c, j}
  + Tripletons: None

[Applications]

1 - Items = products; baskets = sets of products someone bought in one trip to the store.
  - Example application - given that many people buy beer and diapers together:
    + Run a sale on diapers; raise price of beer.
      * Causality is important: if you put beer on sale, people coming in to get beer do
        not necessarily have a baby!
  - Only useful if many buy diapers & beer.
    + This is essential for brick-and-mortar stores, not online stores.
      * Else BM stores are wasting a lot of time optimising for a sale that is rare.
      * Online stores are different, which can cater to the long tail.

2 - Baskets = sentences; items = documents containing those sentences.
  - Items that appear together too often could represent plagiarism.
  - Notice items do not have to be 'in' baskets. # in fact in this example it is the other way around
    + But it is better if baskets have small numbers of items, while items can be in large numbers of baskets.
    + This is because the work it takes to process baskets is quadratic to the number of items contained in it (n * (n-1)/2 ~> n^2/2).

3 - Baskets = documents; items = words.
  - Unusual words appearing together in a large number of documents, e.g. 'Brad' and 'Angelina', may
    indicate an interesting relationship.

[Scale of the Problem]

- WalMart sells 100,000 items and can store billions of baskets.
- The Web has billions of words and many billions of pages.

[Association Rules]

- If-then rules about the contents of baskets.
- {i1, i2, ..., ik} --> means: "if a basket contains all of i1, ..., ik then it is likely to contain j."
  + The degree to which this association is likely is called the confidence of the rule.
- Confidence of this association rule is the probability of j given i1, ..., ik.
  + That is, the fraction of the baskets with i1, ..., ik that also contain j.

[Example: Confidence]

B1 = {m, c, b}    B2 = {m, p, j}
B3 = {m, b}       B4 = {c, j}
B5 = {m, p, b}    B6 = {m, c, b, j}
B7 = {c, b, j}    B8 = {b, c}

- An association rule: {m, b} --> c
  + B1 & B6 do have coke
  + B3 & B5 do not have coke
    * Confidence = 2/4 = 50%

[Finding Association Rules]

- Question: "find all association rules with support >= s and confidence >= c".
  + Note: 'support' of an association rule is the support of the set of items on the left.
- Hard part: finding the frequent itemsets.
  + Note: if {i1, i2, ..., ik} --> has high support and confidence, then both {i1, i2, ..., ik} [s]
    and {i1, i2, ..., ik, j} [c*s] will be 'frequent'.

1. Find all sets with support at least cs.
2. Find all sets with support s. # subset of 1.
3. If {i1, i2, ..., ik, j} has support at least cs, see which subsets
   missing one element have support at least s.
   - Take j to be the missing element
4. {i1, i2, ..., ik} --> j is an acceptable asociation rule if {i1, i2, ..., ik} has
   support s1 >= s, {i1, i2, ..., ik, j} has support s2 >= cs, and s2/s1, the confidence
   of the rule is at least c. # the fraction of the baskets that also contain j

[Computation Model]

- Typically, data is kept in flat files.
- Stored on disk.
- Stored basket-by-basket.
- Expand baskets into pairs, triples, etc. as you read baskets.
  + Use k nested loops to generate all sets of size k.

[File Organization]

Example: items are positive integers, and boundaries between baskets are -1.

[Computational Model (2)]

- The true cost of mining disk-resident data is usually the number of disk I/Os.
- In practice, algorithms for finding frequent itemsets read the data in passes - all baskets read in turn.
- Thus, we measure the cost by the number of passes an algorithm takes.

[Main-Memory Bottleneck]

- For many frequent-itemset algorithms, main memory is the critical resource.
- As we read baskets, we need to count something, e.g., occurrences of pairs.
- The number of different things we can count is limited by main memory.
- Swapping counts in/out is a disaster.

[Finding Frequent Pairs]

- The hardest problem often turns out to be finding frequent pairs.
  + Why? Often frequent pairs are common, frequent triples are rare.
    * Why? Support threshold is usually set high enough that you don not get
      too many frequent itemsets.
- We will concentrate on pairs, then extend to larger sets.

[Naïve Algorithm]

- Read file once, counting in main memory the occurrences of each pair.
  + From each basket of n items, generate its n*(n-1)/2 pairs by two nested loops.
- Fails if (num. items)^2 exceeds main memory.
  + Remember: num. of items can be 100K (Wal-Mart) or 100B (Web Pages).

[Details of Main-Memory Counting]

- Two approaches:
  1. Count all pairs, using a triangular matrix.
  2. Keep a table of triples [i, j, c] = "the count of the pair of items {i, j} is c".
- (1) requires only 4 bytes/pair. # 4 bytes per possible pair
  + Note: always assume integers are 4 bytes.
- (2) requires 12 bytes (4 * 3), but only for those pairs with count > 0. # 12 bytes per existing pair

# Hence if more than 1/3 of possible pairs are present in at least one basket, you prefer the triangular matrix!
# If you expect fewer than 1/3 of possible pairs to be present in the data, then you should go for the tabular approach.

[Triangular-Matrix Approach]

- Number items 1, 2, ...
  + Requires table of size O(n) to convert item names to consecutive integers.
    * Example is a hash table that contains the key (integer) and the original name.
- Count {i, j} only if i < j. # same as an ordered list of length 2
- Keep pairs in the order {1,2}, {1,3}, ..., {1,n}, {2,3}, {2,4}, ..., {2,n},
  {3,4}, ..., {3,4}, ..., {n-1, n}
- Find pair {i,j}, where i<j, at the position:

(i - 1) * (n - (i/2)) + j - i

{3, 5}
n = 10

(3 - 1) * (10 - 3/2) + 5 - 3 = 19
# 9 pairs before it starting with 1 (lowest member)
# 8 pairs before it starting with 2 (lowest member)
# 1 more pair {3, 4} that comes before {3, 5}
# Hence its number is 1 + (9 + 8 + 1) = 19

- Total number of pairs n(n - 1)/2; total bytes about 2*n^2 (4 * n^2/2)

[Details of Tabular Approach]

- Total bytes used is about 12p, where p is the number of pairs that actually occur.
  + Beats triangular marix if at most 1/3 of possible pairs actually occur.
- May require extra space for retrieval structure, e.g. a hash table.

## 3 - 10 - A-Priori Algorithm (13-07).mp4

- Monotonicity of 'Frequent' # An itemset cannot be frequent unless all its subsets are frequent
- Candidate Pairs
- Extension to Larger Itemsets

[The A-Priori Algorithm]

- A two-pass approach called a-priori limits the need for main memory.
- Key idea: monotonicity: if a set of items appears at least s times, so does every subset of s.
- Contrapositive for pairs: if item i does not appear in s baskets, then no pair including i can
  appear in s baskets.

- Pass 1: Read baskets and count in main memory the occurrences of each item.
  + Requires only memory proportional to the number of items.
- Items that appear at least s times are the frequent items.
  + Note that you only need to count up to s, this saves memory!
- Pass 2: Read baskets again and count in main memory only those pairs of which
  were found in Pass 1 to be frequent.
  + If only half the items are frequent, we only have to count ~ 1/4 of all possible pairs.
- Requires memory proportional to square of frequent items only (for counts), plus a list
  of the frequent items (so you know what must be counted).

[Detail for A-Priori]

- You can use the triangular matrix method with n = number of frequent items.
  + May save space compared with storing triples.
- Trick: re-number frequent items, 1, 2, ... and keep a table relating new numbers to
  original item numbers. # To make the triangular matrix work - by recounting it is as small as possible

[Frequent Triples, Etc.]

- For each k, we construct two sets of k-sets (sets of size k):
  + Ck = candidate k-sets = those that might be frequent sets (support >= s) based on
    information from the pass for k-1.
  + Lk = the set of truly frequent k-sets.

[Passes Beyond Two] # Here we go!

- C1 = all items
- In general, Lk = members of Ck with support >= s.
  + Requires one pass.
- Ck+1 = (k+1)-sets, each k of which is in Lk.

1 - Assume k = 3.
2 - We might find {1, 3, 5} is in Lk (L3).
3 - Now look at each item whose number is higher than any in the set. # You don't go down: {1, 3, 5} + 2 would be covered by {1, 2, 3} + 5!
    + The first higher item is 6 in this case: {1, 3, 5, 6}
4 - We already know its subset {1, 3, 5} is frequent / in L3.
    + We have to test 3 other sets:
      * {3, 5, 6}
      * {1, 5, 6}
      * {1, 3, 6}
    + If all three of these are in L3, then {1, 3, 5, 6} would go into C4 # candidate-4
5 - Then we continue with 7: {1, 3, 5, 7}, etc., up to 9!

Memory Requirements

- At the k`th pass, you need space to count each member of Ck.
- In realistic cases, because you need fairly high support, the number of
  candidates of each size drops, once you get beyond pairs.

[Improvements to A-Priori]

- Park-Chen-Yu Algorithm
- Multistage and Multihash
- Single-Pass Approximate Algorithms

[PYC Algorithm]

- During Pass 1 of A-Priori, most memory is idle.
  + The essence of the PCY algorithm is to take advantage of the fact that when the set of
    items is modest in size we have lots of space in memory during the first pass which we
    do not seem to need.
  + We can in fact do some counting that will give us a smaller candidate set of pairs that will
    help us on the second pass where the main memory is often the critical resource since we do
    not have to store as many counts.
- Use that memory to keep counts of buckets into which pairs of items are hashed.
  + JUST THE COUNT, not the pairs themselves!
  + Therefore the space needed for a bucket is small, 4-bytes is surely sufficient.
    * We might get away with 2-bytes if the support threshold is less than 2^16.
- For each basket, enumerate all its pairs, hash them, and increment the resulting
  bucket count by 1.
  + So you do two things:
    * Count items (item + count)
    * Count pairs (count only)
- A bucket is frequent if its count is at least the support threshold.
- If a bucket is not frequent, no pair that hashes to that bucket could possibly be
  a frequent pair.
- On Pass 2, we only count pairs that hash to frequent buckets. # and when both items are frequent
- The counts retrieved in Pass 1, are transformed into a bitmap, and used to cross-check
  in Pass 2. # 1 means frequent, 0 means not frequent
- Since we go from 4-bytes per count to 1 bit, we have a 32/1 compression! # 4 * 8 = 32

[Pass 1: Memory Organization]

- Space to count each item.
  + One (typically) 4-byte integer per item.
- Use the rest of the space for as many integers, representing buckets, as we can.

[PCY Algorithm - Pass 1]

FOR (each basket) {
  FOR (each item in the basket)
    add 1 to item`s count;
  FOR (each pair of items) {
    hash the pair to a bucket;
    add 1 to the count for that bucket
  }
}

[Observations About the Buckets]

1. A bucket that a frequent pair hashes to is surely frequent.
   + In the subsequent pass, its value in the bitmap will be 1, and we will store the count.
   + We cannot use the hash table to eliminate any member of this bucket.
   + Unfortunately many infrequent pairs may hash to the same bucket.
     # e.g. s = 3 and {1, 2} & {3, 4} & {5, 6} all happen to hash to bucket 1;
     # which is now considered frequent and its member counts will all be stored in Pass 2
    * We hope at least one member of the pair is not a frequent item so it
      is excluded by the other role that both i and j must be frequent.
    * However, sometimes you will have a pair that is not frequent, but its items are, and it hashes
      to a frequent bucket, even though it is not frequent.
2. Even without any frequent pair, a bucket can be frequent. # e.g. 5 non-frequent pairs hash to the bucket and s = 5
   + Again, nothing in the bucket can be eliminated.
3. But in the best case, the count for a bucket is less than the support s.
   + Now, all pairs that hash to this bucket can be eliminated as candidates,
     even if the pairs consists of two frequent items.

[PCY Algorithm - Between Passes]

- Replace the buckets by a bit-vector (the 'bitmap'):
  + 1 means the bucket is frequent; 0 means it is not
- 4-byte integers are replaced by bits, so the bit-vector requires 1/32 of memory.
- Also, decide which items are frequent and list them for the second pass.

[PCY Algorithm - Pass 2]

- Count all pairs {i,j} that meet the conditions for being a candidate pair:
  1. Both i and j are frequent items. # For A-Priori this it he ONLY condition
  2. The pair {i,j} hashes to a bucket number whose bit in the bit vector is 1.
  # The second condition is introduced by PCY

[Memory Details]

- Buckets require a few bytes each.
  + Note: we do not have to count past s.
    * If count is less than s, increment with 1, else do nothing.
  + Number of buckets is O(main-memory size).
    * Its size is a reasonable fraction of main memory - typically half or a quarter.
- On a second pass, a table of (item, item, count) triples is essential.
  + The pairs that are exempted because they hashed to an infrequent bucket are scattered
    all over the place and cannot be organised into a nice triangular array.
    * Thus, the hash table must eliminate 2/3 of the candidate pairs for PCY
      to beat A-Priori # as a triangular matrix uses 1/3 the space compared to the tabular approach

[Multistage Algorithm]

- Key idea: After Pass 1 of PCY, rehash only those pairs that qualify for Pass 2 of PCY.
- On middle pass, fewer pairs contribute to buckets, so fewer false-positives - frequent
  buckets with no frequent pairs.

[Multistage - Pass 3]

- Count only those pairs {i, j} that satisfy these candidate pairs conditions:
  1. Both i and j are frequent items. # A-Priori condition
  2. Using the first hash function, the pair hashes to a
     bucket whose bit in the first bit-vector is 1. # PCY condition
  3. Using the second hash function, the pair hashes to a bucket whose bit
     in the second bit-vector is 1. # Multistage condition

[Important Points]

1. The hash functions have to be independent.
# Otherwise each would report the same frequent bucket and we have no advantage
2. We need to check both hashes on the third pass.
   + If not, we would wind up counting pairs of frequent items
     that hashed first to an infrequent bucket but happened to hash
     second to a frequent bucket. # We'd lose the info of the first check!

[Multihash]

- Key idea: use several independent hash tables on the first pass.
- Risk: halving the number of buckets doubles the average count. # You are OoM. E.g. 1 table with 1000 buckets, is 2 tables with 500 each, hence average count goes up.
  We have to be sure most buckets will still not reach count s.
  # Remember that PCY works if the average bucket count is much less than support threshold s.
  # That way many buckets are infrequent and we can eliminate many candidate pairs for the second Pass.
- If so, we can get a benefit like multistage, but in fewer passes.
  + Example:
    * Say we have one table where the average count of items per bucket is s/10.
    * If we then split the table into two half-sized hash tables, the average count will become s/5.
    * Each table will still have mostly infrequent buckets.
    * But now you have 2 chances to eliminate a pair that is not actually frequent.
    * In this case we get much of the benefit of a 2-stage algorithm, but we only use
      2 passes and not 3.

## 3 - 12 - All or Most Frequent Itemsets in 2 Passes (14-40) [Advanced]

[All (or most) Frequent Itemsets In <= 2 Passes]

- Simple Algorithm
- Savasere-Omiecinski-Navathe (SON) Algorithm
- Toivonen`s Algorithm

[Simple Algorithm]

- Take a random sample of the market baskets.
- Run A-Priori or one of its improvements (for sets of all sizes, not just pairs) in main memory,
  so you do not pay for disk I/O each time you increase the size of itemsets.
- Use as your support threshold a suitable, scaled-back number.
  + Example: if your sample is 1/100 of the baskets, use s/100 as your support threshold instead of s.
- Optionally, verify that your guesses are truly frequent in the entire data set by a second pass. # Removing FPs
  + But you do not catch sets frequent in the whole but not in the sample. # cathcing FNs
    * To avoid FNs: use a smaller threshold, e.g. s/125 instead of s/100, which helps catch more truly frequent itemsets.
      - But requires more space.

[SON Algorithm]

- Repeatedly read small subsets of the baskets into main memory and perform the first pass of the
  simple algorithm on each subset.
- An itemset becomes a candidate if it is found to be frequent in any one or more subsets of the baskets.

[SON Algorithm - Pass 2]

- On a second pass, count all the candidate itemsets and determine which are frequent in the entire set.
- Key 'monotonicity' idea: an itemset cannot be frequent in the entire set of baskets unless it is frequent
  in at least one subset.

[SON Algorithm - Distributed Version]

- This idea lends itself to distributed data mining.
- If baskets are distributed among many nodes, compute local frequent itemsets at each node,
  then distribute the candidates from each node (so all nodes have an exhaustive list of all candidates).
- Each node then counts all the candidate itemsets for their respective baskets.
- Finally, accumulate the counts of all candidates in a reduce tasks.

[Toivonens Algorithm]

- Start as in the simple algorithm, but lower the threshold slightly for the sample.
  + Example: if the sample is 1% of the baskets, use s/125 as the support threshold rather than s/100.
  + Goal is to avoid missing any itemset that is frequent in the full set of baskets.
- Add to the itemsets that are frequent in the sample the negative border of these itemsets.
- An itemset is in the negative border if it is not deemed frequent in the sample, but all its
  immediate subsets are.
  + None of these should be frequent in the whole set, since we chose a support threshold that is
    significantly lower (s/125) than the actual threshold (s/100).
  + But if one or more of the item sets in the negative border turn out to be frequent in the whole,
    then we have to believe that the sample was not representative and we need to repeat the entire process.
  + This is why there is no limit to how many passes Toivonens algorithm may need to produce an answer.

[Example: Negative Border]

- {A, B, C, D} is in the negative border if and only if:
  1. It is not frequent in the sample, but
  2. all of its immediate subsets: {A,B,C}, {B,C,D}, {A,C,D} and {A,B,D} are (frequent).
- {A} is in the negative border if and only if it is not frequent in the sample.
  + Because the empty set is always frequent. # Which is the subset of {A}
    * Unless there are fewer baskets than the support threshold (silly case).
- In a second pass, count all candidate frequent itemsets from the first pass, and also count the sets
  that comprise their negative border.
- If no itemset from the negative border turns out to be frequent, then the candidates found to be
  frequent in the whole data are exactly the frequent itemsets.
- What if we find that something in the negative border is actuallyl frequent?
  + We must start over with another sample!
- Try to choose the support threshold so the probability of failure is low, while the number of
  itemsets checked on the second pass fits in main-memory.

[Theorem]

- If there is an itemset that is frequent in the whole, but not in the sample, then there is a member
  of the negative border for the sample that is frequent in the whole.

[Proof]

- Suppose not; i.e.;
  1. There is an itemset S frequent in the whole but not frequent in the sample, and
  2. Nothing in the negative border is frequent in the whole.
- Let T be a smallest subset of S that is not frequent in the sample.
- T is frequent in the whole (S is frequent + monotonicity).
- T is in the negative border (else not 'smallest').

-- Week 3 --

## 4 - 1 - Community Detection in Graphs- Motivation (5-44).mp4

[Community Detection in GraΠphs]

- Facebook friend graph segmentation

[Overlapping Communities]

- Non-overlapping vs. overlapping communities

[Communities as Tiles!]

- What is the structure of community overlaps:
  + Edge density in the overlaps is higher!

## 4 - 2 - The Affiliation Graph Model (10-04).mp4

[The Affiliation Graph Model]

[Plan of Attack]

1 - Given a model, we generate the network: Generative model for network --> network
2 - Given a network, find the 'best' model: network --> Generative model for network

[Model of networks]

- Goal: Define a model that can generate networks (1)
  + The model will have a set of 'parameters' that we will later want
    to estimate (and detect communities)
- Q: Given a set of nodes, how do communities 'generate' edges of the network?

[Community-Affiliation Graph]

- AGM: Affiliation Graph Model: a generative model B(V, C, M, {pc}) for graphs:
  + Set of nodes V, set of communities C, set of memberships M
  + Each community c has a single probability pc (likelihood its member nodes connect with one another)
  + Later we fit the model to networks to detect communities

[AGM: Generative Process]

- AGM generates the links: For each
  + For each pair of nodes in community A, we connect them with prob. pa
  + The overall edge probability is:
    # Probability that node u and v are connected:
    * P(u,v) = 1 - Π(1 - p*c)
              c ∈ Mu ∩ Mv
      - I.e. the product over all the communities they have in common.
      - Mu .. set of communities node u belongs to.
      - Think of this as an 'OR' function: if at least 1 community says 'YES' we create an edge
      - ∈ == element of!
        + x ∈ A ∩ B if and only if:
          * x ∈ A and
          * x ∈ B.
      - The higher the number of shared communities, the higher the probability an
        edge will be created between the two nodes.
      - If u & v share no communities: P(u,v) = ε # some small probability

[AGM: Flexibility]

- AGM can express a variety of community structures:
  + Non-overlapping
  + Overlapping
  + Nested
    * Al we need is the underlying structure (community sets) and AGM can generate the network.

## 4 - 3 - From AGM to BIGCLAM (8-48).mp4

[From AGM to BigCLAM]

- Relaxation: Memberships have strengths
  + FuA: The membership strength of node u to community A
    * FuA = 0: no membership # value always non-negative
  + Each community A links nodes independently:
    * PA(u,v) = 1 - exp(-FuA * FvA)

[Factor Matrix F]

- Community membership strength matrix F
  + Rows: Nodes
  + Columns: Communities
- Basically each row is a vector of community membership strenghts for a certain node.
- PA(u,v) = 1 - exp(-FuA * FvA)
  + Probability of connection is proportional to the product of strengths.
    * Notice: If one node does not belong to the community (FuC = 0) then P(u,v) = 0.
- Prob. that at least one common community C links the nodes:
  + P(u,v) = 1 - Πc(1 - Pc(u,v))
    * 1 - Prob. they connect == prob that they do not connect
       --> product over all communities == probability that they do not connect in any of the communities
       --> 1 - this probability == probability they connect in at least one community.
    * Where Pc == PA(u,v) = 1 - exp(-FuA * FvA)

- Community A links nodes u, v independently:
  + PA(u,v) = exp(-FuA * FvA)
- Then prob. at least one common c links them:
  + P(u,v) = 1 - Πc(1 - Pc(u,v))
           = 1 - exp(-∑c (FuC * FvC))
           = 1 - exp(-Fu * Fv^T) # we replace the ∑ with the dot-product of vectors Fu & Fv

## 4 - 4 - Solving the BIGCLAM (9-19).mp4

[Solving the BigCLAM]

[BiGCLAM: How to find F] # affiliation matrix F

- Task: Given a network G(V,E), estimate F # G == network, V == Nodes, E == Network as found
  + Find F that maximizes the likelihood:
    * arg maxF Π p(u,v) Π (1 - p(u,v))
            (u,v) ∈ E   (u,v) ∉ E
      - Where P(u,v) = 1 - exp(-Fu * Fv^T)
      - Many times we take the logarithm of the likelihood,
        and call it log-likelihood: l(F) = log P(G|F) # Probability of a Graph given F
- Goal: Find F that maximizes l(F):

l(F) = ∑ (log(1 - exp(-Fu * Fv^T))) - ∑ (Fu * Fv^T)
 (u,v) ∈ E                        (u,v) ∉ E

# -> Either you have a model (the Factor Matrix F), and you can generate the network.
# -> Or you have a network, and want to fit the model (you create the Factor Matrix F)!
#    + For this one, you maximize l(F) using gradient ascent.

# Todo:
# - Generate some network in node matrix.
# - Use gradient ascent to generate F (nodes / communities)
# - Use F to generate node matrix --> if similar, success!
# With AGM you have Pc, however with BiGCLAM no more? Replaced with continuousness of degree of affiliation with C's!
# - Therefore the problem is now different... ANyhow re-do video 3 and code it up.

(u,v) ∈ E --> this must be stored in the node graph?





