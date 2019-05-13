# Classical N-Body Problem (Barnes-Hut)

This algorithm is used to quickly approximate the net forces exerted by n bodies on an arbitrary mass in space.  This allows you to quickly solve the n-body problem for every large instances, and my implementation takes advantage of multi-core systems.

## Barnes-Hut Approximation for the N-Body Problem
In classical mechanics, the n-body problem consists of predicting the motion of celestial objects interacting through gravity.

In a computational approach to making these predictions, the forces at a given time are calculated, and then numerical integration is used to predict the positions and velocities a short time in the future.  Repeating this process makes it possible to generate accurate predictions for celestial motion far into the future even when closed-form solutions are not possible.

### Parallel Barnes-Hut Implementation
The 'naive' computational approach to the n-body problem computes the gravitational force induced by every body to every other body. This results in *O(n<sup>2</sup>)* work to compute all forces for a single timestep. The Barnes-Hut algorithm reduces the complexity to *O(nlog(n))*, which makes it tractable to closely approximate solutions for problem sizes that would otherwise be prohibitively expensive. 

The central idea is grouping together bodies that are sufficiently close to each other and sufficiently far from the force calculation point, and approximating the group of bodies as one body with the total group mass located at the group center of mass. 'Sufficiently close' and 'sufficiently far' are defined by a simulation parameter called the multipole acceptance criterion (MAC), usually denoted by &theta;, which is a ratio of the group 'diameter' *d* and the distance *r* between the force calculation point and the group center of mass. In the figure below, diagram *B* represents a Barnes-Hut approximation of diagram *A*.

![MAC](https://upload.wikimedia.org/wikipedia/commons/thumb/e/e2/Barnes_hut.svg/376px-Barnes_hut.svg.png)

Generally, there is a tradeoff between accuracy and speedup. The lower the MAC, the more accurate but more expensive the computation, because only smaller groups (and hence a greater number of total groups) can be used in the approximation.

#### Quadtrees
To find suitable sets of bodies to group together, a Quadtree data structure is used. 

In calculating the force on a particular body, the bodies within a square are grouped together, and if approximating these bodies as a single body meets the MAC (smaller than &theta;), the approximation is used.  Otherwise, the contributions from the four sub-squares are summed together.

![Quadtree example](https://encrypted-tbn1.gstatic.com/images?q=tbn:ANd9GcQzciGQ-YBx_XnetitdFd6x4M91lmMTWbGym9O2U1FpTonfXsc4)

Given the quadtree, parallelizing the force calculation is a simple parallel loop.  Parallelizing the construction of the tree, however, is a little more challenging.  At first, the tree construction may seem like a curious thing to optimize. As the number of bodies grows, however, the share of computational time due to the construction increases.  In addition, for other applications it may be important to compute the gravitational effect of a large number of bodies, not on each other, but on a smaller set of bodies, perhaps even in real-time.  The real-time nature of the problem makes a fast Barnes-Hut approximation necessary, and it would be nice to compute the quadtree quickly so as to be sure to have the tree ready when needed.

The key to fast merging of the trees is to avoid having them overlap.  If a subtree is present in one tree but entirely absent in an another, then the subtree can be added simply by marking it as a child of the appropriate node.  We want to have as much of the merge happen in this way as possible.  To achieve this, we use a MortonKey ordering of the bodies.

#### Morton Ordering
Morton order is a space-filling curve that follows a 'Z' or 'N' pattern microscopically and macroscopically. This order helps us construct our quadtree efficiently.

The Morton order is found by forming Morton 'keys' associated with each point, and then sorting these keys. Each key can be determined by taking integer representations of the x and y coordinates, and interleaving the bits of these integer representations. If our integer representations of x and y for a particular point are 10001 and 01110, respectively, then with interleaving we get the Morton key of 0110101001. After sorting the Morton keys, the points (quadtree leaves) are in the order they would appear in the quadtree.


## Benchmarking
<pre><code>
$ make teragen terametrics
$ if [ ! -e data_bench.dat ]
$ then
$    echo "* Creating datafile"
$    mpirun ./teragen -c 250000 -f data_bench.dat
$    sleep 10
$ else
$    echo "* Reusing existing data file"
$ fi
$ mpirun ./terametrics -f data_bench.dat -c 20
</pre></code>

System: 1 x dual-socket E5-2680v4 node, 16 processors per node

Setup time: 0.923773 seconds
Query time: 0.040223 microseconds per query
