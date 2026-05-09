# A Parallel Framework for the Maximum Common Induced Subgraph Problem
We propose a parallel framework for solving the maximum common induced subgraph problem. The framework parallelizes three state-of-the-art sequential algorithms: McSplitDAL, McSplitDSB, and RRSplit.

## 1.compile
To compile the source code and generate the executable `./mpi/main`, use the following command:
```
mpic++ -O3 -march=native -ffast-math -pedantic -Wall -std=c++17 -I. mpi/main.cpp mpi/Master.cpp mpi/Worker.cpp  mpi/Node.cpp mpi/Tree.cpp McSplitDSB/mcsp.cpp McSplitDSB/graph.cpp  RRSplit/mcsp.cpp RRSplit/graph.cpp McSplitDAL/mcsp.cpp McSplitDAL/graph.cpp -pthread -o mpi/main
```
## 2. Usage
To execute the code, use the following command format:
```
mpirun -np 8 ./mpi/main min_max -P ./images-CVIU11/patterns/pattern1 -T ./images-CVIU11/targets/target1 -l -q -t 1800 -S 3
```
The above parameters are explained as follows:

-   `-np`: the number of processes
-   `-P`: input pattern graph
-   `-T`: input target graph
-   `-l`: format of datasets
-   `-q`: quiet output
-   `-t`: time limit (seconds)
-   `-S`: specifies the algorithm to execute (1 for ParaDAL, 2 for ParaDSB, 3 for ParaRRS)





