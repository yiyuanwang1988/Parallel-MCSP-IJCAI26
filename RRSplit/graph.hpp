#include <limits.h>
#include <stdbool.h>

#include <vector>

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>
using ui = unsigned int;
//using std::cout; using std::endl; using std::vector;

class rrsGraph{
    public:
        struct Graph {
        int n;
        std::vector<std::vector<unsigned int>> adjmat;
        std::vector<unsigned int> label;
        std::vector<unsigned int> degree;                 
        std::vector<std::vector<unsigned int>> adjlist;
        Graph() = default;
        Graph(unsigned int n);
    };
    static Graph induced_subGraph(struct Graph& g, std::vector<int> vv);

    static Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);

    static void set_adjlist(struct Graph & g);

    static void GetEqClass(Graph & g, std::vector<int>& EqClass);

    
};



