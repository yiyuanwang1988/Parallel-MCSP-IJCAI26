#pragma once
#include <limits.h>
#include <stdbool.h>

#include <vector>

class dsbGraph{
    public:
        struct Graph {
            int n;
            std::vector<std::vector<unsigned int>> adjmat;
            std::vector<unsigned int> label;
            Graph() = default;
            Graph(unsigned int n);
        };
        static Graph induced_subgraph(struct Graph& g, std::vector<int> vv);
        static Graph readGraph(char* filename, char format, bool directed, bool edge_labelled, bool vertex_labelled);
};