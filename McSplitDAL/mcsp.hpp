#pragma once
#include <vector>
#include <set>
#include "graph.hpp"   // RRSplitGraph、arguments 等
#include "../mpi/Message.hpp"
#include "../mpi/Solver.hpp"
#include "../mpi/Worker.hpp"
using gtype = double;

class dal:public Solver{
    public:
        dal() = default;
        struct VtxPair {
            int v;
            int w;
            VtxPair(int v, int w): v(v), w(w) {}
        };

        struct Bidomain {
            int l,        r;        // start indices of left and right sets
            int left_len, right_len;
            bool is_adjacent;
            Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):
                    l(l),
                    r(r),
                    left_len (left_len),
                    right_len (right_len),
                    is_adjacent (is_adjacent) { };
        };
        dalGraph::Graph g0_sorted,g1_sorted;

        void init();
        void run(Worker &worker);
        
        struct TaskPackage {
            dalGraph::Graph g0, g1;
            std::vector<int> vertexs;
            std::set<int> unvertexSet;
            unsigned int masterid;
        };
};