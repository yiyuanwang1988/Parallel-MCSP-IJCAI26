#pragma once
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <limits>
#include "graph.hpp"   // RRSplitGraph、arguments 等
#include "../mpi/Message.hpp"
#include "../mpi/Solver.hpp"
#include "../mpi/Worker.hpp"

static std::vector<int> degInC1;
static std::vector<int> degInC2;
static std::vector<int> EqClass;
static std::vector<int> index_right;
static std::vector<int> isfixvertex;


class rrs:public Solver{
    public :
        rrs() = default;
        struct VtxPair {
            int v, w;
            VtxPair(int v, int w):
            v(v),
            w(w){};
        };

        struct Bidomain {
            int l, r, left_len, right_len;
            bool is_adjacent;
            Bidomain(int l, int r, int left_len, int right_len, bool is_adjacent):    
            l(l),
            r(r),
            left_len (left_len),
            right_len (right_len),
            is_adjacent (is_adjacent) { };
        };

        rrsGraph::Graph g0_sorted,g1_sorted;
        void init();
        void run(Worker &worker) override;

        struct TaskPackage {
            rrsGraph::Graph g0, g1;
            std::vector<int> vertexs;
            std::set<int> unvertexSet;
            unsigned int masterid;
        };

        // std::vector<VtxPair> mcs(const rrsGraph::Graph& g0, const rrsGraph::Graph& g1);
    private:
};