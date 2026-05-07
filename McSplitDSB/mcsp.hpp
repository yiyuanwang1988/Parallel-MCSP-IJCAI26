#ifndef MCSP_HPP
#define MCSP_HPP

#include "graph.hpp"   // RRSplitGraph、arguments 等
#include "../mpi/Message.hpp"
#include "../mpi/Solver.hpp"
#include "../mpi/Worker.hpp"
using std::vector;
using std::cout;
using std::endl;

class dsb:public Solver{
    public:
        dsb() = default;
        struct VtxPair {
            int v;
            int w;
            VtxPair(int v, int w): v(v), w(w) {};
        };

        struct Bidomain {
            int l, r;
            int left_len, right_len;
            bool is_adjacent;
            int  size;
            bool is_valid;

            Bidomain(int l, int r, int left_len, int right_len,
                    bool is_adjacent, int size, bool is_valid): 
                    l(l), r(r), left_len(left_len), right_len(right_len),
                    is_adjacent(is_adjacent), size(size), is_valid(is_valid) {};
        };

        dsbGraph::Graph g0_sorted,g1_sorted;
        // std::vector<VtxPair> mcs(dsbGraph::Graph& g0, const dsbGraph::Graph& g1,
        //                         std::vector<int> vertexsSet,
        //                         std::set<int> un);
        struct TaskPackage {
            dsbGraph::Graph g0, g1;
            std::vector<int> vertexs;
            std::set<int> unvertexSet;
            unsigned int masterid;
        };
        
        void init();
        void run(Worker &worker);
    private:
        // void show(const std::vector<VtxPair>& current,
        //         const std::vector<Bidomain>& domains,
        //         const std::vector<int>& left,
        //         const std::vector<int>& right);

        // bool check_sol(const Graph& g0, const Graph& g1,
        //             const std::vector<VtxPair>& solution);

        // void initializeArray();

        // unsigned int calc_bound(const Graph& g0, const Graph& g1,
        //                         std::vector<Bidomain>& domains,
        //                         const std::vector<int>& left,
        //                         const std::vector<int>& right,
        //                         const std::vector<VtxPair>& incumbent,
        //                         const std::vector<VtxPair>& current);

        // int find_min_value(const std::vector<int>& arr, int start_idx, int len);

        // int select_bidomain(const std::vector<Bidomain>& domains,
        //                     const std::vector<int>& left,
        //                     int current_matching_size);

        // int partition(std::vector<int>& all_vv, int start, int len,
        //             const std::vector<unsigned int>& adjrow);

        // std::vector<Bidomain> filter_domains(const std::vector<Bidomain>& d,
        //                                     std::vector<int>& left,
        //                                     std::vector<int>& right,
        //                                     const Graph& g0, const Graph& g1,
        //                                     int v, int w, bool multiway);

        // int index_of_next_smallest(const std::vector<int>& arr,
        //                         int start_idx, int len, int w);

        // void remove_vtx_from_left_domain(std::vector<int>& left,
        //                                 Bidomain& bd, int v);

        // void remove_bidomain(std::vector<Bidomain>& domains, int idx);

        // void solve(const Graph& g0, const Graph& g1,
        //         std::vector<VtxPair>& incumbent,
        //         std::vector<VtxPair>& current,
        //         std::vector<Bidomain>& domains,
        //         std::vector<int>& vertexs,
        //         std::vector<int>& left,
        //         std::vector<int>& right,
        //         unsigned int matching_size_goal);
};
#endif // MCSP_HPP