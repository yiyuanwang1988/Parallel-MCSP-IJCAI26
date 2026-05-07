#include "mcsp.hpp"
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <limits>

using std::vector;
using std::cout;
using std::endl;
static double SolveTime = 0;
static int now_ans = 0;
int isreuse = 0;

static unsigned long long nodes{0};
static unsigned long long cut_branches{0};
static std::queue<rrs::TaskPackage> g_taskQ;
static std::mutex g_qmtx,mtx,data_mtx;
static std::condition_variable g_qcv;
static bool check_sol(const rrsGraph::Graph& g0, const rrsGraph::Graph& g1,
                      const vector<rrs::VtxPair>& solution) {
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (const auto& p : solution) {
        if (used_left[p.v] || used_right[p.w]) {
            cout << "Vertex repeated in the solution" << endl;
            return false;
        }
        used_left[p.v] = true;
        used_right[p.w] = true;
        if (g0.label[p.v] != g1.label[p.w]) {
            cout << "Label mismatch" << endl;
            return false;
        }
        for (size_t j = 0; j < solution.size(); ++j) {
            const auto& q = solution[j];
            if (g0.adjmat[p.v][q.v] != g1.adjmat[p.w][q.w]) {
                cout << "Topological mismatch: ("
                    << p.v << "," << q.v << ") in g0 = " << g0.adjmat[p.v][q.v]
                    << ", (" << p.w << "," << q.w << ") in g1 = "
                    << g1.adjmat[p.w][q.w] << endl;
                return false;
            }
        }
    }
    return true;
}

static void print_current(const vector<rrs::VtxPair>& current, const int id, int flag) {
    if (flag) cout << "in id " << id << " current (" << current.size() << "):";
    else cout << "out id " << id << " current (" << current.size() << "):";
    for (const auto &p : current) {
        cout << ' ' << '(' << p.v << ',' << p.w << ')';
    }
    cout << endl;
}

static void print_domains(const vector<rrs::Bidomain>& domains,
                          const vector<int>& left,
                          const vector<int>& right,
                          const int id) {
    cout << "id " << id << " domains (" << domains.size() << "):" << endl;
    for (size_t di = 0; di < domains.size(); ++di) {
        const auto &bd = domains[di];
        cout << " domain " << di << ": l=" << bd.l << " r=" << bd.r
             << " left_len=" << bd.left_len << " right_len=" << bd.right_len
             << " is_adj=" << (bd.is_adjacent ? 1 : 0) << endl;
        cout << "  left: ";
        for (int i = 0; i < bd.left_len; ++i) cout << left[bd.l + i] << ' ';
        cout << endl;
        cout << "  right: ";
        for (int i = 0; i < bd.right_len; ++i) cout << right[bd.r + i] << ' ';
        cout << endl;
    }
}

static void print_graph_struct(const rrsGraph::Graph& g, const std::string& name, const int id) {
    cout << "id " << id << " Graph `" << name << "` n=" << g.n << "\n";
    cout << " labels: ";
    for (int i = 0; i < g.n; ++i) cout << g.label[i] << ' ';
    cout << "\n";

    cout << " adjacency matrix:\n";
    for (int i = 0; i < g.n; ++i) {
        for (int j = 0; j < g.n; ++j) cout << g.adjmat[i][j] << ' ';
        cout << "\n";
    }

    if (!g.degree.empty() && !g.adjlist.empty()) {
        cout << " degrees and adjacency lists:\n";
        for (int i = 0; i < g.n; ++i) {
            cout << "  " << i << ": deg=" << g.degree[i] << " -> ";
            for (unsigned int k = 0; k < g.degree[i]; ++k) cout << g.adjlist[i][k] << ' ';
            cout << "\n";
        }
    } else {
        cout << " adjlist/degree not set (call set_adjlist first)\n";
    }
}




static int calc_bound(const vector<rrs::Bidomain>& domains) {
    int bound = 0;
    for (const auto& bd : domains)
        bound += std::min(bd.left_len, bd.right_len);
    return bound;
}

static int find_min_value(const vector<int>& arr, int start_idx, int len) {
    int min_v = std::numeric_limits<int>::max();
    for (int i = 0; i < len; ++i)
        if (arr[start_idx + i] < min_v) min_v = arr[start_idx + i];
    return min_v;
}

static int select_bidomain(const vector<rrs::Bidomain>& domains,
                           const vector<int>& left,
                           int current_matching_size) {
    int min_size = std::numeric_limits<int>::max();
    int min_tie  = std::numeric_limits<int>::max();
    int best     = -1;
    for (size_t i = 0; i < domains.size(); ++i) {
        const auto& bd = domains[i];
        if (arguments.connected && current_matching_size > 0 && !bd.is_adjacent) continue;
        int len = (arguments.heuristic == min_max) ?
                  std::max(bd.left_len, bd.right_len) :
                  bd.left_len * bd.right_len;
        if (len < min_size) {
            min_size = len;
            min_tie  = find_min_value(left, bd.l, bd.left_len);
            best     = static_cast<int>(i);
        } else if (len == min_size) {
            int tie = find_min_value(left, bd.l, bd.left_len);
            if (tie < min_tie) {
                min_tie = tie;
                best    = static_cast<int>(i);
            }
        }
    }
    return best;
}

static int partition(vector<int>& arr, int start, int len,
                     const vector<unsigned int>& adjrow) {
    int i = 0;
    for (int j = 0; j < len; ++j) {
        if (adjrow[arr[start + j]]) {
            std::swap(arr[start + i], arr[start + j]);
            ++i;
        }
    }
    return i;
}

static int partition_right(vector<int>& arr, int start, int len,
                            const vector<unsigned int>& adjrow) {
    int i = 0;
    for (int j = 0; j < len; ++j) {
        if (adjrow[arr[start + j]]) {
            std::swap(index_right[arr[start+i]],index_right[arr[start+j]]);
            std::swap(arr[start + i], arr[start + j]);
            ++i;
        }
    }
    return i;
}

static int partition_sparse(vector<int>& arr, int start, int len,
                             int degree, const vector<unsigned int>& adjlist) {
    int i = 0; int pos;
    for (int k = 0; k < degree; ++k) {
        pos = index_right[adjlist[k]]; // 0;   // 你原有逻辑
        if (pos >= start && pos < start + len) {
            std::swap(index_right[arr[start+i]],index_right[arr[pos]]);
            std::swap(arr[start + i], arr[pos]);
            ++i;
        }
    }
    return i;
}


static int index_of_next_smallest(const vector<int>& arr, int start_idx, int len, int w,int v) {
    int idx = -1;
    int smallest = INT_MAX;
    for (int i = 0; i < len; ++i) {
        int val = arr[start_idx + i];
        if (val > w && val < smallest) {
            smallest = val;
            idx = i;
        }
    }
    return idx;
}

static void remove_vtx_from_left_domain(vector<int>& left, rrs::Bidomain& bd, int v) {
    int i = 0;
    while (left[bd.l + i] != v) ++i;
    std::swap(left[bd.l + i], left[bd.l + bd.left_len - 1]);
    --bd.left_len;
}

static void remove_bidomain(vector<rrs::Bidomain>& domains, int idx) {
    domains[idx] = domains.back();
    domains.pop_back();
}




static vector<rrs::Bidomain> filter_domains(const vector<rrs::Bidomain>& d,
                                         vector<int>& left,
                                         vector<int>& right,
                                         const rrsGraph::Graph& g0,
                                         const rrsGraph::Graph& g1,
                                         int v, int w, bool& best_match) {
    vector<rrs::Bidomain> new_d;
    new_d.reserve(d.size());
    unsigned int ccount = 0;
    for (const auto& old_bd : d) {
        int l = old_bd.l, r = old_bd.r;
        int left_len  = partition(left, l, old_bd.left_len, g0.adjmat[v]);
        int right_len = (old_bd.right_len > static_cast<int>(g1.degree[w])) ?
                        partition_sparse(right, r, old_bd.right_len, g1.degree[w], g1.adjlist[w]) :
                        partition_right(right, r, old_bd.right_len, g1.adjmat[w]);
        int left_len_noedge  = old_bd.left_len  - left_len;
        int right_len_noedge = old_bd.right_len - right_len;

        if ((left_len == 0 && right_len == 0) ||
            (left_len_noedge == 0 && right_len_noedge == 0) ||
            old_bd.left_len == 0) ++ccount;

        if (left_len_noedge && right_len_noedge)
            new_d.emplace_back(l + left_len, r + right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent);
        if (left_len && right_len)
            new_d.emplace_back(l, r, left_len, right_len, true);
    }
    best_match = (ccount == d.size());
    return new_d;
}

/* ---------------- solve 核心搜索 ---------------- */
static int solve(const rrsGraph::Graph& g0,
                  const rrsGraph::Graph& g1,
                  vector<rrs::VtxPair>& incumbent,
                  vector<rrs::VtxPair>& current,
                  vector<rrs::Bidomain>& domains,
                  vector<int> vertexs,        // 优先点集，可 pop()
                  vector<int> unvertexs,
                  vector<int> data_unvertexs,
                  vector<int>& left,
                  vector<int>& right,
                  unsigned int matching_size_goal,
                  unsigned int level, const int id, int& current_best)
{
    // -------------- 1. 原有时限/停止检查 --------------
    //if (GetSolveTime() > cutoff_time) return -1;
    if (mcs_should_stop.load()) return -1;

    // -------------- 2. 更新最优解 --------------
    if (current.size() > incumbent.size()) {
        incumbent = current;
        now_ans = current.size();
        SolveTime = GetSolveTime();
    }
    if(current.size() > current_best){
        current_best = current.size();
    }
    // -------------- 3. 剪枝 --------------
    unsigned int bound = current.size() + calc_bound(domains);
    if (bound <= getAns() || bound <= incumbent.size() || bound < matching_size_goal) {
        if(current_best < bound) current_best = bound;
        cut_branches++;
        return -1;
    }
    if(datas.size()&&current.size()){
        std::unique_lock<std::mutex> lock(data_mtx, std::try_to_lock);
        if(lock.owns_lock()){
            vector<int> tmp = data_unvertexs;
            int v = current.front().v,w = current.front().w;
            Data data = {v,w,0,tmp.size(),tmp};
            auto it = lower_bound(datas.begin(),datas.end(),data,cmpdata);
            while(it!=datas.end() && it->v == data.v && it->w == data.w ){
                
                int res = checkdata(it->vs,data.vs);
                if(res >= 0){
                    if(res+it->ans <= max(getAns(),now_ans)){
                        return data.v;
                    }
                }
                it++;
            }
        }
    }
    //if (arguments.big_first && incumbent.size() == matching_size_goal) return;

    // -------------- 4. 选 v --------------
    int v;
    rrs::Bidomain* bd = nullptr;
    int bd_idx = -1;

    if (!vertexs.empty()) {
        /* 优先从 vertexs 弹出第一个顶点 */
        v = vertexs.front();
        vertexs.erase(vertexs.begin());
        /* 找到它所在的 bidomain */
        int found = 0;
        for (int i = 0; i < (int)domains.size() && !found; ++i) {
            for (int j = 0; j < domains[i].left_len; ++j) {
                if (left[domains[i].l + j] == v) {
                    bd = &domains[i];
                    bd_idx = i;
                    found = 1;
                    break;
                }
            }
        }
        if (!found) return -1;          // 顶点已不在任何域，直接回溯
    } else {
        /* 退回原逻辑：先选域再选最小值 */
        bd_idx = select_bidomain(domains, left, current.size());
        if (bd_idx == -1) return -1;
        bd = &domains[bd_idx];
        v = find_min_value(left, bd->l, bd->left_len);
    } 

    // -------------- 5. 把 v 从域里删掉 --------------
    remove_vtx_from_left_domain(left, domains[bd_idx], v);
    // -------------- 6. 原分支逻辑（完全不变） --------------
    int w = -1, idx = -1;
    bd->right_len--;
    for (rrs::VtxPair &a : current)
        if (EqClass[a.v] == EqClass[v] && w < a.w) w = a.w;

    if (w > 0) {
        int count_left = 0, count_right = 0;
        for (int i = bd->left_len; i >= 0; --i)
            if (EqClass[left[bd->l + i]] == EqClass[v]) ++count_left;
        for (int i = bd->right_len; i >= 0; --i)
            if (right[bd->r + i] > w) ++count_right;
        if(bd->left_len<=bd->right_len && count_left>count_right){
            if(bound+count_right-count_left<=incumbent.size() || bound+count_right-count_left<= getAns()){
                return -1;
            }
        }
        if(bd->left_len>bd->right_len && (bd->right_len-count_right)>(bd->left_len-count_left)){
            if(bound+(bd->left_len-count_left)-(bd->right_len-count_right)<=incumbent.size()||bound+(bd->left_len-count_left)-(bd->right_len-count_right) <= getAns()){
                return -1;
            }
        }
        // if ((bd->left_len <= bd->right_len && count_left > count_right) ||
        //     (bd->left_len > bd->right_len && (bd->right_len - count_right) > (bd->left_len - count_left)))
        //     if (bound + std::min(count_right - count_left,
        //                          bd->left_len - count_left - bd->right_len + count_right) <= incumbent.size())
        //         return;
    }

    nodes++;
    bool best_match = false;

    int len = bd->right_len;
    for (int i = len; i >= 0; --i) {
        idx = index_of_next_smallest(right, bd->r+bd->right_len-len, len + 1, w,v);
        if (idx == -1) break;
        w = right[bd->r + idx];
        std::swap(index_right[w], index_right[right[bd->r + bd->right_len]]);
        right[bd->r + idx] = right[bd->r + bd->right_len];
        right[bd->r + bd->right_len] = w;
        
        auto new_domains = filter_domains(domains, left, right, g0, g1, v, w, best_match);

        
        current.emplace_back(v, w);
        int t = solve(g0, g1, incumbent, current, new_domains, vertexs, unvertexs, data_unvertexs, left, right, matching_size_goal, level + 1, id ,current_best);
        current.pop_back();

        if(t!=-1){
            if(t!=v) return t;
        }
        if(current.empty()&&unvertexs.size() < 10){
        // if(new_datas.size()<10){
            if(current_best < bound) current_best = bound;
            Data data = {v,w,current_best,unvertexs.size(),unvertexs};
            {
                std::lock_guard<std::mutex> lock(mtx); // 加锁
                new_datas.push_back(data);
            }
            current_best = 0;           
        }


        if (best_match || bound <= incumbent.size() || bound<= getAns()) return -1;
    }

    bd->right_len++;
    for (int i = 0; i < bd->left_len; ++i)
        if (EqClass[left[bd->l + i]] == EqClass[v]) {
            std::swap(left[bd->l + i], left[bd->l + bd->left_len - 1]);
            --bd->left_len; --i;
        }
    if (bd->left_len == 0) remove_bidomain(domains, bd_idx);

    unvertexs.insert(lower_bound(unvertexs.begin(),unvertexs.end(),v),v);
    int t=solve(g0, g1, incumbent, current, domains, vertexs, unvertexs, data_unvertexs, left, right, matching_size_goal, level + 1, id ,current_best);
    return t;
}


static void erase_g0_vertices(rrsGraph::Graph & g0, const rrsGraph::Graph & g1, 
                              vector<rrs::Bidomain>& domains,
                              vector<int>& left,
                              const std::vector<int> v, const std::set<int>& uv, const int id)
{
   
    const int old_n = g0.n;
    std::vector<int>  remap;
    if (uv.empty()) {
        remap.resize(old_n);
        for (int v = 0; v < old_n; ++v) {
            remap[v] = v;
        }
    }
    else {
        remap.resize(old_n, -1);
        for (int v = 0; v < old_n; ++v) {
            if (uv.count(v) == 0)            // 没被删掉就保留
                remap[v] = v;
        }
    }
 

    /* ---------- 1. 从 left 数组里删除 uv 中的顶点 ---------- */
    auto new_end = std::remove_if(left.begin(), left.end(),
                                  [&](int v) { return uv.count(v); });
    left.erase(new_end, left.end());
    

    /* ---------- 2. 重建 domains：重新统计各域的 left 范围 ---------- */
    vector<rrs::Bidomain> new_domains;
    int pos = 0;                       // 当前扫描位置
    for (auto& bd : domains) {
        int l_old = bd.l;              // 原起始下标
        int r_old = bd.r;              // 右侧起始不变
        int old_left_len = left.size();

        /* 收集该域在 left 中仍存活的顶点 */
        int cnt = 0;
        for (int i = 0; i < old_left_len; ++i) {
            int v = left[l_old + i];   // 注意：此时 left 已删完点
            if (!uv.count(v)) ++cnt;
        }
        if (cnt == 0) continue;        // 空域直接丢弃

        /* 把存活顶点紧凑地放到 [pos, pos+cnt) 区间 */
        int write = pos;
        for (int i = 0; i < old_left_len; ++i) {
            int v = left[l_old + i];
            if (!uv.count(v)) left[write++] = v;
        }
        new_domains.emplace_back(pos, r_old, cnt, bd.right_len, bd.is_adjacent);
        pos += cnt;
    }
    domains.swap(new_domains);


    /* 3. 更新Graph */
    rrsGraph::Graph new_g0(old_n);
    for (int old_v = 0; old_v < old_n; ++old_v) {
        if (remap[old_v] == -1) {//continue;
            new_g0.label[old_v] = 2;
            for (int old_u = 0; old_u < old_n; ++old_u) {
                new_g0.adjmat[old_v][old_u] = 2;
                //cout << new_g0.adjmat[old_v][old_u] << "(" << old_v << "," << old_u << ") ";
            }
            //cout << endl;  
        }
        else {
            new_g0.label[old_v] = g0.label[old_v];
            for (int old_u = 0; old_u < old_n; ++old_u) {
                if (remap[old_u] == -1) new_g0.adjmat[old_v][old_u] = 2;
                else new_g0.adjmat[old_v][old_u] = g0.adjmat[old_v][old_u];
                //cout << new_g0.adjmat[old_v][old_u] << "(" << old_v << "," << old_u << ") ";
            }
            //cout << endl;   
        }
        
    }
    // rebuild degree / adjlist for new_g0
    rrsGraph::set_adjlist(new_g0);

    // move new graph into g0
    g0 = std::move(new_g0);

    rrsGraph::GetEqClass(g0, EqClass);
    index_right.clear();
    for(int i=0;i<g1.n;++i) index_right.push_back(i);
    
}


static vector<rrs::VtxPair> mcs(rrsGraph::Graph & g0, const rrsGraph::Graph & g1, std::vector<int> v, std::set<int> uv, const int id) {
    vector<int> left, right;
    vector<rrs::Bidomain> domains;

    // 标签交集初始化逻辑（你原有）
    std::set<unsigned int> left_labels, right_labels, labels;
    for (unsigned int l : g0.label) left_labels.insert(l);
    for (unsigned int l : g1.label) right_labels.insert(l);
    std::set_intersection(left_labels.begin(), left_labels.end(),
                          right_labels.begin(), right_labels.end(),
                          std::inserter(labels, labels.begin()));

    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();

        for (int i=0; i<g0.n; i++)
            if (g0.label[i]==label)
                left.push_back(i);
        for (int i=0; i<g1.n; i++)
            if (g1.label[i]==label)
                right.push_back(i);

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
    }

   

    erase_g0_vertices(g0, g1, domains, left, v, uv, id);

    vector<rrs::VtxPair> incumbent;
    vector<int> uvs(uv.begin(),uv.end());
    sort(uvs.begin(),uvs.end());
    int current_best = 0;
    if (arguments.big_first) {
        for (unsigned int goal = g0.n; goal > 0; --goal) {
            auto left_copy = left; auto right_copy = right; auto domains_copy = domains;
            vector<rrs::VtxPair> current;
            solve(g0, g1, incumbent, current, domains_copy, v, uvs, uvs, left_copy, right_copy, goal, 1, id ,current_best);
            if (incumbent.size() == goal) break;
            if (!arguments.quiet) cout << "Upper bound: " << goal - 1 << endl;
        }
    } else {
        vector<rrs::VtxPair> current;
        solve(g0, g1, incumbent, current, domains, v, uvs, uvs, left, right, 1, 1, id, current_best);
    }
    //std::cout << "id " << id << " debug rrs   " << check_sol(g0, g1, incumbent) << " " << incumbent.size()  << " " << GetSolveTime() << " isreuse " << isreuse << " totaltime " << GetSolveTime() << " " << thread_should_stop.load() << " " << mcs_should_stop.load() << std::endl;
    return incumbent;
}
static vector<int> calculate_degrees(const rrsGraph::Graph& g) {
    vector<int> degree(g.n, 0);
    for (int v=0; v<g.n; v++) {
        for (int w=0; w<g.n; w++) {
            unsigned int mask = 0xFFFFu;
            if (g.adjmat[v][w] & mask) degree[v]++;
            if (g.adjmat[v][w] & ~mask) degree[v]++;  // inward edge, in directed case
        }
    }
    return degree;
}

static int sum(const vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}
void rrs::init(){
    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    struct rrsGraph::Graph g0 = rrsGraph::readGraph(arguments.filename1, format, arguments.directed,
    arguments.edge_labelled, arguments.vertex_labelled);
    struct rrsGraph::Graph g1 = rrsGraph::readGraph(arguments.filename2, format, arguments.directed,
    arguments.edge_labelled, arguments.vertex_labelled);
    if(g0.n>g1.n){
        struct rrsGraph::Graph gg = g0;
        g0 = g1;
        g1 = gg;
    }  

    vector<int> g0_deg = calculate_degrees(g0);
    vector<int> g1_deg = calculate_degrees(g1);
    for(int i=0;i<g1.n;++i) index_right.push_back(i);

    vv0.resize(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    bool g1_dense = sum(g1_deg) < g1.n*(g1.n-1);
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
    return !g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
    });
    vv1.resize(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    bool g0_dense = sum(g0_deg) < g0.n*(g0.n-1);
    std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
    return !g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
    });

    g0_sorted = rrsGraph::induced_subGraph(g0, vv0);
    g1_sorted = rrsGraph::induced_subGraph(g1, vv1);

    rrsGraph::set_adjlist(g0_sorted);
    rrsGraph::set_adjlist(g1_sorted);
    //EqClass = new ui[g0_sorted.n];
    rrsGraph::GetEqClass(g0_sorted,EqClass);
}


static void worker_thread(int id) {
    while(1){
        rrs::TaskPackage pkg;
        {
            /* 阻塞等任务 */
            std::unique_lock<std::mutex> lk(g_qmtx);
            g_qcv.wait(lk, [] { return !g_taskQ.empty() || thread_should_stop.load();});
            if (thread_should_stop.load() && g_taskQ.empty()) break;
            pkg = std::move(g_taskQ.front());
            g_taskQ.pop();
        }  
        /* ---- 跑 mcs ---- */
        vector<rrs::VtxPair> solution = mcs(pkg.g0, pkg.g1, pkg.vertexs, pkg.unvertexSet, id);
        /* ---- 发回结果 ---- */
        AnsAndTime res{static_cast<int>(solution.size()), SolveTime};
        MPI_Send(&res, 1, MPI_ANS_AND_TIME, pkg.masterid, TAG_NEW_ANS, MPI_COMM_WORLD);
        std::vector<int> vs;
        vs.reserve(res.ans * 2);
        // for (auto &[u, w] : solution) { vs.push_back(vv0[u]); vs.push_back(vv1[w]); }
        for (auto &[u, w] : solution) { vs.push_back(u); vs.push_back(w); }
        MPI_Send(vs.data(), res.ans*2, MPI_INT,pkg.masterid, TAG_RESULT, MPI_COMM_WORLD); 
        // if(thread_should_stop.load()) break;
    }
    // cout<<" debug over "<<id<<endl;
    return;
}

static std::pair<std::vector<int>, std::set<int>> RecvSets() {
    int header[2] = {0, 0};
    MPI_Recv(header, 2, MPI_INT, 0, TAG_HEADER, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int n = header[0];
    int m = header[1];

    std::vector<int> buf(n + m);
    MPI_Recv(buf.data(), n + m, MPI_INT, 0, TAG_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return {
        {buf.begin(), buf.begin() + n}, 
        {buf.begin() + n, buf.end()}      
    };
}
void rrs::run(Worker& worker) { 
    thread_should_stop = false;
    thread wthr(worker_thread, worker.id_);
    MPI_Status status;
    int msg;
    MPI_Irecv(&msg,1,MPI_INT,worker.masterId_,TAG_CMD,MPI_COMM_WORLD,&worker.msgHandlers_.recvRequest);
    while (1) {
        int flag = 0;
        MPI_Test(&worker.msgHandlers_.recvRequest, &flag, &status);
        if (flag) {
            int current_msg = msg;
            if (current_msg == msg::SendInformation) {
                mcs_should_stop = false;
                nodes = 0;
                cut_branches = 0;
                auto [vertexs, unvertexSet] = RecvSets();
                TaskPackage pkg{ g0_sorted,g1_sorted,vertexs,unvertexSet,worker.masterId_ };
                {
                    std::lock_guard<std::mutex> lk(g_qmtx);
                    g_taskQ.push(std::move(pkg));
                }
                g_qcv.notify_all();
            }
            else if (current_msg == msg::UpdateAns) {
                int res;
                MPI_Bcast(&res, 1, MPI_INT, worker.masterId_, MPI_COMM_WORLD);
                ans.store(res, std::memory_order_relaxed);
            }
            else if(current_msg == msg::SendData){
                int res=0;
                MPI_Bcast(&res, 1, MPI_INT, worker.masterId_, MPI_COMM_WORLD);
                if(res){
                    std::lock_guard<std::mutex> lock(data_mtx);
                    std::vector<int> a(res*4);
                    MPI_Bcast(a.data(), res*4, MPI_INT, worker.masterId_, MPI_COMM_WORLD);
                    int lenb = 0;
                    vector<Data> tmp_datas;
                    for(int i = 0; i < res*4; i += 4){
                        Data data = {a[i],a[i+1],a[i+2],a[i+3]};
                        lenb += a[i+3];
                        tmp_datas.push_back(data);
                    }
                    std::vector<int> b(lenb);
                    MPI_Bcast(b.data(), lenb, MPI_INT, worker.masterId_, MPI_COMM_WORLD);
                    lenb=0;
                    for(auto &data:tmp_datas){
                        int len = data.num;
                        while(len--) data.vs.push_back(b[lenb++]);
                        auto it = lower_bound(datas.begin(),datas.end(),data,cmpdata);
                        if(it == datas.end()) {
                            datas.insert(it,data);
                        }
                        else {
                            if(it->v == data.v && it->w == data.w && it->vs == data.vs){
                                if(it->ans < data.ans){
                                    it->ans = data.ans;                         
                                }                           
                            }
                            else {
                                datas.insert(it,data);
                            }
                        } 
                    }                  
                }
            }
            else if (current_msg == msg::UpdateTime) {
                {                  
                    std::lock_guard<std::mutex> lock(mtx); // 加锁
                    int num = new_datas.size();
                    MPI_Send(&num,1,MPI_INT,worker.masterId_,TAG_RECV_DATA1,MPI_COMM_WORLD);
                    if(num){
                        std::vector<int> a(num*4);
                        std::vector<int> b;
                        int lena=0,lenb=0;
                        for(auto &u:new_datas){
                            lenb += u.vs.size();
                        }
                        b.resize(lenb);
                        lena = 0;
                        lenb = 0;
                        for(auto &u:new_datas){
                            a[lena++]=u.v;
                            a[lena++]=u.w;
                            a[lena++]=u.ans;
                            a[lena++]=u.num;
                            for(auto v:u.vs) b[lenb++]=v;
                        }
                        MPI_Send(a.data(), num*4, MPI_INT, worker.masterId_, TAG_RECV_DATA2, MPI_COMM_WORLD);
                        MPI_Send(b.data(), lenb, MPI_INT, worker.masterId_, TAG_RECV_DATA3, MPI_COMM_WORLD);
                        new_datas.clear();
                    }
                }

                AnsAndTime res = { now_ans ,SolveTime };
                MPI_Send(&res, 1, MPI_ANS_AND_TIME, worker.masterId_, TAG_UPDATE_ANS, MPI_COMM_WORLD);
            }
            else if (current_msg == msg::UpdateStatue) {
                mcs_should_stop.store(true, std::memory_order_release);
                g_qcv.notify_all();
            }
            else if (current_msg == msg::Endall) {
                mcs_should_stop.store(true, std::memory_order_release);
                thread_should_stop.store(true, std::memory_order_release);
                g_qcv.notify_all();
                // if (wthr.joinable()) {
                //     g_qcv.notify_all();
                //     cout<<"debug end "<<worker.id_<<endl;
                //     wthr.join();
                // }
                break;
            }
            MPI_Irecv(&msg,1,MPI_INT,worker.masterId_,TAG_CMD,MPI_COMM_WORLD,&worker.msgHandlers_.recvRequest);
        }
        else sleep(1);
    }
    if (wthr.joinable()) {
        g_qcv.notify_all();
        wthr.join();
    }
    // MPI_Barrier(MPI_COMM_WORLD);
    return;
}
