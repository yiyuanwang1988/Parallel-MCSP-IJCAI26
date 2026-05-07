#include"mcsp.hpp"
static double SolveTime = 0;
static int now_ans = 0;

static unsigned long long nodes = 0;
static unsigned long long cut_branches = 0;
static unsigned long long best_count = 0;
unsigned long long calls_for_optimal = 0;



static std::queue<dsb::TaskPackage> g_taskQ;
static std::mutex g_qmtx,mtx,data_mtx;
static std::condition_variable g_qcv;

void show(const vector<dsb::VtxPair>& current, const vector<dsb::Bidomain> &domains,
        const vector<int>& left, const vector<int>& right)
{
    cout << "Nodes: " << nodes << std::endl;
    cout << "Length of current assignment: " << current.size() << std::endl;
    cout << "Current assignment:";
    for (unsigned int i=0; i<current.size(); i++) {
        cout << "  (" << current[i].v << " -> " << current[i].w << ")";
    }
    cout << std::endl;
    for (unsigned int i=0; i<domains.size(); i++) {
        struct dsb::Bidomain bd = domains[i];
        cout << "Left  ";
        for (int j=0; j<bd.left_len; j++)
            cout << left[bd.l + j] << " ";
        cout << std::endl;
        cout << "Right  ";
        for (int j=0; j<bd.right_len; j++)
            cout << right[bd.r + j] << " ";
        cout << std::endl;
    }
    cout << "\n" << std::endl;
}

bool check_sol(const dsbGraph::Graph & g0, const dsbGraph::Graph & g1 , const vector<dsb::VtxPair> & solution) {
    //return true;
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (unsigned int i = 0; i < solution.size(); i++) {
        struct dsb::VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w])
            return false;
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.label[p0.v] != g1.label[p0.w])
            return false;
        for (unsigned int j = i + 1; j < solution.size(); j++) {
            struct dsb::VtxPair p1 = solution[j];
            if (g0.adjmat[p0.v][p1.v] != g1.adjmat[p0.w][p1.w])
                return false;
        }
    }
    return true;
}

const int THRESHOLD = 16;//只给 ≤16 的域做精细化剪枝；再大就用传统 min
const int boxes = 2; //把域按大小差分成两档：≤2 一档，>2 一档。
const int index_partition_size = 2;//每档再按大小差离散成桶，桶宽=2。
const int best_partition_ub = THRESHOLD;//数组下界保护，防止越界。
const int best_partition_diff = 10; //只处理左右大小差 ≤10的域；差再大就不费劲了
int MAX_COUNT = 100;//采样窗口：前 100 次调用做“试用期”。
int MAX_SUCCESS = 30;//如果 100 次里成功≤30 次，认为“费劲不讨好”，关闭精细化剪枝
int called_count = 0;//精细化剪枝被调用的总次数。
int success_count = 0;//其中真正剪掉分支的次数
bool bound_enabled = true;//是否还继续用“费劲”的精细化计算。
bool count_enabled = true;//是否还在统计期

const int ROWS = THRESHOLD - index_partition_size + 2;
//离线训练得到的“期望最大匹配缩减值”表。下标由【域大小】【大小差】二维索引，计算上界时直接查表扣减。
long double arrays[ROWS*boxes];
//临时缓冲区，给精细化计算时统计小端顶点在域内诱导子图的度序列，避免每次都 new vector
int s_degrees[best_partition_ub];

void initializeArray() {
    for (auto & i : arrays) {
        i = 1.0;
    }
}

unsigned int calc_bound(const dsbGraph::Graph &g0, const dsbGraph::Graph &g1, vector<dsb::Bidomain> &domains,
                        const vector<int> &left, const vector<int> &right, const vector<dsb::VtxPair> &incumbent,
                        const vector<dsb::VtxPair> &current) {
    unsigned int incumbent_size = incumbent.size();
    unsigned int current_size = current.size();
    unsigned int mcsplit_bound = current_size;
    long double possible_reductions = 0;
    if (!bound_enabled) {
        for (const dsb::Bidomain &bd: domains) {
            mcsplit_bound += std::min(bd.left_len, bd.right_len);
        }
        return mcsplit_bound;
    }
    for (const dsb::Bidomain &bd: domains) {
        const auto [min_len, max_len] = std::minmax(bd.left_len, bd.right_len);
        const int size_diff = std::abs(bd.left_len - bd.right_len);
        if (min_len > 2 && max_len <= best_partition_ub && size_diff <= best_partition_diff) {
            if (bd.is_valid) {
                //这个域 已经算过
                possible_reductions += min_len - bd.size;
            } else {
                //用 离线统计好的平均 δ 值 立即估算“这个 bidomain 还能再剪多少匹配对”，而 不用实时跑
                possible_reductions += arrays[((max_len - index_partition_size)*index_partition_size) + (size_diff > 2 ? 1 : 0)];
            }
        }
        mcsplit_bound += min_len;
    }
    if (mcsplit_bound <= ans || mcsplit_bound <= incumbent_size || (possible_reductions < 0.5) || ((mcsplit_bound - possible_reductions) > incumbent_size)) {
        return mcsplit_bound;
    }
    if (count_enabled && called_count == MAX_COUNT) {
        count_enabled = false;
        if (success_count < MAX_SUCCESS) {
            bound_enabled = false;
        }
        return mcsplit_bound;
    }
    unsigned int bound = current_size;
    for (dsb::Bidomain &bd: domains) {
        const auto [min_len, max_len] = std::minmax(bd.left_len, bd.right_len);
        const int size_diff = std::abs(bd.left_len - bd.right_len);
        if (bd.is_valid) {
            bound += bd.size;
            continue;
        }
        if (bd.left_len == 2 && bd.right_len == 2) {
            if (g0.adjmat[left[bd.l]][left[bd.l+1]] != g1.adjmat[right[bd.r]][right[bd.r+1]]) {
                bound += 1;
                bd.size = 1;
                int index = (max_len - index_partition_size)*index_partition_size;
                arrays[index] = (arrays[index] + 1) / 2;
            } else {
                bd.size = 2;
                int index = (max_len - index_partition_size)*index_partition_size;
                arrays[index] = (arrays[index]) / 2;
                bound += 2;
            }
            bd.is_valid = true;
            continue;
        }
        if (min_len > 2 && max_len <= best_partition_ub && size_diff <= best_partition_diff) {
            int s_len, b_len;
            int val = min_len;
            int index = (max_len - index_partition_size)*index_partition_size + (size_diff > 2 ? 1 : 0);
            int s_connections = 0, b_connections = 0; // small and big side
//            std::vector<int> s_degrees(min_len);
            std::fill(s_degrees, s_degrees + min_len, 0);//小端度的序列数
            if (bd.left_len <= bd.right_len) {
                s_len = bd.left_len;
                b_len = bd.right_len;
                for (int i = 0; i < s_len; ++i) {
                    const int v1 = left[bd.l +i]; // using references can only make things more expensive for small data types
                    for (int j = i + 1; j < s_len; ++j) {
                        const int v2 = left[bd.l + j];
                        if (g0.adjmat[v1][v2]) {
                            s_degrees[i]++;
                            s_degrees[j]++;
                            s_connections++; // minor trick: we only count undirected edges (everything divided by 2)
                        }
                    }
                }
                // no need to compute degrees for the bigger side
                for (int i = 0; i < b_len; ++i) {
                    const int u1 = right[bd.r + i];
                    for (int j = i + 1; j < b_len; ++j) {
                        const int u2 = right[bd.r + j];
                        if (g1.adjmat[u1][u2])
                            b_connections++;
                    }
                }
            } else {
                s_len = bd.right_len;
                b_len = bd.left_len;
                for (int i = 0; i < s_len; ++i) {
                    const int v1 = right[bd.r + i]; // using references can only make things more expensive for small data types
                    for (int j = i + 1; j < s_len; ++j) {
                        const int v2 = right[bd.r + j];
                        if (g1.adjmat[v1][v2]) {
                            s_degrees[i]++;
                            s_degrees[j]++;
                            s_connections++; // minor trick: we only count undirected edges (everything divided by 2)
                        }
                    }
                }
                // no need to compute degrees for the bigger side
                for (int i = 0; i < b_len; ++i) {
                    const int u1 = left[bd.l + i];
                    for (int j = i + 1; j < b_len; ++j) {
                        const int u2 = left[bd.l + j];
                        if (g0.adjmat[u1][u2])
                            b_connections++;
                    }
                }
            }
            // compute non-connections from connections
            int s_non_connections = s_len * (s_len - 1) / 2 - s_connections;
            int b_non_connections = b_len * (b_len - 1) / 2 - b_connections;
            // optimization starts here
            int diff = s_connections - b_connections;
            if (diff > 0) {
//                std::sort(s_degrees.begin(), s_degrees.end());
                std::sort(s_degrees, s_degrees+s_len);
                diff -= s_degrees[--val];
                while (diff > 0)
                    diff -= s_degrees[--val];
            } else {
                diff = s_non_connections - b_non_connections;
                if (diff > 0) {
//                    std::sort(s_degrees.begin(), s_degrees.end(), std::greater<>());
                    std::sort(s_degrees, s_degrees+s_len, std::greater<>());
                    diff -= s_len - s_degrees[--val];
                    while (diff > 0)
                        diff -= s_len - s_degrees[--val];
                }
            }
            arrays[index] = (arrays[index] + min_len - val) / 2;
            bound += val;
            bd.size = val;
        } else {
            bound += min_len;
            bd.size = min_len;
        }
        bd.is_valid = true;
    }
    if (count_enabled) {
        if (bound <= incumbent_size) {
            success_count++;
        }
        called_count++;
    }
    return bound;
}

int find_min_value(const vector<int>& arr, int start_idx, int len) {
    int min_v = INT_MAX;
    for (int i=0; i<len; i++)
        if (arr[start_idx + i] < min_v)
            min_v = arr[start_idx + i];
    return min_v;
}

int select_bidomain(const vector<dsb::Bidomain>& domains, const vector<int> & left,
        int current_matching_size)
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int best = -1;
    for (unsigned int i=0; i<domains.size(); i++) {
        const dsb::Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size>0 && !bd.is_adjacent) continue;
        int len = arguments.heuristic == min_max ?
                std::max(bd.left_len, bd.right_len) :
                bd.left_len * bd.right_len;
        if (len < min_size) {
            min_size = len;
            min_tie_breaker = find_min_value(left, bd.l, bd.left_len);
            best = i;
        } else if (len == min_size) {
            int tie_breaker = find_min_value(left, bd.l, bd.left_len);
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}

// Returns length of left half of arrays
static int partition(vector<int>& all_vv, int start, int len, const vector<unsigned int> & adjrow) {
    int i=0;
    for (int j=0; j<len; j++) {
        if (adjrow[all_vv[start+j]]) {
            std::swap(all_vv[start+i], all_vv[start+j]);
            i++;
        }
    }
    return i;
}

// multiway is for directed and/or labelled graphs
vector<dsb::Bidomain> filter_domains(const vector<dsb::Bidomain> & d, vector<int> & left,
        vector<int> & right, const dsbGraph::Graph & g0, const dsbGraph::Graph & g1, int v, int w,
        bool multiway)
{
    vector<dsb::Bidomain> new_d;
    new_d.reserve(d.size());
    for (const dsb::Bidomain &old_bd : d) {
        int l = old_bd.l;
        int r = old_bd.r;
        // After these two partitions, left_len and right_len are the lengths of the
        // arrayss of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, old_bd.left_len, g0.adjmat[v]);
        int right_len = partition(right, r, old_bd.right_len, g1.adjmat[w]);
        int left_len_noedge = old_bd.left_len - left_len;
        int right_len_noedge = old_bd.right_len - right_len;
        if (left_len_noedge && right_len_noedge) {
            if (left_len != 0 || right_len != 0) {
                new_d.emplace_back(l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent, -1,
                                   false);
            } else {
                new_d.emplace_back(l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent, old_bd.size, old_bd.is_valid);
            }
        }

        if (multiway && left_len && right_len) {
            auto& adjrow_v = g0.adjmat[v];
            auto& adjrow_w = g1.adjmat[w];
            auto l_begin = std::begin(left) + l;
            auto r_begin = std::begin(right) + r;
            std::sort(l_begin, l_begin+left_len, [&](int a, int b)
                    { return adjrow_v[a] < adjrow_v[b]; });
            std::sort(r_begin, r_begin+right_len, [&](int a, int b)
                    { return adjrow_w[a] < adjrow_w[b]; });
            int l_top = l + left_len;
            int r_top = r + right_len;
            while (l<l_top && r<r_top) {
                unsigned int left_label = adjrow_v[left[l]];
                unsigned int right_label = adjrow_w[right[r]];
                if (left_label < right_label) {
                    l++;
                } else if (left_label > right_label) {
                    r++;
                } else {
                    int lmin = l;
                    int rmin = r;
                    do { l++; } while (l<l_top && adjrow_v[left[l]]==left_label);
                    do { r++; } while (r<r_top && adjrow_w[right[r]]==left_label);
                    new_d.emplace_back(lmin, rmin, l-lmin, r-rmin, true, -1, false);
                }
            }
        } else if (left_len && right_len) {
            if (left_len_noedge != 0 || right_len_noedge != 0) {
                new_d.emplace_back(l, r, left_len, right_len, true, -1, false);
            } else {
                new_d.emplace_back(l, r, left_len, right_len, true, old_bd.size, old_bd.is_valid);
            }
        }
    }
    return new_d;
}

// returns the index of the smallest value in arr that is >w.
// Assumption: such a value exists
// Assumption: arr contains no duplicates
// Assumption: arr has no values==INT_MAX
int index_of_next_smallest(const vector<int>& arr, int start_idx, int len, int w,int v) {
    int idx = -1;
    int smallest = INT_MAX;
    for (int i=0; i<len; i++) {
        int vtx = arr[start_idx+i];
        if (vtx>w && vtx<smallest) 
        {
            smallest = vtx;
            idx = i;
        }

    }
    return idx;
}

void remove_vtx_from_left_domain(vector<int>& left, dsb::Bidomain& bd, int v)
{
    int i = 0;
    while(left[bd.l + i] != v) i++;
    std::swap(left[bd.l+i], left[bd.l+bd.left_len-1]);
    bd.left_len--;
}

void remove_bidomain(vector<dsb::Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size()-1];
    domains.pop_back();
}

static int solve(const dsbGraph::Graph& g0, const dsbGraph::Graph& g1, vector<dsb::VtxPair>& incumbent,
    vector<dsb::VtxPair>& current, vector<dsb::Bidomain>& domains, vector<int> vertexs,vector<int> unvertexs,vector<int> data_unvertexs,
    vector<int>& left, vector<int>& right, unsigned int matching_size_goal,const int id,int& current_best)
{
    //if (GetSolveTime() > cutoff_time) return -1;
    if (mcs_should_stop.load()) {
        return -1;
    }

    if (arguments.verbose) show(current, domains, left, right);
    if (current.size() > incumbent.size()) {
        now_ans = current.size();
        incumbent = current;
        SolveTime = GetSolveTime();
        calls_for_optimal = nodes;
        if (!arguments.quiet) cout << "Incumbent size: " << incumbent.size() << endl;
    }
    if(current.size() > current_best){
        current_best = current.size();
    }
    nodes++;
    //理论最大可能匹配数
    unsigned int bound = calc_bound(g0, g1, domains, left, right, incumbent, current);
    if (bound <= getAns() || bound <= incumbent.size() || bound < matching_size_goal) {
        if(current_best < bound) current_best = bound;
        cut_branches++;
        return -1;
    }
    // if (arguments.big_first && incumbent.size() == matching_size_goal) return -1;

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
    int v;
    dsb::Bidomain* bd = nullptr;
    int bd_idx = -1;
    if (vertexs.size()) {
        v = vertexs.front();
        vertexs.erase(vertexs.begin());
        int flag = 0;
        for (int i = 0; i < domains.size(); i++) {
            for (int j = 0; j < domains[i].left_len; j++) {
                if (left[domains[i].l + j] == v) {
                    bd = &domains[i];
                    bd_idx = i;
                    flag = 1;
                    break;
                }
            }
            if (flag) break;
        }
        if(!flag) {
            return -1;
        }
    }
    else {
        //在所有剩余双向域里挑一个 “最划算” 的域返回下标；
        bd_idx = select_bidomain(domains, left, current.size());
        if (bd_idx == -1) { // In the MCCS case, there may be nothing we can branch on
            return -1;
        }
        bd = &domains[bd_idx];

        //找最小顶点编号（即 min_vertex），作为当前要尝试的变量。
        v = find_min_value(left, bd->l, bd->left_len);
        //把刚才选到的 v 从该双向域的左端抽掉（与区间最后一个元素交换，left_len--），保证后续循环不会重复选它。        
    }
    remove_vtx_from_left_domain(left, domains[bd_idx], v);
    int w = -1;
    bd->right_len--;



    int len = bd->right_len;
    for (int i = len; i >= 0; i --) {
        //找大于当前 w 且编号最小的那个顶点，返回相对下标 idx
        int idx = index_of_next_smallest(right, bd->r+bd->right_len-len, len + 1, w,v);
        if(idx==-1) break;
        //拿到真正的顶点编号 w，准备与左端 v 配对
        w = right[bd->r + idx];
        // swap w to the end of its colour class
        right[bd->r + idx] = right[bd->r + bd->right_len];
        right[bd->r + bd->right_len] = w;

        //把当前匹配对 (v,w) 嵌入到所有双向域
        auto new_domains = filter_domains(domains, left, right, g0, g1, v, w,
            arguments.directed || arguments.edge_labelled);
        current.push_back(dsb::VtxPair(v, w));
        int t = solve(g0, g1, incumbent, current, new_domains, vertexs, unvertexs, data_unvertexs,left, right, matching_size_goal,id,current_best);
        current.pop_back();
        
        if(t!=-1){
            if(t!=v) {
                return t;
            }
        }
        //if(GetSolveTime() > cutoff_time) return -1;
        if(current.empty()&&unvertexs.size() < 10){
            if(current_best < bound) current_best = bound;
            Data data = {v,w,current_best,unvertexs.size(),unvertexs};
            {
                std::lock_guard<std::mutex> lock(mtx); // 加锁
                new_datas.push_back(data);
            }
            current_best = 0;           
        }        
    }
    bd->right_len++;
    if (bd->left_len == 0)
        remove_bidomain(domains, bd_idx);

    unvertexs.insert(lower_bound(unvertexs.begin(),unvertexs.end(),v),v);
    int t = solve(g0, g1, incumbent, current, domains, vertexs, unvertexs, data_unvertexs, left, right, matching_size_goal,id,current_best);
    return t;
}



static void erase_g0_vertices(dsbGraph::Graph& g0,
                       std::vector<dsb::Bidomain>& domains,
                       std::vector<int>& left,
                       const std::set<int>& uv)
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

    /* Rebuild `left` and `domains` safely from the original `left` snapshot.
       The previous code erased elements from `left` first and then used
       the old domain indices (bd.l, bd.left_len) against the shortened
       vector, which leads to out-of-bounds accesses and heap corruption.
    */
    vector<int> old_left = left; // preserve original ordering and indices
    vector<dsb::Bidomain> new_domains;
    vector<int> new_left;
    new_left.reserve(old_left.size());

    int pos = 0; // write position in new_left
    for (const dsb::Bidomain &bd : domains) {
        int l_old = bd.l;
        int left_len = bd.left_len;

        /* count surviving vertices and append them into new_left in one pass */
        int cnt = 0;
        for (int i = 0; i < left_len; ++i) {
            int vv = old_left[l_old + i];
            if (!uv.count(vv)) {
                new_left.push_back(vv);
                ++cnt;
            }
        }
        if (cnt == 0) continue; // drop empty domain

        new_domains.emplace_back(pos, bd.r, cnt, bd.right_len, bd.is_adjacent, -1, false);
        pos += cnt;
    }

    // swap in the rebuilt arrays
    left.swap(new_left);
    domains.swap(new_domains);

    /* 4. 更新Graph */
    dsbGraph::Graph new_g0(old_n);
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
    g0 = std::move(new_g0);
}

static vector<dsb::VtxPair> mcs(dsbGraph::Graph & g0, const dsbGraph::Graph & g1, std::vector<int> v, std::set<int> uv,const int id) {
    
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions

    auto domains = vector<dsb::Bidomain> {};

    std::set<unsigned int> left_labels;
    std::set<unsigned int> right_labels;
    for (unsigned int label : g0.label) left_labels.insert(label);
    for (unsigned int label : g1.label) right_labels.insert(label);
    std::set<unsigned int> labels;  // labels that appear in both graphs
    std::set_intersection(std::begin(left_labels),
                          std::end(left_labels),
                          std::begin(right_labels),
                          std::end(right_labels),
                          std::inserter(labels, std::begin(labels)));

    // Create a bidomain for each label that appears in both graphs
    for (unsigned int label : labels) {
        int start_l = left.size();
        int start_r = right.size();
        for (int i=0; i<g0.n; i++)
            if (g0.label[i]==label){
                left.push_back(i);
            }     
        for (int i=0; i<g1.n; i++)
            if (g1.label[i]==label){
                right.push_back(i);
            }
        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false, -1, false});
    }
    erase_g0_vertices(g0, domains, left, uv);
        

    vector<dsb::VtxPair> incumbent;
    vector<int> uvs(uv.begin(),uv.end());
    sort(uvs.begin(),uvs.end());
    initializeArray();
    nodes = 0;
    cut_branches = 0;
    int current_best = 0;
    if (arguments.big_first) {
        for (int k=0; k<g0.n; k++) {
            unsigned int goal = g0.n - k;
            auto left_copy = left;
            auto right_copy = right;
            auto domains_copy = domains;
            vector<dsb::VtxPair> current;
            solve(g0, g1, incumbent, current, domains_copy, v, uvs,uvs,left_copy, right_copy, goal,id,current_best);
            if (incumbent.size() == goal) break;
            if (!arguments.quiet) cout << "Upper bound: " << goal-1 << std::endl;
        }

    } else {
        vector<dsb::VtxPair> current;
        solve(g0, g1, incumbent, current, domains, v, uvs,uvs,left, right, 1,id, current_best);
    }
    //std::cout << "id " << id << "debug dsb   " <<  check_sol(g0, g1, incumbent) << " incumbent " << incumbent.size() << " time " << SolveTime << " totaltime " << GetSolveTime() << " " << thread_should_stop.load() << "  " << mcs_should_stop.load() << std::endl;
    return incumbent;
}

inline int sum(const std::vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}
inline std::vector<int> calculate_degrees(const dsbGraph::Graph & g) {
    std::vector<int> degree(g.n, 0);
    for (int v=0; v<g.n; v++) {
        for (int w=0; w<g.n; w++) {
            unsigned int mask = 0xFFFFu;
            if (g.adjmat[v][w] & mask) degree[v]++;
            if (g.adjmat[v][w] & ~mask) degree[v]++;  // inward edge, in directed case
        }
    }
    return degree;
}
void dsb::init(){
    initializeArray();
    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    struct dsbGraph::Graph g0 = dsbGraph::readGraph(arguments.filename1, format, arguments.directed,
                            arguments.edge_labelled, arguments.vertex_labelled);
    struct dsbGraph::Graph g1 = dsbGraph::readGraph(arguments.filename2, format, arguments.directed,
                            arguments.edge_labelled, arguments.vertex_labelled);

    std::vector<int> g0_deg = calculate_degrees(g0);
    std::vector<int> g1_deg = calculate_degrees(g1);
    
    vv0.resize(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    bool g1_dense = sum(g1_deg) > g1.n*(g1.n-1);//稠密图 
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
        return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
    });

    vv1.resize(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    bool g0_dense = sum(g0_deg) > g0.n*(g0.n-1);
    std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {
        return g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
    });

    g0_sorted = dsbGraph::induced_subgraph(g0, vv0);
    g1_sorted = dsbGraph::induced_subgraph(g1, vv1);
}


static void worker_thread(int id) {
    while(1){
        dsb::TaskPackage pkg;
        {
            /* 阻塞等任务 */
            std::unique_lock<std::mutex> lk(g_qmtx);
            g_qcv.wait(lk, [] { return !g_taskQ.empty() || thread_should_stop.load();});
            if (thread_should_stop.load() && g_taskQ.empty()) break;
            pkg = std::move(g_taskQ.front());
            g_taskQ.pop();
        }     
        /* ---- 跑 mcs ---- */
        vector<dsb::VtxPair> solution = mcs(pkg.g0, pkg.g1, pkg.vertexs, pkg.unvertexSet,id);
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
void dsb::run(Worker& worker) {
    thread_should_stop = false;
    thread wthr(worker_thread, worker.id_);

    int msg;
    MPI_Irecv(&msg, 1, MPI_INT, worker.masterId_, TAG_CMD, MPI_COMM_WORLD, &worker.msgHandlers_.recvRequest);
    MPI_Status status;
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
                //cout << "id " << worker.id_ << " UpdateStatue mcs_should_stop " << mcs_should_stop.load()  << endl;
                g_qcv.notify_all();
            }
            else if (current_msg == msg::Endall) {
                mcs_should_stop.store(true, std::memory_order_release);
                thread_should_stop.store(true, std::memory_order_release);
                //cout << "id " << worker.id_ << " UpdateStatue mcs_should_stop " << mcs_should_stop.load()  << " thread_should_stop " << thread_should_stop.load() << endl;
                g_qcv.notify_all();
                // if (wthr.joinable()) {
                //     g_qcv.notify_all();
                //     cout<<"debug end "<<worker.id_<<endl;
                //     wthr.join();
                // }
                break;
            }
            MPI_Irecv(&msg, 1, MPI_INT, worker.masterId_, TAG_CMD, MPI_COMM_WORLD, &worker.msgHandlers_.recvRequest);
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
