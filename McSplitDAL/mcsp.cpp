#include"mcsp.hpp"

#define VSCORE
//#define DEBUG
#define Best
using std::vector;
using std::cout;
using std::endl;

static double SolveTime = 0;
static int now_ans = 0;

static unsigned long long nodes{ 0 };
static unsigned long long cutbranches{0};
static unsigned long long conflicts=0;
static unsigned long long bestnodes=0,bestcount=0;
static int dl=0;
static int Maxnum=0;
static int M=1;
static int num = 0;
const int short_memory_threshold = 1e5;
const int long_memory_threshold = 1e9;
// static vector<int> scores[2000];
static std::queue<dal::TaskPackage> g_taskQ;
static std::mutex g_qmtx,mtx,data_mtx;
static std::condition_variable g_qcv;
void show(const vector<dal::VtxPair>& current, const vector<dal::Bidomain> &domains,
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
        struct dal::Bidomain bd = domains[i];
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

bool check_sol(const dalGraph::Graph & g0, const dalGraph::Graph & g1 , const vector<dal::VtxPair> & solution) {
    // return true;
    vector<bool> used_left(g0.n, false);
    vector<bool> used_right(g1.n, false);
    for (unsigned int i=0; i<solution.size(); i++) {
        struct dal::VtxPair p0 = solution[i];
        if (used_left[p0.v] || used_right[p0.w]){
            cout<<"Vertex repeated in the solution "<<p0.v<<" "<<used_left[p0.v]<<" "<<p0.w<<" "<<used_right[p0.w]<<endl;
            return false;
        }
            
        used_left[p0.v] = true;
        used_right[p0.w] = true;
        if (g0.label[p0.v] != g1.label[p0.w]){
            cout<<"Label mismatch"<<endl;
            return false;
        }   
        for (unsigned int j=i+1; j<solution.size(); j++) {
            struct dal::VtxPair p1 = solution[j];
            if (g0.adjmat[p0.v][p1.v] != g1.adjmat[p0.w][p1.w]){
                cout<<"Topological mismatch"<<endl;
                return false;
            }     
        }
    }
    return true;
}

int calc_bound(const vector<dal::Bidomain>& domains) {
    int bound = 0;
    for (const dal::Bidomain &bd : domains) {
        bound += std::min(bd.left_len, bd.right_len);
    }
    return bound;
}
int selectV_index(const vector<int> &arr, const vector<gtype> &lgrade, int start_idx, int len) {
    int idx = -1;
    gtype max_g = -1;
    int vtx, best_vtx = INT_MAX;
    for(int i = 0; i < len; i++) {
        vtx = arr[start_idx + i];
        if(lgrade[vtx] > max_g) {
            idx = i;
            best_vtx = vtx;
            max_g = lgrade[vtx];
        }
        else if(lgrade[vtx] == max_g) {
            if(vtx < best_vtx) {
                idx = i;
                best_vtx = vtx;
            }
        }
    }
    return idx;
}

int select_bidomain(const vector<dal::Bidomain>& domains, const vector<int> & left, const vector<gtype> &lgrade,
        int current_matching_size)
{
    // Select the bidomain with the smallest max(leftsize, rightsize), breaking
    // ties on the smallest vertex index in the left set
    int min_size = INT_MAX;
    int min_tie_breaker = INT_MAX;
    int tie_breaker;
    unsigned int i;  int len;
    int best = -1;
    for (i=0; i<domains.size(); i++) {
        const dal::Bidomain &bd = domains[i];
        if (arguments.connected && current_matching_size>0 && !bd.is_adjacent) continue;
            len = arguments.heuristic == min_max ?
                std::max(bd.left_len, bd.right_len) :
                bd.left_len * bd.right_len;//最多有这么多种连线
        if (len < min_size) {//优先计算可能连线少的情况,当domain中可能配对数最小时，选择domain中含有最大度的点那个domain
            min_size = len;//此时进行修改将度最大的点进行优化
            min_tie_breaker = left[bd.l + selectV_index(left,lgrade,bd.l, bd.left_len)];//序号越小度越大
            best = i;
        } else if (len == min_size) {
            tie_breaker = left[bd.l + selectV_index(left,lgrade,bd.l, bd.left_len)];
            if (tie_breaker < min_tie_breaker) {
                min_tie_breaker = tie_breaker;
                best = i;
            }
        }
    }
    return best;
}

// Returns length of left half of array
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

int remove_matched_vertex(vector<int> &arr, int start, int len, const vector<int> &matched) {
    int p = 0;
    for(int i = 0; i < len; i++) {
        if(!matched[arr[start + i]]) {
            std::swap(arr[start + i], arr[start + p]);
            p++;
        }
    }
    return p;
}
vector<dal::Bidomain> rewardfeed_RL(const vector<dal::Bidomain> & d,int bd_idx,vector<dal::VtxPair> & current, vector<int> & left,
                               vector<int> & right,vector<gtype> & lgrade,vector<gtype> & rgrade,const dalGraph::Graph & g0, const dalGraph::Graph & g1, int v, int w,
                               vector<int> &g0_matched, vector<int> &g1_matched, bool multiway)
    {//每个domain均分成2个new_domain分别是与当前对应点相连或者不相连
       vector<dal::Bidomain> new_d;
       new_d.reserve(d.size());
       current.push_back(dal::VtxPair(v,w));
       g0_matched[v] = 1;
       g1_matched[w] = 1;
       //unsigned int old_bound = current.size() + calc_bound(d)+1;//这里的domain是已经去掉选了的点，但是该点还没纳入current所以+1
       //unsigned int new_bound(0);
       int l,r,j=-1; unsigned int i;
       int temp=0,total=0;
       for (const dal::Bidomain &old_bd : d) {
           j++;
           l = old_bd.l;
           r = old_bd.r;
           // Remove already matched vertices from consideration, then partition
           int unmatched_left_len = remove_matched_vertex(left, l, old_bd.left_len, g0_matched);
           int unmatched_right_len = remove_matched_vertex(right, r, old_bd.right_len, g1_matched);
           // After these two partitions, left_len and right_len are the lengths of the
           // arrays of vertices with edges from v or w (int the directed case, edges
           // either from or to v or w)
           int left_len = partition(left, l, unmatched_left_len, g0.adjmat[v]);//将与选出来的顶点相连顶点数返回，
           int right_len = partition(right, r, unmatched_right_len, g1.adjmat[w]);//并且将相连顶点依次交换到数组前面
           int left_len_noedge = old_bd.left_len - left_len;
           int right_len_noedge = old_bd.right_len - right_len;
    //这里传递v,w选取的下标到bd_idx，j用来计算当前所分domain的下标，这样就能将bound更精确地计算
           temp=std::min(old_bd.left_len,old_bd.right_len)-std::min(left_len,right_len)-std::min(left_len_noedge,right_len_noedge);
           total+=temp;

    #ifdef DEBUG

           printf("adj=%d ,noadj=%d ,old=%d ,temp=%d \n",std::min(left_len,right_len),
                  std::min(left_len_noedge,right_len_noedge),std::min(old_bd.left_len,old_bd.right_len),temp);
           cout<<"j="<<j<< "  idx="<<bd_idx<<endl;
           cout<<"gl="<<lgrade[v]<<" gr="<<rgrade[w]<<endl;
    #endif
           if (left_len_noedge && right_len_noedge)//new_domain存在的条件是同一domain内需要存在与该对应点同时相连，或者同时不相连的点的点
               new_d.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
               //Question:when the value of 'left_len_noedge' will be 0,
               //what does it mean if it happens?
           if (multiway && left_len && right_len) {//不与顶点相连的顶点
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
                       new_d.push_back({lmin, rmin, l-lmin, r-rmin, true});
                   }
               }
           } else if (left_len && right_len) {
               new_d.push_back({l, r, left_len, right_len, true});//与顶点相连的顶点 标志为Ture
           }
       }
       if(total>0){
           conflicts++;
           lgrade[v]+=total;
           if(lgrade[v]>short_memory_threshold){
               for(i=0;i<lgrade.size();i++)
                   lgrade[i]=lgrade[i]/2;
           }
           rgrade[w]+=total;
           if(rgrade[w]>short_memory_threshold){
               for(i=0;i<rgrade.size();i++)
                   rgrade[i]=rgrade[i]/2;
           }
       }
    #ifdef DEBUG
     cout<<"new domains are "<<endl;
     for (const dal::Bidomain &testd : new_d){
         l = testd.l;
         r = testd.r;
         for(j=0;j<testd.left_len;j++)
             cout<<left[l+j] <<" ";
         cout<<" ; ";
         for(j=0;j<testd.right_len;j++)
             <<right[r+j]<<" " ;
         cout<<endl;
     }
    #endif
       return new_d;
    }





// multiway is for directed and/or labelled graphs
vector<dal::Bidomain> rewardfeed(const vector<dal::Bidomain> & d,int bd_idx,vector<dal::VtxPair> & current, vector<int> &g0_matched, vector<int> &g1_matched,
        vector<int> & left, vector<int> & right,vector<gtype> & lgrade,vector<gtype> & rgrade,
        const dalGraph::Graph & g0, const dalGraph::Graph & g1, int v, int w,
        bool multiway)
{//每个domain均分成2个new_domain分别是与当前对应点相连或者不相连
//    assert(g0_matched[v] == 0);
//    assert(g1_matched[w] == 0);
    current.push_back(dal::VtxPair(v, w));
    g0_matched[v] = 1;
    g1_matched[w] = 1;

    int leaves_match_size = 0, v_leaf, w_leaf;
    for(unsigned int i = 0, j = 0; i < g0.leaves[v].size() && j < g1.leaves[w].size(); ) {
        if(g0.leaves[v][i].first < g1.leaves[w][j].first) i++;
        else if(g0.leaves[v][i].first > g1.leaves[w][j].first) j++;
        else {
            const vector<int> &leaf0 = g0.leaves[v][i].second;
            const vector<int> &leaf1 = g1.leaves[w][j].second;
            for(unsigned int p = 0, q = 0; p < leaf0.size() && q < leaf1.size(); ) {
                if(g0_matched[leaf0[p]]) p++;
                else if(g1_matched[leaf1[q]]) q++;
                else {
                    v_leaf = leaf0[p], w_leaf = leaf1[q];
                    p++, q++;
                    current.push_back(dal::VtxPair(v_leaf, w_leaf));
                    g0_matched[v_leaf] = 1;
                    g1_matched[w_leaf] = 1;
                    leaves_match_size++;
                }
            }
            i++, j++;
        }
    }


    vector<dal::Bidomain> new_d;
    new_d.reserve(d.size());
    //unsigned int old_bound = current.size() + calc_bound(d)+1;//这里的domain是已经去掉选了的点，但是该点还没纳入current所以+1
    //unsigned int new_bound(0);
    int l,r,j=-1;
    int temp=0,total=0;
    int unmatched_left_len, unmatched_right_len;
    for (const dal::Bidomain &old_bd : d) {
        j++;
        l = old_bd.l;
        r = old_bd.r;
        if(leaves_match_size > 0 && old_bd.is_adjacent == false) {
            unmatched_left_len = remove_matched_vertex(left, l, old_bd.left_len, g0_matched);
            unmatched_right_len = remove_matched_vertex(right, r, old_bd.right_len, g1_matched);
        }
        else {
            unmatched_left_len = old_bd.left_len;
            unmatched_right_len = old_bd.right_len;
        }
        // After these two partitions, left_len and right_len are the lengths of the
        // arrays of vertices with edges from v or w (int the directed case, edges
        // either from or to v or w)
        int left_len = partition(left, l, unmatched_left_len, g0.adjmat[v]);//将与选出来的顶点相连顶点数返回，
        int right_len = partition(right, r, unmatched_right_len, g1.adjmat[w]);//并且将相连顶点依次交换到数组前面
        int left_len_noedge = unmatched_left_len - left_len;
        int right_len_noedge = unmatched_right_len - right_len;
//这里传递v,w选取的下标到bd_idx，j用来计算当前所分domain的下标，这样就能将bound更精确地计算
        temp=std::min(old_bd.left_len,old_bd.right_len)-std::min(left_len,right_len)-std::min(left_len_noedge,right_len_noedge);
        total+=temp;

#ifdef DEBUG

        printf("adj=%d ,noadj=%d ,old=%d ,temp=%d \n",std::min(left_len,right_len),
               std::min(left_len_noedge,right_len_noedge),std::min(old_bd.left_len,old_bd.right_len),temp);
        cout<<"j="<<j<< "  idx="<<bd_idx<<endl;
        cout<<"gl="<<lgrade[v]<<" gr="<<rgrade[w]<<endl;
 #endif
        if (left_len_noedge && right_len_noedge)//new_domain存在的条件是同一domain内需要存在与该对应点同时相连，或者同时不相连的点的点
            new_d.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge, old_bd.is_adjacent});
        if (multiway && left_len && right_len) {//不与顶点相连的顶点
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
                    new_d.push_back({lmin, rmin, l-lmin, r-rmin, true});
                }
            }
        } else if (left_len && right_len) {
            new_d.push_back({l, r, left_len, right_len, true});//与顶点相连的顶点 标志为Ture
        }
    }
    int domgrade = new_d.size();
    total=total+domgrade;
    if(total>0){
        conflicts++;

        lgrade[v] += total;
        rgrade[w] += total;

        if(lgrade[v] > short_memory_threshold)
            for(int i = 0; i < g0.n; i++)
                lgrade[i] = lgrade[i] / 2;
        if(rgrade[w] > long_memory_threshold)
            for(int i = 0; i < g1.n; i++)
                rgrade[i] = rgrade[i] / 2;
    }
#ifdef DEBUG
  cout<<"new domains are "<<endl;
  for (const dal::Bidomain &testd : new_d){
      l = testd.l;
      r = testd.r;
      for(j=0;j<testd.left_len;j++)
          cout<<left[l+j] <<" ";
      cout<<" ; ";
      for(j=0;j<testd.right_len;j++)
          cout<<right[r+j]<<" " ;
      cout<<endl;
  }
#endif
    return new_d;
}

int selectW_index(const vector<int> &arr, const vector<gtype> &rgrade, int start_idx, int len, const vector<int> &wselected,int v) { 
    int idx = -1;
    gtype max_g = -1;
    int vtx, best_vtx = INT_MAX;
    for(int i = 0; i < len; i++) {
        vtx = arr[start_idx + i];
        if(wselected[vtx] == 0) {
            if(rgrade[vtx] > max_g) {
                idx = i;
                best_vtx = vtx;
                max_g = rgrade[vtx];
            }
            else if(rgrade[vtx] == max_g) {
                if(vtx < best_vtx) {
                    idx = i;
                    best_vtx = vtx;
                }
            }         
        }
    }
    return idx;
}

void remove_vtx_from_array(vector<int> &arr, int start_idx, int &len, int remove_idx) {
    len--;
    std::swap(arr[start_idx + remove_idx], arr[start_idx + len]);
}

void remove_bidomain(vector<dal::Bidomain>& domains, int idx) {
    domains[idx] = domains[domains.size()-1];
    domains.pop_back();
}

static int solve(const dalGraph::Graph & g0, const dalGraph::Graph & g1,
                  vector<gtype> &V, vector<gtype> &lgrade,
                  vector<gtype> &rgrade, vector<vector<gtype>> &Q,
                  vector<dal::VtxPair> & incumbent,
                  vector<dal::VtxPair> & current,
                  vector<int> &g0_matched, vector<int> &g1_matched,
                  vector<dal::Bidomain> & domains,
                  vector<int>& vertexs,        // 优先点集，可 pop()
                  vector<int> unvertexs,
                  vector<int> data_unvertexs,                  
                  vector<int> & left, vector<int> & right,
                  unsigned int matching_size_goal,int& current_best)
{
    //if (GetSolveTime() > cutoff_time) return -1;
    if (mcs_should_stop.load()) return -1;

    nodes++;
    if (current.size() > incumbent.size()) {
        num = 0;
        incumbent = current;
        now_ans = current.size();
        SolveTime = GetSolveTime();
        bestcount = cutbranches + 1;
        bestnodes = nodes;
        if (!arguments.quiet)  cout << "Incumbent size: " << incumbent.size() << endl;
    }
    if(current.size() > current_best){
        current_best = current.size();
    }

    unsigned int bound = current.size() + calc_bound(domains);
    if (bound <= incumbent.size() || bound < matching_size_goal || bound <= getAns()) {
        if(current_best < bound) current_best = bound;
        cutbranches++;
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
    /* ---------- 1. 选 v ---------- */
    int v, w;
    int tmp_idx = -1;
    dal::Bidomain *bd_ptr = nullptr;
    int bd_idx = -1;

    if (!vertexs.empty()) {
        /* 优先从 vertexs 弹出第一个顶点 */
        v = vertexs.front();
        vertexs.erase(vertexs.begin());
        /* 找到它所在的 bidomain */
        for (int i = 0; i < (int)domains.size() && !~tmp_idx; ++i)
            for (int j = 0; j < domains[i].left_len; ++j)
                if (left[domains[i].l + j] == v) {
                    bd_ptr = &domains[i];
                    bd_idx = i;
                    tmp_idx = j;
                    break;
                }
        if (!~tmp_idx) return -1;          // 该顶点已不在任何域，直接回溯
    } else {
        /* 原逻辑：先选域再选最小顶点 */
        num++;
        if (num > Maxnum) {
            M = M % 2 + 1;
            num = 0;
        }
        switch (M) {
        case 1: {bd_idx = select_bidomain(domains, left, lgrade, current.size()); break;}
        case 2: {bd_idx = select_bidomain(domains, left, V,   current.size()); break;}
        default: break;
        }

        if (bd_idx == -1) return -1;
        bd_ptr = &domains[bd_idx];

        switch (M) {
        case 1: {tmp_idx = selectV_index(left, lgrade, bd_ptr->l, bd_ptr->left_len); break;}
        case 2: {tmp_idx = selectV_index(left, V,     bd_ptr->l, bd_ptr->left_len); break;}
        default: break;
        }
        v = left[bd_ptr->l + tmp_idx];
    }
    /* ---------- 2. 把 v 从域里删掉 ---------- */
    if (vertexs.empty())
        remove_vtx_from_array(left, bd_ptr->l, bd_ptr->left_len, tmp_idx);
    else
        remove_vtx_from_array(left, bd_ptr->l, bd_ptr->left_len,
                              std::find(left.begin() + bd_ptr->l,
                                        left.begin() + bd_ptr->l + bd_ptr->left_len,
                                        v) - (left.begin() + bd_ptr->l));

    /* ---------- 3. 原分支逻辑（完全不变） ---------- */
    vector<int> wselected(g1.n, 0);
    bd_ptr->right_len--;

    int len = bd_ptr->right_len;
    for (int i = len; i >= 0; --i) {
        if (i != 0) {
            num++;
            if (num > Maxnum) { M = M % 2 + 1; num = 0; }
        }
        switch (M) {
        case 1: tmp_idx = selectW_index(right, rgrade, bd_ptr->r+bd_ptr->right_len-len, len + 1, wselected, v); break;
        case 2: tmp_idx = selectW_index(right, Q[v],   bd_ptr->r+bd_ptr->right_len-len, len + 1, wselected, v); break;
        default: break;
        }
        w = right[bd_ptr->r + tmp_idx];
        wselected[w] = 1;
        std::swap(right[bd_ptr->r + tmp_idx], right[bd_ptr->r + bd_ptr->right_len]);
       
        unsigned int cur_len = current.size();
        vector<dal::Bidomain> new_domains;
        switch (M) {
        case 1: new_domains = rewardfeed_RL(domains, bd_idx, current, left, right, lgrade, rgrade, g0, g1, v, w,
                                            g0_matched, g1_matched, arguments.directed || arguments.edge_labelled); break;
        case 2: new_domains = rewardfeed(domains, bd_idx, current, g0_matched, g1_matched, left, right, V, Q[v], g0, g1, v, w,
                                         arguments.directed || arguments.edge_labelled); break;
        default: break;
        }

        dl++;
        int t = solve(g0, g1, V, lgrade, rgrade, Q, incumbent, current, g0_matched, g1_matched, new_domains, vertexs, unvertexs, data_unvertexs, left, right, matching_size_goal,current_best);
        while (current.size() > cur_len) {
            dal::VtxPair pr = current.back();
            current.pop_back();
            g0_matched[pr.v] = 0;
            g1_matched[pr.w] = 0;
        }
        if(t!=-1){
            if(t!=v) return t;
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
    bd_ptr->right_len++;
    if (bd_ptr->left_len == 0) remove_bidomain(domains, bd_idx);
    unvertexs.insert(lower_bound(unvertexs.begin(),unvertexs.end(),v),v);
    int t = solve(g0, g1, V, lgrade, rgrade, Q, incumbent, current, g0_matched, g1_matched, domains, vertexs, unvertexs, data_unvertexs, left, right, matching_size_goal,current_best);
    return t;
}

static void erase_g0_vertices(dalGraph::Graph& g0,
                              const dalGraph::Graph& g1,
                              vector<dal::Bidomain>& domains,
                              vector<int>& left,
                              std::vector<int>& v,
                              const std::set<int>& uv,
                              vector<int>& g0_matched)
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
    vector<dal::Bidomain> new_domains;
    vector<int> new_left;
    new_left.reserve(old_left.size());

    int pos = 0; // write position in new_left
    for (const dal::Bidomain &bd : domains) {
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

        new_domains.emplace_back(pos, bd.r, cnt, bd.right_len, bd.is_adjacent);
        pos += cnt;
    }

    // swap in the rebuilt arrays
    left.swap(new_left);
    domains.swap(new_domains);

     /* 3. 更新Graph */
    dalGraph::Graph new_g0(old_n);
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

    dalGraph::pack_leaves(new_g0);

    // move new graph into g0
    g0 = std::move(new_g0);


}

static vector<dal::VtxPair> mcs(dalGraph::Graph & g0, const dalGraph::Graph & g1, std::vector<int> v, std::set<int> uv, const int id) {
    vector<int> left;  // the buffer of vertex indices for the left partitions
    vector<int> right;  // the buffer of vertex indices for the right partitions

    vector<int> g0_matched(g0.n, 0);
    vector<int> g1_matched(g1.n, 0);
    //*******
    Maxnum=2*std::min(g0.n,g1.n);
    //******
    vector<gtype> lgrade(g0.n,0);
    vector<gtype> V(g0.n, 0);
    vector<gtype> rgrade(g1.n,0);
    vector<vector<gtype>> Q(g0.n, vector<gtype> (g1.n, 0));


    auto domains = vector<dal::Bidomain> {};

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
            if (g0.label[i]==label)
                left.push_back(i);
        for (int i=0; i<g1.n; i++)
            if (g1.label[i]==label)
                right.push_back(i);

        int left_len = left.size() - start_l;
        int right_len = right.size() - start_r;
        domains.push_back({start_l, start_r, left_len, right_len, false});
    }



    erase_g0_vertices(g0, g1, domains, left, v, uv, g0_matched);

    vector<dal::VtxPair> incumbent;
    nodes = 0;
    cutbranches = 0;
    vector<int> uvs(uv.begin(),uv.end());
    sort(uvs.begin(),uvs.end());
    int current_best = 0;
    if (arguments.big_first) {
        for (int k=0; k<g0.n; k++) {
            unsigned int goal = g0.n - k;
            auto left_copy = left;
            auto right_copy = right;
            auto domains_copy = domains;
            vector<dal::VtxPair> current;
            solve(g0, g1, V,lgrade,rgrade, Q, incumbent, current, g0_matched, g1_matched, domains_copy, v, uvs, uvs, left_copy, right_copy, goal,current_best);
            if (incumbent.size() == goal) break;
            if (!arguments.quiet) cout << "Upper bound: " << goal-1 << std::endl;
        }

    } else {
        vector<dal::VtxPair> current;
        solve(g0, g1, V, lgrade,rgrade, Q, incumbent, current, g0_matched, g1_matched, domains, v, uvs, uvs, left, right, 1,current_best);
    }
    //std::cout << "id " << id << "debug dal   " <<  check_sol(g0, g1, incumbent) << " incumbent " << incumbent.size() << " time " << SolveTime << " totaltime " << GetSolveTime() << " " << thread_should_stop.load() << " " << mcs_should_stop.load() << std::endl;
    return incumbent;
}

vector<int> calculate_degrees(const dalGraph::Graph & g) {
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

int sum(const vector<int> & vec) {
    return std::accumulate(std::begin(vec), std::end(vec), 0);
}

void dal::init() {
    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    struct dalGraph::Graph g0 = dalGraph::readGraph(arguments.filename1, format, arguments.directed,
            arguments.edge_labelled, arguments.vertex_labelled);
    struct dalGraph::Graph g1 = dalGraph::readGraph(arguments.filename2, format, arguments.directed,
            arguments.edge_labelled, arguments.vertex_labelled);

    vector<int> g0_deg = calculate_degrees(g0);
    vector<int> g1_deg = calculate_degrees(g1);

    vv0.resize(g0.n);
    std::iota(std::begin(vv0), std::end(vv0), 0);
    bool g1_dense = sum(g1_deg) > g1.n*(g1.n-1);
    std::stable_sort(std::begin(vv0), std::end(vv0), [&](int a, int b) {
        return g1_dense ? (g0_deg[a]<g0_deg[b]) : (g0_deg[a]>g0_deg[b]);
    });

    vv1.resize(g1.n);
    std::iota(std::begin(vv1), std::end(vv1), 0);
    bool g0_dense = sum(g0_deg) > g0.n*(g0.n-1);
    std::stable_sort(std::begin(vv1), std::end(vv1), [&](int a, int b) {//????????????????????????????????????????????????????
        return g0_dense ? (g1_deg[a]<g1_deg[b]) : (g1_deg[a]>g1_deg[b]);
    });

    g0_sorted = dalGraph::induced_subgraph(g0, vv0);
    g1_sorted = dalGraph::induced_subgraph(g1, vv1);


    dalGraph::pack_leaves(g0_sorted);
    dalGraph::pack_leaves(g1_sorted);
}


static void worker_thread(int id) {
    while(1){
        dal::TaskPackage pkg;
        {
            /* 阻塞等任务 */
            std::unique_lock<std::mutex> lk(g_qmtx);
            g_qcv.wait(lk, [] { return !g_taskQ.empty() || thread_should_stop.load();});
            if (thread_should_stop.load() && g_taskQ.empty()) break;
            pkg = std::move(g_taskQ.front());
            g_taskQ.pop();
        }      
        /* ---- 跑 mcs ---- */
        vector<dal::VtxPair> solution = mcs(pkg.g0, pkg.g1, pkg.vertexs, pkg.unvertexSet, id);
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
void dal::run(Worker &worker){
    thread_should_stop = false;
    thread wthr(worker_thread,worker.id_);

    int msg;
    MPI_Irecv(&msg,1,MPI_INT,worker.masterId_,TAG_CMD,MPI_COMM_WORLD,&worker.msgHandlers_.recvRequest);
    MPI_Status status;
    while(1){
        int flag = 0;
        MPI_Test(&worker.msgHandlers_.recvRequest, &flag, &status);
        if (flag) {
            int current_msg = msg;
            if (current_msg == msg::SendInformation) {
                mcs_should_stop = false;
                nodes = 0;
                cutbranches = 0;

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