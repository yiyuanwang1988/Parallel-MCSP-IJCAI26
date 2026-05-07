#pragma once
#include <mpi.h>
#include<bits/stdc++.h>
using namespace std;
const int N = 1e5+7;

const int TAG_CMD = 0;
const int TAG_RESULT = 1;
const int TAG_DATA = 2; // 用于发送节点信息
const int TAG_NEW_ANS = 3;
const int TAG_HEADER = 4;
const int TAG_UPDATE_ANS = 5;
const int TAG_RECV_DATA1 = 6;
const int TAG_RECV_DATA2 = 7;
const int TAG_RECV_DATA3 = 8;
inline int g0n;
inline int cutoff_time = 0;
inline std::vector<int> vv0;
inline std::vector<int> vv1;
extern double start_time;
extern std::atomic<int> ans;
inline int getAns() {
    return ans.load(std::memory_order_relaxed);
}
extern int FreeWorkerNum;
extern int score[N];
inline std::atomic<bool> mcs_should_stop{ false };
inline std::atomic<bool> thread_should_stop{false};

inline double GetTime() {
    using namespace std::chrono;
    return duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count() / 1000.0;
}
inline double GetSolveTime(){
    return GetTime() - start_time;
}
inline int* change(std::vector<int> v){
    int *a = new int[v.size()];
    for(int i=0;i<v.size();i++) a[i]=v[i];
    return a;
}
enum WorkerType {
    DAL = 1,
    DSB,
    RRS,
};

enum msg {
    SendInformation = 1,
    UpdateAns,
    UpdateStatue,
    UpdateTime,
    Endall,
    RecvData,
    SendData,
};

enum MessageType {
	SOLVE_TASK = 1,
	SHARE_INFO,
	TERMINATE,
};
struct Task {
    std::vector<int> vertexs;
    std::set<int> unvertexSet;
};
struct Data{
    int v,w,ans,num;
    vector<int> vs;
};
extern vector<Data> datas;
extern vector<Data> new_datas;
inline bool cmpdata(const Data &x, const Data &y){
    if (x.v != y.v) return x.v < y.v;
    return x.w < y.w;
}
inline int checkdata(vector<int>a, vector<int>b){ // a带复用信息 b带data
    set<int>s;
    set_intersection(a.begin(),a.end(),b.begin(),b.end(),inserter(s,begin(s)));
    if(s.size() == a.size()) return 0;
    if(s.size() == b.size()) return a.size()-b.size();
    return -1;
}

// struct Gadjmat{
//     vector<vector<int>> adjmat;
//     vector<int> scores;
// };
// extern Gadjmat adjmatg0,adjmatg1;
// inline bool check_v_w(int v,int w){
//     if(adjmatg0.scores[v] > adjmatg1.scores[w]) return false;
//     int i = 0;
//     for(int u:adjmatg0.adjmat[v]){
//         if(adjmatg0.scores[u]>adjmatg1.scores[i]) return false;
//         i++;
//     }
//     return true;
// }



struct AnsAndTime {
    int ans;
    double tt;
};
inline MPI_Datatype MPI_ANS_AND_TIME;
inline MPI_Datatype MPI_DATA;

struct MessageHandler {
	MPI_Request sendRequest;
	MPI_Request recvRequest;
	int idWorker;
};

enum ProblemStatus
{
  Optimal,   // Optimal
  Cut,  // Cut
  Solve,  // Solve
  UnSolve // UnSolve
};

enum NodeStatus
{
  New,
  Waiting,
  Running,
  BranchedWaiting,
  BranchedRunning,
  End
};

enum Heuristic { min_max, min_product };
struct Arguments{
    bool quiet;
    bool verbose;
    bool dimacs;
    bool lad;
    bool connected;
    bool directed;
    bool edge_labelled;
    bool vertex_labelled;
    bool big_first;
    Heuristic heuristic;
    char *filename1;
    char *filename2;
    int timeout;
    int arg_num;
    int solver;
};
extern Arguments arguments;

inline void set_default_arguments() {
    arguments.quiet = false;
    arguments.verbose = false;
    arguments.dimacs = false;
    arguments.lad = false;
    arguments.connected = false;
    arguments.directed = false;
    arguments.edge_labelled = false;
    arguments.vertex_labelled = false;
    arguments.big_first = false;
    arguments.filename1 = NULL;
    arguments.filename2 = NULL;
    arguments.timeout = 1800;
    arguments.arg_num = 0;
}
static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

inline void read_params(int argc, char **argv){
	for (int i = 1; i < argc; i++){
      
        if (std::string(argv[i]) == "min_max") {
            arguments.heuristic = min_max;
            continue;
        }  
        else if (std::string(argv[i]) == "min_product") {
            arguments.heuristic = min_product;
            continue;
        }
        else if (std::string(argv[i]) == "-P") {
            i++;
            arguments.filename1 = argv[i];
            continue;
        }
        else if (std::string(argv[i]) == "-T") {
            i++;
            arguments.filename2 = argv[i];
            continue;
        }
        else if (std::string(argv[i]) == "-d") {
            if (arguments.lad)
                fail("The -d and -l options cannot be used together.\n");
            arguments.dimacs = true;
            continue;
        }
        else if (std::string(argv[i]) == "-l") {
            if (arguments.dimacs)
                fail("The -d and -l options cannot be used together.\n");
            arguments.lad = true;
            continue;
        }
        else if (std::string(argv[i]) == "-q") {
            arguments.quiet = true;
            continue;
        }
        else if (std::string(argv[i]) == "-v") {
            arguments.verbose = true;
            continue;
        }
        else if (std::string(argv[i]) == "-c") {
            if (arguments.directed)
                fail("The connected and directed options can't be used together.");
            arguments.connected = true;
            continue;
        }
        else if (std::string(argv[i]) == "-i") {
            if (arguments.connected)
                fail("The connected and directed options can't be used together.");
            arguments.directed = true;
            continue;
        }
        else if (std::string(argv[i]) == "-a") {
            if (arguments.vertex_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.edge_labelled = true;
            arguments.vertex_labelled = true;
            continue;
        }
        else if (std::string(argv[i]) == "-x") {
            if (arguments.edge_labelled)
                fail("The -a and -x options can't be used together.");
            arguments.vertex_labelled = true;
            continue;
        }
        else if (std::string(argv[i]) == "-b") {
            arguments.big_first = true;
            continue;
        }
        else if (std::string(argv[i]) == "-t") {
            i++;
            arguments.timeout = std::stoi(argv[i]);
            continue;
        }
        else if (std::string(argv[i]) == "-S") {
            i++;
            arguments.solver = std::stoi(argv[i]);
            // if (std::stoi(argv[i]) == 1) arguments.solver = "Mcsplit-DSB";
            // else if (std::stoi(argv[i]) == 2) arguments.solver = "RRSplit";
            // else if (std::stoi(argv[i]) == 3) arguments.solver = "Mcsplit-DAL";
            continue;
        }
	}  
}


inline void read(int argc,char** argv){
        set_default_arguments();
		read_params(argc,argv);
}