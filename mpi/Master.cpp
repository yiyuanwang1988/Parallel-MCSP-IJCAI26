#include"Master.hpp"
Tree tree;
static double SolveTime = 0;
static std::vector<int> WorkerStatus;
static std::vector<Node*> Worker;
Master::Master(unsigned id, unsigned worldSize){
    id_ = id;
    worldSize_ = worldSize;
    WorkerStatus.assign(worldSize_, 0);
    Worker.assign(worldSize_, nullptr);
    msgHandlers_.resize(worldSize_);
    SendInformationHead.resize(worldSize);
    SendInformationBody.resize(worldSize);
    SendAns.resize(worldSize);
    SendStatus.resize(worldSize);
}

static std::vector<std::array<int, 2>> head;
static std::vector<std::vector<int>> body;
void Master::Mpi_SendInformation(Node* node, int recv) {
    if (head.size() != worldSize_) {
        head.resize(worldSize_);
    }
    if (body.size() != worldSize_) {
        body.resize(worldSize_);
    }
    int n = node->Vertexs.size();
    int m = node->Unvertexs.size();
    head[recv][0] = n;
    head[recv][1] = m;
    MPI_Send(head[recv].data(), 2, MPI_INT, recv, TAG_HEADER, MPI_COMM_WORLD);

    body[recv].clear();
    body[recv].reserve(n + m);
    if(n) body[recv].insert(body[recv].end(), node->Vertexs.begin(), node->Vertexs.end());
    if(m) body[recv].insert(body[recv].end(), node->Unvertexs.begin(), node->Unvertexs.end());

    MPI_Send(body[recv].data(), n + m, MPI_INT, recv, TAG_DATA, MPI_COMM_WORLD);
}

void Master::Mpi_UpdateAns() {
    for (int i = 1; i < worldSize_; i++) {
        MPI_Isend(&ans, 1, MPI_INT, i, TAG_NEW_ANS, MPI_COMM_WORLD, &SendAns[i]);
    }
}

void Master::Mpi_UpdateStatus(Node *node){
    Node *leftnode = node->leftNode;
    Node *rightnode = node->rightNode;
    if (leftnode != nullptr && rightnode != nullptr) {
        if (leftnode->nodeStatus != NodeStatus::End) {
            Mpi_Msg(msg::UpdateStatue, leftnode->worker);
            tree.SetNodeStatus(leftnode, NodeStatus::End);
            Mpi_UpdateStatus(leftnode);
        }
        if (rightnode->nodeStatus != NodeStatus::End){
            Mpi_Msg(msg::UpdateStatue, rightnode->worker); 
            tree.SetNodeStatus(rightnode, NodeStatus::End);
            Mpi_UpdateStatus(rightnode);
        } 
    }
}



void Master::Mpi_Msg(int order, int recv) {
    if (recv) {
        if(recv == -1) return;
        //MPI_Isend(&order, 1, MPI_INT, recv, 0, MPI_COMM_WORLD,&msgHandlers_[recv].sendRequest);
        MPI_Send(&order, 1, MPI_INT, recv, TAG_CMD, MPI_COMM_WORLD);
    }
    else 
        for(int i = 1; i < worldSize_; i ++) {
            //MPI_Isend(&order, 1, MPI_INT, i, 0, MPI_COMM_WORLD,&msgHandlers_[i].sendRequest);
            MPI_Send(&order, 1, MPI_INT, i, TAG_CMD, MPI_COMM_WORLD);
        }
} 

void Master::Mpi_UpdateTime() {
    new_datas.clear();
    vector<Data> tmp_datas;
    std::vector<AnsAndTime> res(worldSize_);
    for(int i = 1; i < worldSize_; i ++) {
        tmp_datas.clear();
        //获取信息并处理
        int num = 0;
        MPI_Recv(&num,1,MPI_INT,i,TAG_RECV_DATA1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        if(num){
            std::vector<int> a(num*4);
            MPI_Recv(a.data(), 4 * num, MPI_INT, i, TAG_RECV_DATA2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int lenb = 0;
            for(int j = 0; j < 4*num; j += 4){
                Data data = {a[j+0],a[j+1],a[j+2],a[j+3]};
                data.vs.clear();
                lenb += a[j+3];
                tmp_datas.push_back(data);
            }
            std::vector<int> b(lenb);
            MPI_Recv(b.data(), lenb, MPI_INT, i, TAG_RECV_DATA3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            lenb = 0;
            for(auto &data:tmp_datas){
                int len = data.num;
                while(len--) data.vs.push_back(b[lenb++]);
                auto it = lower_bound(datas.begin(),datas.end(),data,cmpdata);
                if(it == datas.end()) {
                    datas.insert(it,data);
                    new_datas.push_back(data);
                }
                else {
                    if(it->v == data.v && it->w == data.w && it->vs == data.vs){
                        if(it->ans < data.ans){
                            it->ans = data.ans;
                            new_datas.push_back(data);                            
                        }                           
                    }
                    else {
                        datas.insert(it,data);
                        new_datas.push_back(data);
                    }
                }
            }                  
        } 
        MPI_Recv(&res[i], 1, MPI_ANS_AND_TIME, i, TAG_UPDATE_ANS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(res[i].ans > ans){
			ans = res[i].ans;
			SolveTime = res[i].tt;
		}

	}
    if (GetSolveTime() < cutoff_time) {
        Mpi_Msg(msg::UpdateAns, 0);
        MPI_Bcast(&ans, 1, MPI_INT, id_, MPI_COMM_WORLD);

        Mpi_Msg(msg::SendData,0);
        int res =new_datas.size();
        MPI_Bcast(&res,1,MPI_INT,id_,MPI_COMM_WORLD);
        if(res){
            std::vector<int> a(res * 4);
            std::vector<int> b;
            int lena=0,lenb=0;
            for(auto &u : new_datas){
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
                for(auto &k:u.vs){
                    b[lenb++]=k;
                }
            }
            MPI_Bcast(a.data(), res * 4, MPI_INT, id_, MPI_COMM_WORLD);
            MPI_Bcast(b.data(), lenb, MPI_INT, id_, MPI_COMM_WORLD);
        }

    }
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
void readnum(){
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
    g0n = g0.n;
}

void Master::run(){
    readnum();
    tree.BuildInitNodes();
    std::vector<AnsAndTime> res(worldSize_);
    MPI_Status status;
    int should_stop = 0;
    std::vector<std::pair<int,int>> ansVertexs;
    while(GetSolveTime() < cutoff_time || FreeWorkerNum+1 < worldSize_){
        if(!should_stop && GetSolveTime() < cutoff_time){
            for(int i=1;i<worldSize_;i++){
                if(WorkerStatus[i] == 0){
                    Node* node = tree.GetNodeToRun();
                    while(node!=nullptr&&g0n - node->Unvertexs.size() <= ans){
                        //cout << "run not node End " << node->nodeId << endl;
                        Mpi_UpdateStatus(node);
                        tree.SetNodeStatus(node,NodeStatus::End);
                        node = tree.GetNodeToRun();                      
                    }
                    if(node == nullptr) {
                        break;
                    }
                    //node->time = GetSolveTime();
                    node->worker = i;
                    Worker[i] = node;
                    WorkerStatus[i] = 1;
                    FreeWorkerNum--;
                    Mpi_Msg(msg::SendInformation,i);
                    Mpi_SendInformation(node,i);
                    MPI_Irecv(&res[i], 1, MPI_ANS_AND_TIME, i, TAG_NEW_ANS, MPI_COMM_WORLD, &msgHandlers_[i].recvRequest);
                }
            }
        }

        int nowtime = int(GetSolveTime());
        if (!should_stop && nowtime % 10 == 0 && nowtime) {
            Mpi_Msg(msg::UpdateTime,0);
            Mpi_UpdateTime();
            sleep(1);
        }
        for(int i=1;i<worldSize_;i++){
            int flag = 0;
            MPI_Test(&msgHandlers_[i].recvRequest,&flag,&status);
            // MPI_Wait(&msgHandlers_[i].recvRequest,&status);
            if(flag && WorkerStatus[i]){
                FreeWorkerNum++;
                WorkerStatus[i] = 0; 
                vector<int> vs(res[i].ans*2);
                MPI_Recv(vs.data(), res[i].ans *2, MPI_INT, i, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
           
                if(res[i].ans >= ans){
                    ans.store(res[i].ans, std::memory_order_relaxed);
                    for (int j = 0; j < res[i].ans * 2; j += 2) score[vs[j]]++;
                    // if (ans == g0_sorted.n) break;
                    ansVertexs.clear();
                    for (int j = 0; j < res[i].ans * 2; j += 2) {
                        ansVertexs.push_back({ vs[j],vs[j + 1] });
                    }  

                    SolveTime = min(SolveTime,res[i].tt);
                    if(!should_stop && GetSolveTime() < cutoff_time){
                        Mpi_Msg(msg::UpdateAns,0);
                        MPI_Bcast(&ans, 1, MPI_INT, id_, MPI_COMM_WORLD);
                    }
                }
                if(!should_stop && GetSolveTime() < cutoff_time){
                    
                    Mpi_UpdateStatus(Worker[i]);
                    tree.SetNodeStatus(Worker[i],NodeStatus::End);

                    Node* parentnode = Worker[i]->parentNode;
                    if (parentnode == nullptr) {
                        continue;
                    }
                    if (parentnode->leftNode->nodeStatus == NodeStatus::End && parentnode->rightNode->nodeStatus == NodeStatus::End) {
                        if (parentnode->nodeStatus != NodeStatus::End) {
                            Mpi_Msg(msg::UpdateStatue, parentnode->worker);
                            tree.SetNodeStatus(parentnode, NodeStatus::End); 
                        }
                    }                        
                }
            }
        }
        if (!should_stop && GetSolveTime() > cutoff_time) {
            Mpi_Msg(msg::Endall,0);
            should_stop = -1;
           
        }
        else if (!should_stop &&  (getAns() >= g0n||(tree.waitingNodes.empty()&&tree.runningNodes.empty()))) {
            Mpi_Msg(msg::Endall,0);
            should_stop = 1;
        }
        if(should_stop && FreeWorkerNum+1 == worldSize_) break;
        
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << arguments.filename1 <<" "<< arguments.filename2;
    cout<<" ans " << getAns() << " time " << SolveTime << " totaltime " << GetSolveTime() << " solved " << should_stop << std::endl;
    //for (auto& [v, w] : ansVertexs) {
    //    std::cout<< vv0[v] << " " << vv1[w] << std::endl;
        //cout<<v<<" "<<w<<endl;
    //}
    return;
}
