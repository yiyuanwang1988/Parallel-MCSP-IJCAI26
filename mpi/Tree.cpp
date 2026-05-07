#include"Tree.hpp"


Tree::Tree(){}
Tree::~Tree(){}

bool Cmp::operator()(const Node* a, const Node* b) const {
    if (a->time == b->time) return a->nodeId < b->nodeId;
    return a->time < b->time;
}

void Tree::BuildInitNodes(){
    rootNode_ = new Node(nullptr,1);
    BranchNode(rootNode_);
}

void Tree::SetNodeStatus(Node* node,NodeStatus status){
    if(status == NodeStatus::Waiting) {
        waitingNodes.insert(node);
    }
    else if(status == NodeStatus::Running){
        waitingNodes.erase(node);
        runningNodes.insert(node);
        if(node->leftNode == nullptr && node->rightNode == nullptr)
        BranchNode(node);
    }
    else if(status == NodeStatus::End){
        waitingNodes.erase(node);
        runningNodes.erase(node);
        
    }
    node->nodeStatus = status;
}

void Tree::BranchNode(Node* node) {
    Node* leftnode = new Node(node,node->nodeId*2);
    leftnode->Vertexs = node->Vertexs;
    leftnode->Unvertexs = node->Unvertexs;

    Node* rightnode = new Node(node,node->nodeId*2+1);
    rightnode->Vertexs = node->Vertexs;
    rightnode->Unvertexs = node->Unvertexs;

    int num = -1, v = -1;
    for (int i = 0; i < g0n; i++) {
        int flag = 0;
        for (auto vv : node->Vertexs) {
            if (i == vv) {
                flag = 1;
                break;
            }
        }
        for (auto vv : node->Unvertexs) {
            if (i == vv) {
                flag = 1;
                break;
            }
        }
        if (flag) continue;
        if (score[i] > num) {
            num = score[i];
            v = i;
        }
    }
    if(v==-1) return;
    leftnode->Vertexs.push_back(v);
    rightnode->Unvertexs.push_back(v);

    node->leftNode = leftnode;
    node->rightNode = rightnode;

    //cout << "BranchNode Waiting " << leftnode->nodeId  << " " << rightnode->nodeId  << endl;
                        
    SetNodeStatus(leftnode, NodeStatus::Waiting);
    SetNodeStatus(rightnode, NodeStatus::Waiting);

}

Node *Tree::GetNodeToRun(){
  Node *resNode = nullptr;
  if (!waitingNodes.empty())
  {
    resNode = *(waitingNodes.begin());
    //cout << "GetNodeToRun Running " << resNode->nodeId  << endl;
      
    SetNodeStatus(resNode, NodeStatus::Running);
  }
  return resNode;
}