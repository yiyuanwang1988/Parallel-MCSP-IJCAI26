#include"Node.hpp"
int NODEID = 1;

Node::Node(Node* parentNode_,int id){
    nodeId = id;
	time = GetTime();
    nodeStatus = NodeStatus::Waiting;
    parentNode = parentNode_;
    worker = -1;
    leftNode = nullptr;
    rightNode = nullptr;
}
Node::~Node(){}

