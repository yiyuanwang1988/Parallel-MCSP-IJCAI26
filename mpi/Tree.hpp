#pragma once
#include"Node.hpp"
#include"../mpi/Message.hpp"

class Node;
class Master;

struct Cmp {
    bool operator()(const Node* a, const Node* b) const; // 只声明
};

class Tree
  {
  private:
      
  public:
      std::set<Node *, Cmp> waitingNodes;   // 等待队列：存储待处理的节点
      std::set<Node *, Cmp> runningNodes;  // 运行队列：存储正在处理的节点
    Tree();
    ~Tree();
    Node *rootNode_;

    void BuildInitNodes();
    Node* GetNodeToRun();
    void SetNodeStatus(Node* node,NodeStatus status);
    void BranchNode(Node *node);
};