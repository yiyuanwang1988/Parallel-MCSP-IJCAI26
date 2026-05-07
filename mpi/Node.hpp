#pragma once

extern double GetTime();
#include"Message.hpp"
class Node{
    public:
        Node(Node *parentNode_,int id);
        ~Node();
        long long nodeId;
        std::vector<int> Vertexs;
        std::vector<int> Unvertexs;
        int worker;
        int time;
        int nodeStatus;
        Node* parentNode;
        Node* leftNode, *rightNode;
    private:
};