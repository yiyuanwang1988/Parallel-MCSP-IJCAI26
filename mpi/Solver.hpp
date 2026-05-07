#pragma once
class Worker;          // 前向声明
class Solver {
public:
    Solver() = default; 
    ~Solver() = default;
    virtual void init() = 0;
    virtual void run(Worker& worker) = 0;   // 引用传参
};