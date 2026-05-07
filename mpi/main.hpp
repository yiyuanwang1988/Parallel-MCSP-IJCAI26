#pragma once
#include"Master.hpp"
#include"Worker.hpp"
#include"Message.hpp"
Arguments arguments;
std::atomic<int> ans{0};
double start_time = GetTime();
int FreeWorkerNum;
int score[N];
extern void read(int argc,char** argv);
vector<int> visited[6000];
vector<int> need_visited[6000];
vector<Data> datas,new_datas;
