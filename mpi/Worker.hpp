#pragma once

#include"Message.hpp"
#include"../McSplitDSB/mcsp.hpp"
#include"../RRSplit/mcsp.hpp"
#include"../McSplitDAL/mcsp.hpp"
extern int cutoff_time;
extern double start_time;
extern double GetTime();

inline int CurrentSolveTime = 0x3f3f3f3f;

class Worker{
    private:

    public:
        unsigned id_;
        unsigned masterId_;
        unsigned type_;
        MessageHandler msgHandlers_;
        Worker(unsigned id_, unsigned type_);
        Worker() = default;
        void run();
};
