#include"Worker.hpp"

Worker::Worker(unsigned id,unsigned type){
    id_ = id;
    masterId_ = 0;
    type_ = type;
}


void Worker::run(){
    if(type_ == WorkerType::DSB){
        dsb* solver = new dsb();
        solver->init();
        solver->run(*this);
    }
    else if(type_ == WorkerType::RRS){
        rrs* solver = new rrs();
        solver->init();
        solver->run(*this);
    }
    else if(type_ == WorkerType::DAL){
        dal* solver = new dal();
        solver->init();
        solver->run(*this);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}