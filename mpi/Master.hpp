#pragma once

#include"Message.hpp"
#include"../mpi/Tree.hpp"
#include"../mpi/Node.hpp"
#include"Worker.hpp"
extern int cutoff_time;
extern double start_time;
extern double GetTime();

class Tree;

class Master {
	private :
		unsigned id_;
		unsigned worldSize_;
		std::vector<MessageHandler> msgHandlers_;
		std::vector<MPI_Request> SendInformationHead;
		std::vector<MPI_Request> SendInformationBody;
		std::vector<MPI_Request> SendAns;
		std::vector<MPI_Request> SendStatus;
	public :
		Master(unsigned id, unsigned worldSize);
		void run();
		void Mpi_Msg(int order, int recv);
		void Mpi_SendInformation(Node* node,int recv);
		void Mpi_UpdateAns();
		void Mpi_UpdateStatus(Node* node);
		void Mpi_UpdateTime();
		void Mpi_visited(Node* parentnode,Node* node);
};

