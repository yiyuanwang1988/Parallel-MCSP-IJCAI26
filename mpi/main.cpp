#include"main.hpp"

int main(int argc, char** argv) {
	int provided = MPI_THREAD_SINGLE;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	if (provided < MPI_THREAD_MULTIPLE) {
		std::cerr << "MPI implementation does not support MPI_THREAD_MULTIPLE, aborting." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		return 1;
	}

	// MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	read(argc,argv);
	cutoff_time = arguments.timeout;
	start_time = GetTime();

	int          blocklen[] = { 1, 1 };
	MPI_Datatype types[] = { MPI_INT, MPI_DOUBLE };
	MPI_Aint     disp[] = { offsetof(AnsAndTime, ans),offsetof(AnsAndTime, tt) };
	MPI_Type_create_struct(2, blocklen, disp, types, &MPI_ANS_AND_TIME);
	MPI_Type_commit(&MPI_ANS_AND_TIME);

	if (rank == 0) {
		FreeWorkerNum = size-1;
		Master master(rank, size);
		master.run();
	}
	else {
		// Worker worker(rank,(rank+2)%3+1);
		Worker worker(rank,arguments.solver);
		worker.run();
	}


	MPI_Type_free(&MPI_ANS_AND_TIME);
	MPI_Finalize();
}
