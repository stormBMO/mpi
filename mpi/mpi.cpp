#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include <iostream>

#define N 1000

using namespace std;

void multiplyMatrix(double* a, double* b, double* c, int rank, int size) {
	for (int i = 0; i < N; i++) {
		for (int j = 0 + rank; j < N; j += size) {
			for (int k = 0; k < N; k++) {
				c[i * N + j] += a[i * N + k] * b[k * N + j];
			}
		}
	}
}


int main(int argc, char* argv[]) {
	double* a, * b, * c;
	double* result = 0;
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	a = b = c = new double[N * N];

	if (rank == 0) {
		result = new double[N * N];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				a[i * N + j] = rand() / 10;
				b[i * N + j] = rand() / 10;
			}
		} 
		for (int i = 0; i < N * N; i++) {
			c[i] = 0;
			result[i] = 0;
		}
		double startTime = MPI_Wtime();
		for (int count = 1; count < size; count++) {
			MPI_Send(a, N * N, MPI_DOUBLE, count, 1, MPI_COMM_WORLD);
			MPI_Send(b, N * N, MPI_DOUBLE, count, 2, MPI_COMM_WORLD);
			MPI_Send(c, N * N, MPI_DOUBLE, count, 3, MPI_COMM_WORLD);
		}
		multiplyMatrix(a, b, c, rank, size);
		MPI_Reduce(c, result, N * N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		double endTime = MPI_Wtime();
		cout << "Time result: " << endTime - startTime << endl;
	}
	else {
		MPI_Recv(a, N * N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(b, N * N, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(c, N * N, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		multiplyMatrix(a, b, c, rank, size);
		MPI_Reduce(c, result, N * N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize();
}