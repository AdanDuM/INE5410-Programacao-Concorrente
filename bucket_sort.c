#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void print_vector(unsigned *vector, unsigned size);
void generate_vector(unsigned *vector, unsigned size);
void create_buckets(struct Bucket buckets[], unsigned size, unsigned nbuckets);
void distribute(unsigned *vector, unsigned size);
void quicksort();

struct Bucket {
	unsigned begin_range;
	unsigned end_range;
};

int main(int argc, char **argv) {
	unsigned tam_vet;
	unsigned nbuckets;
	unsigned nprocs;
	unsigned imprimir;

	unsigned *vector;

	unsigned size, rank;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (argc < 4) {
		printf("Numero de entradas invalida.\n");
		exit(1);
	}

	tam_vet = atoi(argv[1]);
	nbuckets = atoi(argv[2]);
	nprocs = atoi(argv[3]);

	if (nbuckets > tam_vet) {
		printf("Erro: Numero de buckets maior que tamanho do vetor.\n");
		exit(1);
	}

	if (argc == 5 && argv[4] == 1) {
		imprimir = 1;
	}

	if (rank == 0) {
		unsigned i;

		struct Bucket buckets[nbuckets];

		//buckets = malloc(sizeof(unsigned) * nbuckets);

		vector = (unsigned*)malloc(sizeof(unsigned) * tam_vet);

		generate_vector(vector, tam_vet);

		if (imprimir) {
			print_vector(vector, tam_vet);
		}

		create_buckets(buckets, tam_vet, nbuckets);
		distribute(vector, tam_vet);

		for (i = 0; i < buckets; ++i) {
			MPI_Send();
		}

		for (i = 0; i < buckets; ++i) {
			MPI_Recv();
		}

		if (imprimir) {
			print_vector(vector, tam_vet);
		}
	} else {
		MPI_Recv();

		vector = malloc(sizeof(unsigned) * tam_vet);

		quicksort();

		MPI_Send();
	}

	free(vector);

	MPI_Finalize();

	return 0;
}

void print_vector(unsigned *vector, unsigned size) {
	unsigned i;
	for (i = 0; i < size; ++i) {
		printf("%10d", vector[i]);
	}
}

void generate_vector(unsigned *vector, unsigned size) {
	srand(time(NULL));

	unsigned i;
	for (i = 0; i < size; ++i) {
		vector[i] = rand() % size;
	}
}

void create_buckets(struct Bucket *buckets, unsigned size, unsigned nbuckets) {
	unsigned i;
	unsigned range = size / nbuckets;
	unsigned rest = size % nbuckets;

	buckets = (struct Bucket*)malloc(sizeof(struct Bucket) * nbuckets);

	if (rest == 0) {
		buckets[0].begin_range = 0;
		buckets[0].end_range = buckets[0].begin_range + range - 1;
		++rest;
	} else {
		buckets[0].begin_range = 0;
		buckets[0].end_range = buckets[0].begin_range + range;
	}

	for (i = 1; i < rest; ++i) {
		buckets[i].begin_range = buckets[i - 1].end_range + 1;
		buckets[i].end_range = buckets[i].begin_range + range + 1;
	}
	for (i = rest; i < size; ++i) {
		buckets[i].begin_range = buckets[i - 1].end_range + 1;
		buckets[i].end_range = buckets[i].begin_range + range;
	}
}

void distribute(unsigned *vector, unsigned size) {
	unsigned i;
	for (i = 0; i < size; ++i) {
		if () {}
	}
}

void quicksort() {}
