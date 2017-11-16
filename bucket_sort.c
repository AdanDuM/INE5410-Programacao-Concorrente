
//  !	Bucket Sort MPI
/*	!
 *  \Copyright (C) 2017  Adan Pereira Gomes
 *  \Released under the GNU General Public License 2.0
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define DEBUG 0	// Mude para 1 caso queira ver os numeros impressos no terminal

typedef struct {
	unsigned *vector;
	unsigned size;
	unsigned begin_range;
	unsigned end_range;
} Bucket;

void print_vector(const unsigned vector[], unsigned size);
void populate_vector(unsigned vector[], unsigned size);
void create_buckets(Bucket buckets[], unsigned size, unsigned nbuckets);
void destroy_buckets(Bucket buckets[], unsigned nbuckets);
void insert(Bucket *bucket, unsigned number);
void distribute_numbers(Bucket buckets[], unsigned nbuckets, unsigned vector[], unsigned size);
void distribute_buckets(Bucket buckets[], unsigned nbuckets);
void receive_sort_send();
void join_buckets(Bucket buckets[], unsigned nbuckets, unsigned vector[], unsigned size);
void quicksort(unsigned vector[], int low, int high);
int partition(unsigned vector[], int low, int high);
void swap(unsigned *a, unsigned *b);

unsigned nprocs, rank;

int main(int argc, char **argv) {
	unsigned tam_vet;
	unsigned nbuckets;
	unsigned imprimir;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, (int*)&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, (int*)&rank);

	if (argc < 3) {
		printf("Numero de entradas invalida.\n");
		exit(1);
	}

	tam_vet = atoi(argv[1]);
	nbuckets = atoi(argv[2]);

	if (nbuckets > tam_vet) {
		printf("Erro: Numero de buckets maior que tamanho do vetor.\n");
		exit(1);
	}

	if (argc == 4 && atoi(argv[3]) == 1) {
		imprimir = 1;
	}

	if (rank == 0) {
		unsigned vector[tam_vet];
		Bucket buckets[nbuckets];

		populate_vector(vector, tam_vet);
		create_buckets(buckets, tam_vet, nbuckets);

		if (imprimir) {
			print_vector(vector, tam_vet);
		}

		distribute_numbers(buckets, nbuckets, vector, tam_vet);

		MPI_Barrier(MPI_COMM_WORLD);
		distribute_buckets(buckets, nbuckets);
		join_buckets(buckets, nbuckets, vector, tam_vet);
		destroy_buckets(buckets, nbuckets);

		if (imprimir) {
			print_vector(vector, tam_vet);
		}
	} else {
		MPI_Barrier(MPI_COMM_WORLD);
		receive_sort_send();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}

void print_vector(const unsigned vector[], unsigned size) {
	unsigned i;
	printf("\nVector = { ");
	for (i = 0; i < size; ++i) {
		printf("%d,", vector[i]);
	}
	printf("\b }");
}

void populate_vector(unsigned vector[], unsigned size) {
	srand(time(NULL));

	unsigned i;
	for (i = 0; i < size; ++i) {
		vector[i] = rand() % size;
	}
}

void create_buckets(Bucket buckets[], unsigned size, unsigned nbuckets) {
	unsigned i = 0;
	unsigned range = size / nbuckets;
	unsigned rest = size % nbuckets;

	if (rest == 0) {
		buckets[0].vector = (unsigned*)malloc(sizeof(unsigned) * 1);
		buckets[0].size	= 0;
		buckets[0].begin_range	= 0;
		buckets[0].end_range	= buckets[0].begin_range + range - 1;
		++rest;
	} else {
		buckets[0].vector = (unsigned*)malloc(sizeof(unsigned) * 1);
		buckets[0].size	= 0;
		buckets[0].begin_range	= 0;
		buckets[0].end_range	= buckets[0].begin_range + range;
	}

	for (i = 1; i < rest; ++i) {
		buckets[i].vector = (unsigned*)malloc(sizeof(unsigned) * 1);
		buckets[i].size	= 0;
		buckets[i].begin_range	= buckets[i - 1].end_range + 1;
		buckets[i].end_range	= buckets[i].begin_range + range;
	}
	for (i = rest; i < nbuckets; ++i) {
		buckets[i].vector = (unsigned*)malloc(sizeof(unsigned) * 1);
		buckets[i].size	= 0;
		buckets[i].begin_range	= buckets[i - 1].end_range + 1;
		buckets[i].end_range	= buckets[i].begin_range + range - 1;
	}
#if DEBUG
	for (i = 0; i < nbuckets; ++i) {
		printf("\nBucket %d = [%d, %d]", i, buckets[i].begin_range, buckets[i].end_range);
	}
#endif
}

void destroy_buckets(Bucket buckets[], unsigned nbuckets) {
	unsigned i;
	for (i = 0; i < nbuckets; ++i) {
		free(buckets[i].vector);
	}
}

void insert(Bucket *bucket, unsigned number) {
	unsigned n = (*bucket).size++;
	(*bucket).vector[n++] = number;
	(*bucket).vector = (unsigned *)realloc((*bucket).vector, sizeof(unsigned) * (n + 1));
}

void distribute_numbers(Bucket buckets[], unsigned nbuckets, unsigned vector[], unsigned size) {
	unsigned i, j;
	for (i = 0; i < size; ++i) {
		for (j = 0; j < nbuckets; ++j) {
			if (vector[i] <= buckets[j].end_range && vector[i] >= buckets[j].begin_range) {
				insert(&buckets[j], vector[i]);
				break;
			}
		}
	}
#if DEBUG
	for (i = 0; i < nbuckets; ++i) {
		printf("\nBucket vector %d (%d) = { ", i, buckets[i].size);
		for (j = 0; j < buckets[i].size; ++j) {
			printf("%d,", buckets[i].vector[j]);
		}
		printf("\b }");
	}
	printf("\n");
#endif
}

void distribute_buckets(Bucket buckets[], unsigned nbuckets) {
	unsigned buckets_received = 0, buckets_sended = 0;
	unsigned process_rank;
	unsigned terminate = 0;

	while (1) {
		if (buckets_sended < nbuckets && buckets[buckets_sended].size < 2) {
			++buckets_sended;
			++buckets_received;
			continue;
		}

		if (buckets_received == nbuckets && buckets_sended == nbuckets) {
			unsigned remaining_process = nprocs;
			terminate = 1;

			while (remaining_process-- > 1) {
				MPI_Recv(&process_rank, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Send(&terminate, 1, MPI_UNSIGNED, process_rank, 666, MPI_COMM_WORLD);
			}

			break;
		} else {
			unsigned bucket_buffer[2];	// ID, Size
			unsigned option;

			MPI_Recv(&process_rank, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&terminate, 1, MPI_UNSIGNED, process_rank, 666, MPI_COMM_WORLD);
			MPI_Recv(&option, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			if (option == 1 && buckets_sended < nbuckets) {	// Send Buckets
#if DEBUG
				printf("\nSending Bucket vector %d (%d)", buckets_sended, buckets[buckets_sended].size);
#endif
				bucket_buffer[0] = buckets_sended;
				bucket_buffer[1] = buckets[buckets_sended].size;

				MPI_Send(&bucket_buffer[0], 1, MPI_UNSIGNED, process_rank, 1, MPI_COMM_WORLD);
				MPI_Send(&bucket_buffer[1], 1, MPI_UNSIGNED, process_rank, 2, MPI_COMM_WORLD);
				MPI_Send(buckets[buckets_sended++].vector, bucket_buffer[1], MPI_UNSIGNED, process_rank, 3, MPI_COMM_WORLD);
			}
			if (option == 2 && buckets_received < nbuckets) {	// Receive Buckets
				MPI_Recv(&bucket_buffer[0], 1, MPI_UNSIGNED, process_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&bucket_buffer[1], 1, MPI_UNSIGNED, process_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#if DEBUG
				printf("\nReceived Bucket vector %d (%d)", bucket_buffer[0], bucket_buffer[1]);
#endif
				MPI_Recv(buckets[bucket_buffer[0]].vector, bucket_buffer[1], MPI_UNSIGNED, process_rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				buckets[bucket_buffer[0]].size = bucket_buffer[1];
				++buckets_received;
			}
		}
	}
}

void receive_sort_send() {
	unsigned bucket_buffer[2];	// ID, Size
	unsigned option, terminate = 0;

	while (1) {
		if (terminate == 1) {
			break;
		} else {
			option = 1;	// Receive Buckets
			MPI_Send(&rank, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(&terminate, 1, MPI_UNSIGNED, 0, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (terminate) continue;
			MPI_Send(&option, 1, MPI_UNSIGNED, 0, 10, MPI_COMM_WORLD);

			MPI_Recv(&bucket_buffer[0], 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&bucket_buffer[1], 1, MPI_UNSIGNED, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			unsigned bucket_vector[bucket_buffer[1]];

			MPI_Recv(&bucket_vector, bucket_buffer[1], MPI_UNSIGNED, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			quicksort(bucket_vector, 0, bucket_buffer[1]-1);	// Sort vector

			option = 2;	// Send sorted Buckets
			MPI_Send(&rank, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(&terminate, 1, MPI_UNSIGNED, 0, 666, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Send(&option, 1, MPI_UNSIGNED, 0, 10, MPI_COMM_WORLD);

			MPI_Send(&bucket_buffer[0], 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);
			MPI_Send(&bucket_buffer[1], 1, MPI_UNSIGNED, 0, 2, MPI_COMM_WORLD);
			MPI_Send(bucket_vector, bucket_buffer[1], MPI_UNSIGNED, 0, 3, MPI_COMM_WORLD);
		}
	}
}

void join_buckets(Bucket buckets[], unsigned nbuckets, unsigned vector[], unsigned size) {
	unsigned i, j, bucket_size, vector_index = 0;

	// for (i = 0; i < nbuckets; ++i) printf("\nBUCKET SIZE = %d\n", buckets[i].size);

	for (i = 0; i < nbuckets; ++i) {
		bucket_size = buckets[i].size;
		if (bucket_size == 0) continue;
		for (j = 0; j < buckets[i].size; ++j) {
			if (buckets[i].size == 0) break;
			vector[vector_index] = buckets[i].vector[j];
			vector_index++;
		}

	}
}

void quicksort(unsigned vector[], int low, int high) {
	if (low < high) {
		int part = partition(vector, low, high);
		quicksort(vector, low, part-1);
		quicksort(vector, part+1, high);
	}
}

int partition(unsigned vector[], int low, int high) {
	int i, j;
	unsigned pivot = vector[high];
	i = low - 1;

	for (j = low; j < high; ++j) {
		if (vector[j] <= pivot) {
			++i;
			swap(&vector[i], &vector[j]);
		}
	}
	swap(&vector[i + 1], &vector[high]);
	return i + 1;
}

void swap(unsigned *a, unsigned *b) {
	unsigned aux = *a;
	*a = *b;
	*b = aux;
}
