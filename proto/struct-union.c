/*
 * @brief   Example of custom data type with union.
 *
 * @author  Maksims Abalenkovs
 * @email   maksims.abalenkovs@stfc.ac.uk
 * @date    Oct 17, 2024
 * @version 0.2
 */

#include <stdio.h>
#include <stdlib.h>

// Matrix block
typedef struct {
    union {
        int    *i;  
        double *d;
    } a;     // pointer to first element in matrix block
    int  k;  // block index
    int  p;  // no. of block rows
    int  q;  // no. of block columns
    int  m;  // no. of matrix rows
    int  n;  // no. of matrix columns
    int  r;  // MPI rank of block owner
} mtrx_blk;

int main(int argc, char *argv[]) {

    // no. of matrix rows
    int const M = 10;

    // no. of matrix columns
    int const N = 10;

    // no. of blocks
    int const NBLK = 2;

    // no. of block rows
    int const P = M/NBLK;

    // no. of block columns
    int const Q = N;

    // allocate memory for matrix of doubles
    double *A = (double*) malloc(M*N*sizeof(double));

    // allocate memory for matrix of integers
    int    *B = (int*)    malloc(M*N*sizeof(int));

    // initialise matrices with consecutive natural numbers
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {

            A[i*N+j] = (double) i*N+j;
            B[i*N+j] = i*N+j;
        }
    }

    // @test print out matrix of doubles to screen
    printf("A(%d x %d):\n", M, N);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {

            printf("%f ", A[i*N+j]);
        }
        printf("\n");
    }

    // @test print out matrix of integers to screen
    printf("B(%d x %d):\n", M, N);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {

            printf("%d ", B[i*N+j]);
        }
        printf("\n");
    }

    // allocate array of blocks of doubles
    mtrx_blk *blk_d = (mtrx_blk*) malloc(NBLK*sizeof(mtrx_blk));

    // allocate array of blocks of integers
    mtrx_blk *blk_i = (mtrx_blk*) malloc(NBLK*sizeof(mtrx_blk));

    // for each block of doubles AND
    // for each block of integers
    for (int k = 0; k < NBLK; k++) {

        blk_d[k] = (mtrx_blk) { .a.d = &A[k*P*Q], .k = k, .p = P, .q = Q, .m = M, .n = N, .r = k };
        blk_i[k] = (mtrx_blk) { .a.i = &B[k*P*Q], .k = k, .p = P, .q = Q, .m = M, .n = N, .r = k };
    }

    // @test print (double) blocks out to screen
    // @todo STOPPED HERE!!!
    for (int k = 0; k < NBLK; k++) {

        printf("\na_%d(%d x %d):\n", k, blk[k].p, blk[k].q);

        for (int i = 0; i < blk[k].p; i++) {
            for (int j = 0; j < blk[k].q; j++) {

                printf("%f ", blk[k].a.d[i*N+j]);
            }
            printf("\n");
        }
    }

    // free memory of arrays of blocks
    free(blk_d);
    free(blk_i);

    // free memory of matrix
    free(A);

    return 0;
}

// @eof struct-union.c
