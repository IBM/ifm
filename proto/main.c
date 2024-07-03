/*
 * @brief   Prototypes loops for IFM's diffusive routing skipping last row and
 *          skipping last column.
 *
 * @author  Maksims Abalenkovs
 * @email   maksims.abalenkovs@stfc.ac.uk
 * @date    Jul 3, 2024
 * @version 0.1
 */

#include <stdio.h>

int main() {

    int const _nrow = 6;
    int const _ncol = 4;

    int cur = 0;

    printf("Entering loop 1 (Skip last row)\n");

    for (int i = 0; i < _nrow-1; i++) {
        for (int j = 0; j < _ncol; j++) {

            cur = i*(_ncol)+j;
            printf("%d ", cur);
        }
    }

    printf("\n\nEntering loop 2 (Skip last column)\n");

    for (int i = 0; i < _nrow; i++) {
        for (int j = 0; j < _ncol-1; j++) {

            cur = i*(_ncol)+j;
            printf("%d ", cur);
        }
    }
}

// @eof main.c
