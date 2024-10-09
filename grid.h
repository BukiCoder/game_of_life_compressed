#include <stdio.h>
#include <stdlib.h>
#include "bin_utils.h"

typedef struct grid
{
    int size_x;
    int size_y;
    int** state;
    int** new_state;
} GRID;

void saveFild(const char* filename, GRID* grid);
GRID* loadFild(const char* filename);
int freeFild(GRID* grid);

