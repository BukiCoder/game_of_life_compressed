#include <stdio.h>
#include <stdlib.h>
#inlude "bin_utils.h"

typedef struct grid
{
    int size_x;
    int size_y;
    int** state;
    int** new_state;
} GRID;

void saveFild(const char* filename, GRID* grid)
{
    FILE* file = fopen(filename, "w");
    if(file == NULL)
    {
        fprintf(stderr, "%s", "Cannot open file");
        exit(1);
    }
    fprintf(file, "%d %d\n", grid->size_x, grid->size_y);
    for(int i = 0; i < fild->size_y; i++)
    {
        for(int j = 0; j < grid->size_x; j++)
        {
            fprintf(file, "%d", get(grid->state,i,j));
        }
        fprintf(file, "\n");
    }
}

GRID* loadFild(const char* filename)
{
    FILE* file = fopen(filename, "r");
    if(file != NULL)
    {
        GRID* fild = malloc(sizeof(*fild));

        int x, y = 0;
        fscanf(file, "%d %d", &y, &x);

        int** state     = calloc(y, sizeof(int*));
        int** new_state = calloc(y, sizeof(int*));
        fgetc(file);

        for(int i = 0; i < y; i++)
        {
            state[i] = calloc(x /32 + (x % 32 == 0? 0 : 1), sizeof(int));
            new_state[i] = calloc(x/32 + (x % 32 == 0? 0 : 1), sizeof(int));
            for(int j = 0; j < x; j++)
            {
                set(state, fgetc(file) - '0', i, j) ;
            }
            fgetc(file);
        }
        grid->size_x = x;
        grid->size_y = y;
        grid->state = state;
        grid->new_state = new_state;
        fclose(file);
        return grid;
    }
    else
    {
        return NULL;
    }
}

int freeFild(GRID* grid)
{
    free(grid->state);
    free(grid->new_state);
    free(grid);
    return 0;
}
