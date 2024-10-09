#include <stdio.h>
typedef struct fild
{
    int size_x;
    int size_y;
    int** state;
    int** new_state;
} FILD;
int get(int** arr, int y, int x)
{
    return (arr[y][x / 32] >> (31 - (x % 32))) & 0x00000001;
}
void saveFild(const char* filename, FILD* fild)
{
    FILE* file = fopen(filename, "w");
    fprintf(file, "%d %d\n", fild->size_x, fild->size_y);
    for(int i = 0; i < fild->size_y; i++)
    {
        for(int j = 0; j < fild->size_x; j++)
        {
            fprintf(file, "%d", get(fild->state,i,j));
        }
        fprintf(file, "\n");
    }
}
FILD* loadFild(const char* filename)
{
    FILE* file = fopen(filename, "r");
    if(file != 0)
    {
        FILD* fild = malloc(sizeof(FILD));

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
                state[i][j / 32] = (state[i][j / 32] & (~((1) << (31 -  (j%32)))))  + ((fgetc(file) - '0') << (31 - (j % 32)));
               // set(state, fgetc(file) - '0', i, j) ;
            }
            fgetc(file);
        }
            fild->size_x = x;
            fild->size_y = y;
            fild->state = state;
            fild->new_state = new_state;
            fclose(file);
            return fild;
        }
        else
        {
            return 0;
        }



    }
    int freeFild(FILD* fild)
    {
        free(fild->state);
        free(fild->new_state);
        free(fild);
        return 0;
    }

