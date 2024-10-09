int get(int** arr, int y, int x)
{
    return (arr[y][x / 32] >> (31 - (x % 32))) & 0x00000001;
}

int set(int** arr, int val, int y, int x)
{
    arr[y][x / 32] = (arr[y][x / 32] & (~((1) << (31 -  (x% 3 2)))))  + ((val) << (31 - (x % 32)));
}
