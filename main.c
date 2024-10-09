#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include "fild.h"
#include "binaryUtils.h"
#include "dirent.h"

#define ls 0b00100100100100100100100100100100
#define s2 0b01001001001001001001001001001001
#define s3 0b01101101101101101101101101101101

#define lsL 0b0010010010010010010010010010010010010010010010010010010010010010
#define s2L 0b0100100100100100100100100100100100100100100100100100100100100100
#define s3L 0b0110110110110110110110110110110110110110110110110110110110110110

void process_cell(FILD* fild, int i, int j)
{
	int iu = (i + 1) % fild->size_y;
	int id = (i - 1 + fild->size_y) % fild->size_y;

	int ju = (j + 1) % fild->size_x;
	int jd = (j - 1 + fild->size_x) % fild->size_x;

	int L = get(fild->state, iu, j);
	int R = get(fild->state, id, j);

	int LD = get(fild->state, i, ju);
	int D = get(fild->state, id, ju);
	int RD = get(fild->state, iu, ju);

	int LU = get(fild->state, i, jd);
	int U = get(fild->state, id, jd);
	int RU = get(fild->state, iu, jd);

	int sum = L + R + LD + D + RD + LU + U + RU;

	if (sum == 3 && get(fild->state, i, j) == 0)
	{
		set(fild->new_state, 1, i, j);
	}
	else if (get(fild->state, i, j) == 1 && !(sum == 2 || sum == 3))
	{
		set(fild->new_state, 0, i, j);
	}
	else
	{
		set(fild->new_state, get(fild->state, i, j), i, j);
	}
}

void next_step(FILD* fild)
{
	for (int i = 0; i < fild->size_y; i++)
	{
		for (int j = 0; j < fild->size_x; j++)
		{
			process_cell(fild, i, j);
		}
	}
	int** tmp = fild->state;
	(fild->state) = fild->new_state;
	fild->new_state = tmp;
}

long long getChL(unsigned long long block1, unsigned long long block2, unsigned long long block3)
{
	unsigned long long sum = 0;
	sum += (block1 & lsL) + ((block1 >> 1) & lsL) + ((block1 >> 2) & lsL); //sum higher
	sum += (block3 & lsL) + ((block3 >> 1) & lsL) + ((block3 >> 2) & lsL); //sum lower

	unsigned long long sumLR = (block2 & lsL) + ((block2 >> 2) & lsL); //sum left + right
	unsigned long long kostil = (sum ^ lsL) & (sumLR ^ (lsL | (lsL << 2))); //sumLR == 10(2) and sum == 110(6)
	sum += sumLR - (kostil & (kostil >> 1) & (kostil >> 2) & lsL); //if kostil, then sub 1

	unsigned long long sxs2 = (sum ^ ~s2L);
	unsigned long long sxs3 = (sum ^ ~s3L);
	unsigned long long is2 = sxs2 & (sxs2 >> 1) & (sxs2 >> 2) & lsL; //sum == 2
	unsigned long long is3 = sxs3 & (sxs3 >> 1) & (sxs3 >> 2) & lsL; //sum == 3

	return ((is3 & ((~(block2 >> 1)))) | ((~(is2 | is3)) & ((block2 >> 1)))) & lsL;
}

int getCh(int block1, int block2, int block3)
{
	unsigned int sum = 0;
	sum += (block1 & ls) + ((block1 >> 1) & ls) + ((block1 >> 2) & ls); //sum higher
	sum += (block3 & ls) + ((block3 >> 1) & ls) + ((block3 >> 2) & ls); //sum lower

	unsigned int sumLR = (block2 & ls) + ((block2 >> 2) & ls); //sum left + right
	int kostil = (sum ^ ls) & (sumLR ^ (ls | (ls << 2))); //sumLR == 10(2) and sum == 110(6)
	sum += sumLR - (kostil & (kostil >> 1) & (kostil >> 2) & ls); //if kostil, then sub 1

	unsigned int sxs2 = (sum ^ ~s2);
	unsigned int sxs3 = (sum ^ ~s3);
	unsigned int is2 = sxs2 & (sxs2 >> 1) & (sxs2 >> 2) & ls; //sum == 2
	unsigned int is3 = sxs3 & (sxs3 >> 1) & (sxs3 >> 2) & ls; //sum == 3

	return ((is3 & ((~(block2 >> 1)))) | ((~(is2 | is3)) & ((block2 >> 1)))) & ls;
}

inline void process_block(FILD* fild, int b1, int b2, int b3, int j, int start, int end)
{
	int block1 = fild->state[b1][j];
	int block2 = fild->state[b2][j];
	int block3 = fild->state[b3][j];

	int a1 = getCh(block1, block2, block3) << 1;
	int a2 = getCh(block1 << 1, block2 << 1, block3 << 1);
	int a3 = getCh(block1 << 2, block2 << 2, block3 << 2) >> 1;

	int change = a1 ^ a2 ^ a3;

	fild->new_state[b2][j] = fild->state[b2][j] ^ change;

	process_cell(fild, b2, start);
	process_cell(fild, b2, end);
}

inline void process_blockL(FILD* fild, int b1, int b2, int b3, int j, int start, int end)
{
	int kk = ((fild->state[b1]))[j];
	long long block1 = ((long long*)(fild->state[b1]))[j];
	long long block2 = ((long long*)fild->state[b2])[j];
	long long block3 = ((long long*)fild->state[b3])[j];

	long long a1 = getChL(block1, block2, block3) << 1;
	long long a2 = getChL(block1 << 1, block2 << 1, block3 << 1);
	long long a3 = getChL(block1 << 2, block2 << 2, block3 << 2) >> 1;

	long long change = a1 ^ a2 ^ a3;

	((long long*)fild->new_state[b2])[j] = ((long long*)fild->state[b2])[j] ^ change;

	process_cell(fild, b2, start);
	process_cell(fild, b2, end);
}

void nextStepL(FILD* fild)
{
	for (int j = 0; j < fild->size_x / 64; j++)
	{
		process_blockL(fild, fild->size_y - 1, 0, 1, j, j * 64, j * 64 + 63);
	}

	for (int i = 1; i < fild->size_y - 1; i++)
	{
		for (int j = 0; j < fild->size_x / 64; j++)
		{
			process_blockL(fild, i - 1, i, i + 1, j, j * 64, j * 64 + 63);
		}
	}

	for (int j = 0; j < fild->size_x / 64; j++)
	{
		process_blockL(fild, fild->size_y - 2, fild->size_y - 1, 0, j, j * 64, j * 64 + 63);
	}

	if ((fild->size_x) % 64 != 0)
	{
		int j = fild->size_x / 64;
		process_blockL(fild, fild->size_y - 1, 0, 1, j, j * 64, fild->size_x - 1);
		for (int i = 1; i < fild->size_y - 1; i++)
		{
			process_blockL(fild, i - 1, i, i + 1, j, j * 64, fild->size_x - 1);
		}
		process_blockL(fild, fild->size_y - 2, fild->size_y - 1, 0, j, j * 64, fild->size_x - 1);
	}

	int** tmp = fild->state;
	(fild->state) = fild->new_state;
	fild->new_state = tmp;
}

void nextStepInt(FILD* fild)
{
	for (int j = 0; j < fild->size_x / 32; j++)
	{
		process_block(fild, fild->size_y - 1, 0, 1, j, j * 32, j * 32 + 31);
	}

	for (int i = 1; i < fild->size_y - 1; i++)
	{
		for (int j = 0; j < fild->size_x / 32; j++)
		{
			process_block(fild, i - 1, i, i + 1, j, j * 32, j * 32 + 31);
		}
	}

	for (int j = 0; j < fild->size_x / 32; j++)
	{
		process_block(fild, fild->size_y - 2, fild->size_y - 1, 0, j, j * 32, j * 32 + 31);
	}

	if ((fild->size_x) % 32 != 0)
	{
		int j = fild->size_x / 32;
		process_block(fild, fild->size_y - 1, 0, 1, j, j * 32, fild->size_x - 1);
		for (int i = 1; i < fild->size_y - 1; i++)
		{
			process_block(fild, i - 1, i, i + 1, j, j * 32, fild->size_x - 1);
		}
		process_block(fild, fild->size_y - 2, fild->size_y - 1, 0, j, j * 32, fild->size_x - 1);
	}

	int** tmp = fild->state;
	(fild->state) = fild->new_state;
	fild->new_state = tmp;
}

int test(FILD* fild1, FILD* fild2, void (*step1) (FILD*), void (*step2) (FILD*), int rounds_count, int align)
{
	for (int i = 0; i < rounds_count; i++)
	{
		saveFild("f1.txt", fild1);
		step1(fild1);
		saveFild("f1s.txt", fild1);
		saveFild("f2.txt", fild2);
		step2(fild2);
		saveFild("f2s.txt", fild2);
		for (int j = 0; j < fild1->size_y; j++)
		{
			if (memcmp(fild1->state[j], fild2->state[j], (fild1->size_x / align * 8)) != 0)
			{
				return 0;
			}
			if (fild1->size_x % align != 0 && ((fild1->state[j][fild1->size_x / align - 1] & ((~((long long)0)) << (align - 1 - fild1->size_x % align)))
				!= (fild2->state[j][fild1->size_x / align - 1] & ((~((long long)0)) << (align - 1 - fild1->size_x % align)))))
			{
				return 0;
			}
		}
	}
	return 1;
}
