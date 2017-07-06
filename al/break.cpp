#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <queue>
#include <sys/time.h>
#include <vector>
using namespace std;
#define fst first
#define snd second
#define mp make_pair
typedef long long ll;
typedef unsigned long long ull;
typedef unsigned int uint;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;

const int NB_SBOXES = 8;

const int E[48] = {
	32, 1, 2, 3, 4, 5,
	4, 5, 6, 7, 8, 9,
	8, 9, 10, 11, 12, 13,
	12, 13, 14, 15, 16, 17,
	16, 17, 18, 19, 20, 21,
	20, 21, 22, 23, 24, 25,
	24, 25, 26, 27, 28, 29,
	28, 29, 30, 31, 32, 1
};

const int P[32] = {
	16, 7, 20, 21, 29, 12, 28, 17,
	1, 15, 23, 26, 5, 18, 31, 10,
	2, 8, 24, 14, 32, 27, 3, 9,
	19, 13, 30, 6, 22, 11, 4, 25
};

const int S[8][4][16] = {
	{ 
		{ 14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7 },
		{ 0, 15, 7, 4, 14, 2, 13, 1, 10, 6, 12, 11, 9, 5, 3, 8 },
		{ 4, 1, 14, 8, 13, 6, 2, 11, 15, 12, 9, 7, 3, 10, 5, 0 },
		{ 15, 12, 8, 2, 4, 9, 1, 7, 5, 11, 3, 14, 10, 0, 6, 13 }
	},
	{
		{ 15, 1, 8, 14, 6, 11, 3, 4, 9, 7, 2, 13, 12, 0, 5, 10 },
		{ 3, 13, 4, 7, 15, 2, 8, 14, 12, 0, 1, 10, 6, 9, 11, 5 },
		{ 0, 14, 7, 11, 10, 4, 13, 1, 5, 8, 12, 6, 9, 3, 2, 15 },
		{ 13, 8, 10, 1, 3, 15, 4, 2, 11, 6, 7, 12, 0, 5, 14, 9 }
	},
	{
		{ 10, 0, 9, 14, 6, 3, 15, 5, 1, 13, 12, 7, 11, 4, 2, 8 },
		{ 13, 7, 0, 9, 3, 4, 6, 10, 2, 8, 5, 14, 12, 11, 15, 1 },
		{ 13, 6, 4, 9, 8, 15, 3, 0, 11, 1, 2, 12, 5, 10, 14, 7 },
		{ 1, 10, 13, 0, 6, 9, 8, 7, 4, 15, 14, 3, 11, 5, 2, 12 }
	},
	{
		{ 7, 13, 14, 3, 0, 6, 9, 10, 1, 2, 8, 5, 11, 12, 4, 15 },
		{ 13, 8, 11, 5, 6, 15, 0, 3, 4, 7, 2, 12, 1, 10, 14, 9 },
		{ 10, 6, 9, 0, 12, 11, 7, 13, 15, 1, 3, 14, 5, 2, 8, 4 },
		{ 3, 15, 0, 6, 10, 1, 13, 8, 9, 4, 5, 11, 12, 7, 2, 14 }
	},
	{
		{ 2, 12, 4, 1, 7, 10, 11, 6, 8, 5, 3, 15, 13, 0, 14, 9 },
		{ 14, 11, 2, 12, 4, 7, 13, 1, 5, 0, 15, 10, 3, 9, 8, 6 },
		{ 4, 2, 1, 11, 10, 13, 7, 8, 15, 9, 12, 5, 6, 3, 0, 14 },
		{ 11, 8, 12, 7, 1, 14, 2, 13, 6, 15, 0, 9, 10, 4, 5, 3 }
	},
	{
		{ 12, 1, 10, 15, 9, 2, 6, 8, 0, 13, 3, 4, 14, 7, 5, 11 },
		{ 10, 15, 4, 2, 7, 12, 9, 5, 6, 1, 13, 14, 0, 11, 3, 8 },
		{ 9, 14, 15, 5, 2, 8, 12, 3, 7, 0, 4, 10, 1, 13, 11, 6 },
		{ 4, 3, 2, 12, 9, 5, 15, 10, 11, 14, 1, 7, 6, 0, 8, 13 }
	},
	{
		{ 4, 11, 2, 14, 15, 0, 8, 13, 3, 12, 9, 7, 5, 10, 6, 1 },
		{ 13, 0, 11, 7, 4, 9, 1, 10, 14, 3, 5, 12, 2, 15, 8, 6 },
		{ 1, 4, 11, 13, 12, 3, 7, 14, 10, 15, 6, 8, 0, 5, 9, 2 },
		{ 6, 11, 13, 8, 1, 4, 10, 7, 9, 5, 0, 15, 14, 2, 3, 12 }
	},
	{
		{ 13, 2, 8, 4, 6, 15, 11, 1, 10, 9, 3, 14, 5, 0, 12, 7 },
		{ 1, 15, 13, 8, 10, 3, 7, 4, 12, 5, 6, 11, 0, 14, 9, 2 },
		{ 7, 11, 4, 1, 9, 12, 14, 2, 0, 6, 10, 13, 15, 3, 5, 8 },
		{ 2, 1, 14, 7, 4, 10, 8, 13, 15, 12, 9, 0, 3, 5, 6, 11 }
	}
};

const int rotations[16] = {
	1, 1, 2, 2, 2, 2, 2, 2,
	1, 2, 2, 2, 2, 2, 2, 1
};

const int PC1[56] = {
	57, 49, 41, 33, 25, 17, 9,
	1, 58, 50, 42, 34, 26, 18,
	10, 2, 59, 51, 43, 35, 27,
	19, 11, 3, 60, 52, 44, 36,
	63, 55, 47, 39, 31, 23, 15,
	7, 62, 54, 46, 38, 30, 22,
	14, 6, 61, 53, 45, 37, 29,
	21, 13, 5, 28, 20, 12, 4
};

const int PC2[48] = {
	14, 17, 11, 24, 1, 5,
	3, 28, 15, 6, 21, 10,
	23, 19, 12, 4, 26, 8,
	16, 7, 27, 20, 13, 2,
	41, 52, 31, 37, 47, 55,
	30, 40, 51, 45, 33, 48,
	44, 49, 39, 56, 34, 53,
	46, 42, 50, 36, 29, 32
};

const int IP[64] = {
	58, 50, 42, 34, 26, 18, 10, 2,
	60, 52, 44, 36, 28, 20, 12, 4,
	62, 54, 46, 38, 30, 22, 14, 6,
	64, 56, 48, 40, 32, 24, 16, 8,
	57, 49, 41, 33, 25, 17, 9, 1,
	59, 51, 43, 35, 27, 19, 11, 3,
	61, 53, 45, 37, 29, 21, 13, 5,
	63, 55, 47, 39, 31, 23, 15, 7
};

const int FP[64] = {
	40, 8, 48, 16, 56, 24, 64, 32,
	39, 7, 47, 15, 55, 23, 63, 31,
	38, 6, 46, 14, 54, 22, 62, 30,
	37, 5, 45, 13, 53, 21, 61, 29,
	36, 4, 44, 12, 52, 20, 60, 28,
	35, 3, 43, 11, 51, 19, 59, 27,
	34, 2, 42, 10, 50, 18, 58, 26,
	33, 1, 41, 9, 49, 17, 57, 25
};

ll cut(ll x, ll n) {
	return x & ((1ll << n) - 1);
}

void print_bits(ll x, int n) {
	for (int i = n - 1; i >= 0; --i) {
		printf("%lld", (x >> i) & 1);
		if ((n - i) % 4 == 0)
			printf(" ");
	}
}

vector<int> get_com_pos(int sbox1, int sbox2) {
//	assert(abs(min(sbox1 - sbox2, min(8 - sbox1 + sbox2, 8 - sbox2 + sbox1))) == 1);
	vector<int> ans;
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			if (E[6 * sbox1 + i] == E[6 * sbox2 + j])
				ans.push_back(6 - i - 1);
	return ans;
}

int merge_bits(vector<int> com_pos, int rem_val, int com_val, int nb_bits) {
	int p = 1;
	int ans = 0;
	for (int i = 0; i < nb_bits; ++i) {
		bool rem = find(com_pos.begin(), com_pos.end(), i) == com_pos.end(); 
		int& val = rem ? rem_val : com_val;
		ans += p * (val % 2);
		val /= 2;
		p *= 2;
	}
	return ans;
}

int do_sbox(int sbox, int in) {
	return S[sbox][(in & 1) | ((in >> 4) & 2)][(in >> 1) & ((1 << 4) - 1)];
}

ull shiftrot(ull x, ull shift, int nb_bits) {
		return cut(x << shift, nb_bits) | ((x >> (nb_bits - shift)) & cut(~0ull, shift));
}

int xor_all(int x) {
	int ans = 0;
	while (x > 0) {
		ans ^= (x % 2);
		x /= 2;
	}
	return ans;
}

ull do_perm(ull in, const int *a, int nin, int nout) {
	ull out = 0;
	for (int i = 0; i < nout; ++i)
		out = (out << 1) | ((in >> (nin - a[i])) & 1);
	return out;
}

ull encrypt_one(ull in, ull key) {
	ull l_in = in >> 32;
	ull r_in = cut(in, 32);
	ull in_sbox = do_perm(r_in, E, 32, 48) ^ key;
	ull out = 0;
	for (int i = 0; i < 8; ++i)
		out = (out << 4) | do_sbox(i, cut(in_sbox >> (6 * (7 - i)), 6));
	out = do_perm(out, P, 32, 32);
	return (r_in << 32) | (l_in ^ out);
}

// without initial&final permutations 
ull encrypt(ull in, ull key, int levelmax, int level = 1) {
	if (level > levelmax)
		return in;
	ull lkey = shiftrot(key >> 28, rotations[level - 1], 28);
	ull rkey = shiftrot(cut(key, 28), rotations[level - 1], 28);
	ull new_key = (lkey << 28) | rkey;
	ull out = encrypt_one(in, do_perm(new_key, PC2, 56, 48));
	return encrypt(out, new_key, levelmax, level + 1);
}

void test() {
	ull in = 0x0123456789ABCDEFull;
	ull key = in;
//	ull key = 0x133457799BBCDFF1ull;
	in = do_perm(in, IP, 64, 64);
	key = do_perm(key, PC1, 64, 56);
	ull out = encrypt(in, key, 16);
	out = do_perm((cut(out, 32) << 32) | (out >> 32), FP, 64, 64);
	printf("out = 0x%llX\n", out);
}

// number of plaintext/ciphertext couples
const int N = 1 << 10;

void extract() {
	FILE *fp = fopen("approx", "r");
	assert(fp != NULL);
	char line[1024];
	while (fgets(line, sizeof line, fp) != NULL) {
		char *p = line;
		double score;
		int sz = 0;
		sscanf(line, "%f%n", &score, sz);
		vector<ull> mask(16);
		for (int i = 0; i < 16; ++i) {
			p += sz;
			sscanf(line, "%llu%n", &mask[i], &sz);
		}
	}
	fclose(fp);
}

int main() {
	extract();	
	return 0;
}
