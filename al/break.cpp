#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <queue>
#include <set>
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

ull shiftrot_left_1(ull x, int nb_bits) {
	return cut(x << 1, nb_bits) | ((x >> (nb_bits - 1)) & 1); 
}

ull shiftrot_left(ull x, int shift, int nb_bits) {
	ull ans = x;
	for (int i = 0; i < shift; ++i)
		ans = shiftrot_left_1(ans, nb_bits);
	return ans;
}

ull shiftrot_right_1(ull x, int nb_bits) {
	return (x >> 1) | ((x & 1) << (nb_bits - 1));
}

ull shiftrot_right(ull x, int shift, int nb_bits) {
	ull ans = x;
	for (int i = 0; i < shift; ++i)
		ans = shiftrot_right_1(ans, nb_bits);
	return ans;
}

int xor_all(int x) {
	int ans = 0;
	while (x > 0) {
		ans ^= (x % 2);
		x /= 2;
	}
	return ans;
}

ull nth_bit(ull x, int i) {
	return (x >> i) & 1;
}

ull do_perm(ull in, const int *a, int nin, int nout) {
	ull out = 0;
	for (int i = 0; i < nout; ++i)
		out = (out << 1) | ((in >> (nin - a[i])) & 1);
	return out;
}

ull do_inv_perm(ull in, const int *a, int nin, int nout) {
	ull out = 0;
	for (int i = 0; i < nin; ++i)
		out |= nth_bit(in, i) << (nout - a[nin - i - 1]);
	return out;
}

ull F_func(ull in, ull key) {
	ull ans = 0;
	ull in_sbox = do_perm(in, E, 32, 48) ^ do_perm(key, PC2, 56, 48);
	for (int i = 0; i < 8; ++i)
		ans = (ans << 4) | do_sbox(i, cut(in_sbox >> (6 * (7 - i)), 6));
	return do_perm(ans, P, 32, 32);
}

ull encrypt_one(ull in, ull key) {
	ull l_in = in >> 32;
	ull r_in = cut(in, 32);
	return (r_in << 32) | (l_in ^ F_func(r_in, key));
}

ull swap_halves(ull in) {
	return (cut(in, 32) << 32) | (in >> 32);
}

ull rotkey_left(ull key, int shift) {
	return (shiftrot_left(key >> 28, shift, 28) << 28) | shiftrot_left(cut(key, 28), shift, 28);
}

ull rotkey_right(ull key, int shift) {
	return (shiftrot_right(key >> 28, shift, 28) << 28) | shiftrot_right(cut(key, 28), shift, 28);
}

// without initial&final permutations 
ull encrypt(ull in, ull key, int levelmax, int level = 1) {
	if (level > levelmax)
		return in;
	ull new_key = rotkey_left(key, rotations[level - 1]);
	ull out = encrypt_one(in, new_key);
	return encrypt(out, new_key, levelmax, level + 1);
}

void test() {
	ull in = 0x0DFEBD56523FAA85ull;
	ull key = in;
	in = do_perm(in, IP, 64, 64);
	key = do_perm(key, PC1, 64, 56);
	ull out = encrypt(in, key, 16);
	out = do_perm(swap_halves(out), FP, 64, 64);
	printf("out = 0x%llX\n", out);
	/*
	ull tmp = encrypt(in, key, 15);
	print_bits(tmp, 64);
	puts("");
	ull out = encrypt(in, key, 16);
//	print_bits(out, 64);
	out = (out >> 32) | (cut(out, 32) << 32);
//	print_bits(out, 64);
	ull lkey = key >> 28, rkey = cut(key, 28);
	out = encrypt(out, (shiftrot_right_1(lkey, 28) << 28) | shiftrot_right_1(rkey, 28), 1);
	print_bits((out >> 32) | (cut(out, 32) << 32), 64);
	puts("");
	*/
}

// number of plaintext/ciphertext couples
const ull N = 1 << 22;
const int L = 12;
const int M = 10;
const int LEVEL = 6;
//vector<pair<double, <ull>>> approx[8];

ull one(int nb_bits) {
	return (1ull << nb_bits) - 1;
}

int popcount(ull x) {
	int ans = 0;
	while (x > 0) {
		ans += x & 1;
		x >>= 1;
	}
	return ans;
}

vector<pair<double, vector<ull>>> approxs;
map<pair<ull, ull>, vector<int>> mmask;

void extract() {
	FILE *fp = fopen("approx6", "r");
	assert(fp != NULL);
	char line[1024];
	map<vector<ull>, double> mapprox;
	int nblines = 0;
	while (fgets(line, sizeof line, fp) != NULL) {
		char *p = line;
		double score;
		int sz = 0;
		sscanf(line, "%le%n", &score, &sz);
		vector<ull> mask(LEVEL);
		for (int i = 0; i < LEVEL; ++i) {
			p += sz;
			sscanf(p, "%llu%n", &mask[i], &sz);
		}
		mapprox[mask] = score;
		nblines++;
	}
	printf("%d lines\n", nblines);
	for (auto it: mapprox) {
		vector<ull> mask = it.fst;
		double score = it.snd;
		uint mask_pl = cut(mask[0], 32) ^ (mask[1] >> 32);
//		uint mask_ph = mask[0] >> 32;
		uint mask_cl = cut(mask[LEVEL - 1], 32) ^ (mask[LEVEL - 2] >> 32);
//		uint mask_ch = mask[13] >> 32;
		ull key_mask_in = 0, key_mask_out = 0;

		mask_pl = do_inv_perm(mask_pl, P, 32, 32);
		for (int i = 7; i >= 0; --i) {
			ull lmask = cut(mask_pl >> (4 * i), 4) != 0 ? one(6) : 0;
			key_mask_in = (key_mask_in << 6) | lmask;
		}
		key_mask_in = do_inv_perm(key_mask_in, PC2, 48, 56);
		key_mask_in = rotkey_right(key_mask_in, 1);

		mask_cl = do_inv_perm(mask_cl, P, 32, 32);
		for (int i = 7; i >= 0; --i) {
			ull lmask = cut(mask_cl >> (4 * i), 4) != 0 ? one(6) : 0;
			key_mask_out = (key_mask_out << 6) | lmask;
		}
		key_mask_out = do_inv_perm(key_mask_out, PC2, 48, 56);
		assert(rotkey_left(rotkey_right(key_mask_out, 14), 14) == key_mask_out);
		key_mask_out = rotkey_right(key_mask_out, 14); // delete for LEVEL == 14
//		ull key_mask = key_mask_in | key_mask_out;
//		if (popcount(key_mask) <= L) {
		if (fabs(score - 3.81e-3) <= 1e-4) {
			approxs.push_back({score, mask});
			auto key_pair = make_pair(key_mask_in, key_mask_out);
			if (mmask.find(key_pair) == mmask.end())
				mmask[key_pair] = vector<int>(0);
			mmask[key_pair].push_back(int(approxs.size()) - 1);
			/*
			printf("%e ", score);
			puts("");
			print_bits(key_mask_in, 56);
			puts("");
			print_bits(key_mask_out, 56);
			puts("");
			print_bits(key_mask, 56);
			puts("");
			puts("-----------------------");
			*/
		}
	}
	printf("size = %d\n", int(mmask.size()));
	for (auto it: mmask) {
		print_bits(it.fst.fst | it.fst.snd, 56);
		sort(it.snd.begin(), it.snd.end(), [&](int x, int y) { return approxs[x].fst > approxs[y].fst; });
		/*
		unique(it.snd.begin(), it.snd.end(), [&](int x, int y) {
			return approxs[x].snd[0] == approxs[y].snd[0] &&
				   approxs[x].snd[1] == approxs[y].snd[1];
		});
		*/
		// std::unique ?
		printf("-> %d %e\n", int(it.snd.size()), approxs[it.snd[0]].fst);
		/*
		for (int x: it.snd) {
			print_bits(approxs[x].snd[0], 64);
			print_bits(approxs[x].snd[1], 64); puts("");
		}
		*/
	}
	fclose(fp);
}

double square(double x) { return x * x; }

//ull X[N], Y[N];
ull key;

ull rand64() {
	return rand() + (1ull << 32) * rand();
}
/*
void generate_pairs() {
	for (int i = 0; i < N; ++i) {
		X[i] = rand64();
		Y[i] = encrypt(X[i], key, LEVEL + 2);
	}
}
*/

ull fill_with_mask(ull mask, ull val, int nb_bits) {
	ull ans = 0;
	for (int i = 0; i < nb_bits; ++i)
		if (nth_bit(mask, i)) {
			ans |= (val & 1) << i;
			val >>= 1;
		}
	return ans;
}

ull empty_with_mask(ull mask, ull val, int nb_bits) {
	ull ans = 0;
	for (int i = nb_bits - 1; i >= 0; --i)
		if (nth_bit(mask, i))
			ans = (ans << 1) | nth_bit(val, i); 
	return ans;
}

int F[1 << L][1 << M];

void break_cipher() {
	for (auto it: mmask) {
		auto key_mask = it.fst.fst | it.fst.snd;
		// on-line phase of Matsui's alg 2 in multiple dimension
		for (int i = 0; i < (1 << L); ++i)
			fill(F[i], F[i] + (1 << M), 0);
		int m = min(M, int(it.snd.size()));
		vector<uint> mask_pl(m), mask_ph(m), mask_cl(m), mask_ch(m);
		for (int j = 0; j < m; ++j) {
			mask_pl[j] = cut(approxs[it.snd[j]].snd[0], 32) ^ (approxs[it.snd[j]].snd[1] >> 32);
			mask_ph[j] = approxs[it.snd[j]].snd[0] >> 32;
			mask_cl[j] = cut(approxs[it.snd[j]].snd[LEVEL - 1], 32) ^ (approxs[it.snd[j]].snd[LEVEL - 2] >> 32);
			mask_ch[j] = approxs[it.snd[j]].snd[LEVEL - 1] >> 32;
		}
		ull real_k = empty_with_mask(key_mask, key, 56);
		ull rand_k = rand() % (1ull << popcount(key_mask));
		printf("%llu %llu %d\n", N, 1ull << popcount(key_mask), m);
		map<vector<uint>, int> pairtr;
		for (ull i = 0; i < N; ++i) {
			ull X = rand64();
			ull Y = swap_halves(encrypt(X, key, LEVEL + 2));
			vector<uint> v;
			ull mph = 0, mpl = 0, mch = 0, mcl = 0;
			for (int j = 0; j < m; ++j) {
				mph |= mask_ph[j];
				mpl |= mask_pl[j];
				mch |= mask_ch[j];
				mcl |= mask_cl[j];
			}
			v.push_back(cut(X, 32) & mph);
			v.push_back((X >> 32) & mpl);
			v.push_back((Y >> 32) & mcl);
			v.push_back(cut(Y, 32) & mch);
			mpl = do_inv_perm(mpl, P, 32, 32);
			ull mx = 0;
			for (int i = 7; i >= 0; --i) {
				ull lmask = cut(mpl >> (4 * i), 4) != 0 ? one(6) : 0;
				mx = (mx << 6) | lmask;
			}
			v.push_back(do_inv_perm(mx, E, 48, 32) & cut(X, 32)); // prendre le bon bit ?
			ull my = 0;
			mcl = do_inv_perm(mcl, P, 32, 32);
			for (int i = 7; i >= 0; --i) {
				ull lmask = cut(mcl >> (4 * i), 4) != 0 ? one(6) : 0;
				my = (my << 6) | lmask;
			}
			v.push_back(do_inv_perm(my, E, 48, 32) & cut(Y, 32)); // pareil
			if (pairtr.find(v) == pairtr.end()) pairtr[v] = 0;
			pairtr[v]++;
//			for (ull _k = 0; _k < (1ull << popcount(key_mask)); ++_k) {
			for (int _it = 0; _it < 2; ++_it)  {
				ull _k = _it == 0 ? real_k : rand_k;
//				ull k = _it == 0 ? key : fill_with_mask(key_mask, _k, 56);
				ull k = fill_with_mask(key_mask, _k, 56);
				ull x = F_func(cut(X, 32), rotkey_left(k, 1));
				ull y = F_func(cut(Y, 32), rotkey_left(k, 14)); // change when LEVEL = 14
				ull xp = F_func(v[4], rotkey_left(k, 1));
				ull yp = F_func(v[5], rotkey_left(k, 14));
//				assert(encrypt((cut(X, 32) << 32) | (x ^ (X >> 32)), rotkey_left(k, 1), LEVEL + 2, 2) == swap_halves(Y));
				/*
				assert(encrypt(encrypt(X, key, 1), rotkey_left(key, 1), 4, 2) == encrypt(X, key, 4));
				assert((((x ^ (X >> 32))) | (cut(X, 32) << 32)) == encrypt(X, key, 1));
				assert(encrypt(encrypt(X, key, 1), rotkey_left(key, 1), LEVEL + 2, 2) == encrypt(X, key, LEVEL + 2));
				assert(encrypt(X, key, LEVEL + 1) == ((cut(Y, 32) << 32) | (y ^ (Y >> 32))));
				assert(swap_halves(encrypt((x ^ (X >> 32)) | (cut(X, 32) << 32), rotkey_left(key, 1), LEVEL + 1, 2)) == 
					((y ^ (Y >> 32)) | (cut(Y, 32) << 32)));
				*/
				int eta = 0;
				for (int j = 0; j < m; ++j) {
//					print_bits(it.fst, 56); puts("");
//					print_bits(key, 56); puts("");
//					print_bits(k, 56); puts("");
//					assert((mask_pl & (x ^ (X >> 32))) == (mask_pl & cut(encrypt(X, key, 1), 32)));
					/*
					print_bits(key, 56); puts("");
					print_bits(k, 56); puts("");
					print_bits(it.fst, 56); puts("");
					print_bits(mask_ch & (y ^ (Y >> 32)), 64); puts("");
					print_bits(encrypt(Y, k, 1), 64); puts("");
					print_bits(encrypt(Y, rotkey_left(k, 1), 1), 64); puts("");
					print_bits(mask_ch & cut(encrypt(Y, rotkey_left(k, 11), 1), 32), 64); puts("");
					puts("-------");
					*/
//					assert((mask_cl & encrypt(0, rotkey_left(key, 13), 1)) == (mask_cl & encrypt(0, rotkey_left(k, 13), 1)));
//					assert((mask_cl & (y ^ (Y >> 32))) == (mask_cl & encrypt(Y, rotkey_left(k, 13), 1)));
					assert((mask_pl[j] & x) == (mask_pl[j] & xp));
					assert((mask_pl[j] & (X >> 32)) == (mask_pl[j] & v[1]));
					assert((mask_ph[j] & cut(X, 32)) == (mask_ph[j] & v[0]));
					assert((mask_cl[j] & y) == (mask_cl[j] & yp));
					assert((mask_cl[j] & (Y >> 32)) == (mask_cl[j] & v[2]));
					assert((mask_ch[j] & cut(Y, 32)) == (mask_ch[j] & v[3]));
					eta = (eta << 1) | xor_all(
						(mask_pl[j] & (x ^ (X >> 32))) ^
						(mask_ph[j] & cut(X, 32)) ^
						(mask_cl[j] & (y ^ (Y >> 32))) ^
						(mask_ch[j] & cut(Y, 32))
					);
					/*
					eta |= xor_all(
						(mask_pl & cut(X, 32)) ^
						(mask_ph & (X >> 32)) ^
						(mask_cl & cut(Y, 32)) ^
						(mask_ch & (Y >> 32))
					) << j;
					*/
				}
				F[_k][eta]++;
			}
			if (i % 1000 == 0) {
				printf("%d %.10e\n", F[real_k][0], 1. * F[real_k][0] / (i + 1) - 1. / (1 << m));
				printf("%d %.10e\n", F[rand_k][0], 1. * F[rand_k][0] / (i + 1) - 1. / (1 << m));
				puts("");
			}
		}
		// off-line phase of alg 2 using khi^2-method
		vector<pair<double, ull>> S;
		for (ull k = 0; k < (1ull << popcount(key_mask)); ++k) {
			double s = 0;
			for (int eta = 0; eta < (1 << m); ++eta)
				s += square(1. * F[k][eta] / N - 1. / (1 << m));
///			if (k == _k)
//				printf("%.10e\n", square(1. * F[k][0] / N - 1. / (1 << m)));
			S.push_back({s, k});
		}
		sort(S.begin(), S.end());
		/*
		print_bits(it.fst, 56);
		puts("");
		print_bits(key, 56);
		puts("");
		print_bits(key&it.fst, 56);
		puts("");
		print_bits(empty_with_mask(it.fst, key), 56);
		puts("");
		print_bits(fill_with_mask(it.fst, empty_with_mask(it.fst, key), 56), 56);
		puts("");
		puts("----------");
		*/
		int p = 0;
		for (auto it2: S) {
			if (it2.snd == empty_with_mask(key_mask, key, 56))
				printf("found pos %d/%d (%e .. %e .. %e)\n", int(S.size()) - p, int(S.size()), S[0].fst, it2.fst, S.back().fst);
			p++;
		}
	}
}

int main() {
	time_t seed = time(NULL);
	printf("seed = %u\n", (unsigned)seed);
//	srand(seed);
	assert(RAND_MAX == INT_MAX);
	key = rand() + ((1ull * rand() % (1 << 24)) << 32);
	extract();
	break_cipher();
	return 0;
}
