#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <vector>
using namespace std;
#define fst first
#define snd second
typedef pair<int, int> pii;
typedef long long ll;

ll gcd(ll a, ll b) {
	if (b == 0) return a;
	return gcd(b, a % b);
}

struct Frac {
	ll num, den;

	Frac() {}
	Frac(ll num, ll den) : num(num), den(den) {}

	void printfloat() const {
		printf("%6.3f", 1. * num / den);
	}

	double floatval() const {
		return 1. * num / den;
	}

	void reduce() {
		ll g = gcd(llabs(num), llabs(den));
		num /= g;
		den /= g;
		if (den < 0) { num = -num; den = -den; }
	}

	Frac operator + (const Frac& f) const {
		Frac ans(num * f.den + f.num * den, den * f.den);
		ans.reduce();
		return ans;
	}

	Frac operator - (const Frac& f) const {
		Frac ans(num * f.den - f.num * den, den * f.den);
		ans.reduce();
		return ans;
	}

	Frac operator * (const Frac& f) const {
		Frac ans(num * f.num, den * f.den);
		ans.reduce();
		return ans;
	}

	Frac operator / (const ll& x) const {
		Frac ans(num, den * x);
		ans.reduce();
		return ans;
	}

	bool operator < (const Frac& f) const {
		return num * f.den < den * f.num;
	}
};

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

vector<pii> get_common_input(int sbox1, int sbox2) {
	assert(abs(min(sbox1 - sbox2, min(8 - sbox1 + sbox2, 8 - sbox2 + sbox1))) == 1);
	vector<pii> ans(2);
	ans[0].fst = -1;
	for (int i = 6 * sbox1; i < 6 * (sbox1 + 1); ++i)
		for (int j = 6 * sbox2; j < 6 * (sbox2 + 1); ++j) {
			if (E[i] == E[j]) {
				if (ans[0].fst == -1) { ans[0].fst = i - 6 * sbox1; ans[1].fst = j - 6 * sbox2; }
				else { ans[0].snd = i - 6 * sbox1; ans[1].snd = j - 6 * sbox2; }
			}
		}
	return ans;
}

int build_val_with(pii common, int rem_val, int common_val, int nb_bits) {
	int p = 1;
	int ans = 0;
	for (int i = 0; i < nb_bits; ++i) {
		if (common.fst == i) ans += p * (common_val % 2);
		else if (common.snd == i) ans += p * (common_val / 2);
		else { ans += p * (rem_val % 2); rem_val /= 2; }
		p *= 2;
	}
	return ans;
}

int do_sbox(int sbox, int in) {
	return S[sbox][(in & 1) | ((in >> 4) & 2)][(in >> 1) & ((1 << 4) - 1)];
}

int xor_all(int x) {
	int ans = 0;
	while (x > 0) {
		ans ^= (x % 2);
		x /= 2;
	}
	return ans;
}

// [sbox][in][key][out][val]
Frac papprox[NB_SBOXES][1 << 6][1 << 6][1 << 4][4];
vector<Frac> l1[1 << 16], l2[1 << 16];

void compute_best_approx_one(int sbox, pii common, int in_mask, int key_mask, int out_mask, int common_val) {
	for (int rem_in_val = 0; rem_in_val < (1 << 4); ++rem_in_val) {
		int in_val = build_val_with(common, rem_in_val, common_val, 6);
		for (int key_val = 0; key_val < (1 << 6); ++key_val) {
			assert(key_mask != in_mask || ((in_val & in_mask) ^ (key_val & key_mask)) == ((in_val ^ key_val) & in_mask));
			int in_xor = xor_all((in_val & in_mask) ^ (key_val & key_mask));
			int out_val = do_sbox(sbox, in_val ^ key_val);
			papprox[sbox][in_mask][key_mask][out_mask][common_val].den++;
			papprox[sbox][in_mask][key_mask][out_mask][common_val].num += in_xor == xor_all(out_val & out_mask);
		}
	}
}

Frac get_dist(vector<Frac> f1, vector<Frac> f2) {
	Frac ans(0, 1);
	for (int i = 0; i < 4; ++i)
		ans = ans + f1[i] * f2[i] + (Frac(1, 1) - f1[i]) * (Frac(1, 1) - f2[i]);
	return ans / 4;
}

int cut(int x, int n) {
	return x & ((1 << n) - 1);
}

void compute_best_approx(int sbox1, int sbox2) {
	vector<pii> common = get_common_input(sbox1, sbox2);	
	printf("Sboxes %d and %d have input %d and %d in common\n", sbox1, sbox2, common[0].fst, common[0].snd);
	double M = 0;
	for (int in_mask = 0; in_mask < (1 << 6); ++in_mask)
		for (int key_mask = 0; key_mask < (1 << 6); ++key_mask)
			for (int out_mask = 0; out_mask < (1 << 4); ++out_mask) {
				double x = 0;
				for (int common_val = 0; common_val < (1 << 2); ++common_val) {
					compute_best_approx_one(sbox1, common[0], in_mask, key_mask, out_mask, common_val);
					compute_best_approx_one(sbox2, common[1], in_mask, key_mask, out_mask, common_val);
					x += papprox[sbox1][in_mask][key_mask][out_mask][common_val].floatval();
				}
				l1[in_mask | (key_mask << 6) | (out_mask << 12)] = vector<Frac>(papprox[sbox1][in_mask][key_mask][out_mask], papprox[sbox1][in_mask][key_mask][out_mask] + 4);
				l2[in_mask | (key_mask << 6) | (out_mask << 12)] = vector<Frac>(papprox[sbox2][in_mask][key_mask][out_mask], papprox[sbox2][in_mask][key_mask][out_mask] + 4);
				if (in_mask != 0) {
					M = max(M, fabs(x / 4 - 0.5));
				}
			}
	printf("%.4f\n", M);
	double ans = 0;
	priority_queue<pair<double, pii>, vector<pair<double, pii>>, greater<pair<double, pii>>> q;
	const int qsz = 1 << 20;
//	for (int com_mask = 0; com_mask < (1 << 4); ++com_mask) {
//		if ((com_mask & (com_mask >> 2) & 1) || ((com_mask >> 1) & (com_mask >> 3) & 1))
//			continue;
	for (int in_mask = 0; in_mask < (1 << 12); ++in_mask) {
		for (int key_mask = 0; key_mask < (1 << 12); ++key_mask)
			for (int out_mask = 0; out_mask < (1 << 8); ++out_mask) {
				int i = cut(in_mask, 6) | (cut(key_mask, 6) << 6) | (cut(out_mask, 4) << 12);
				int j = (in_mask >> 6) | ((key_mask >> 6) << 6) | ((out_mask >> 4) << 12);
				if (l1[i].empty() || l2[j].empty())
					continue;
				double x = fabs(get_dist(l1[i], l2[j]).floatval() - 0.5);
				q.push({x, {i, j}});
				ans = max(ans, x);
				if (int(q.size()) > qsz)
					q.pop();
			}
		if (in_mask % 4 == 0)
			printf("%.4f%%\n", 100. * in_mask / (1 << 12));
	}
	printf("together = %.4f\n", ans);
	char name[64];
	sprintf(name, "approx%d%d", sbox1, sbox2);
	FILE *fp = fopen(name, "w+");
	while (!q.empty()) {
		fprintf(fp, "%d %d %f\n", q.top().snd.fst, q.top().snd.snd, q.top().fst);
		q.pop();
	}
	fclose(fp);
/*
	for (int in_mask = 1; in_mask < (1 << 6); ++in_mask) {
		int key_mask = in_mask;
		for (int out_mask = 1; out_mask < (1 << 4); ++out_mask) {
			Frac f(0, 1);
			for (int i = 0; i < 4; ++i)
				f = f + papprox[4][in_mask][key_mask][out_mask][i];
			f = (f / 4) - Frac(1, 2);
			f.printfloat();
			printf(" ");
		}
		puts("");
	}
*/
}

int main() {
	for (int i = 0; i < 8; ++i)
		compute_best_approx(i, (i + 1) % 8);	
//	compute_best_approx(0, 1);
/*
	for (int in_mask = 1; in_mask <= (1 << 5); ++in_mask) {
		for (int out_mask = 1; out_mask < (1 << 4); ++out_mask) {
			int c = 0;
			for (int in_val = 0; in_val < (1 << 6); ++in_val) {
				int out_val = do_sbox(4, in_val);
				if (xor_all(out_val & out_mask) == xor_all(in_val & in_mask))
					c++;
			}
			printf("%3d ", c - 32);
				
		}
		puts("");
	}
*/
	return 0;
}
