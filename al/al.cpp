#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <queue>
#include <vector>
using namespace std;
#define fst first
#define snd second
#define mp make_pair
typedef long long ll;
typedef pair<int, int> pii;
typedef pair<ll, ll> pll;

ll gcd(ll a, ll b) {
	if (b == 0) return a;
	return gcd(b, a % b);
}

struct Frac {
	ll num, den;

	Frac() { num = 0; den = 1; }
	Frac(ll num, ll den) : num(num), den(den) {}

	void printfloat() const {
		printf("%.8f", 1. * num / den);
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

	Frac operator * (const ll& x) const {
		Frac ans(num * x, den);
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

	bool operator == (const Frac& f) const {
		return num * f.den == den * f.num;
	}
	
	bool operator >= (const Frac& f) const {
		return f < *this || f == *this;
	}
};

struct Matrix {
	vector<vector<ll>> m;

	Matrix() {}
	Matrix(int n) {
		m = vector<vector<ll>>(n);
		for (int i = 0; i < n; ++i) {
			m[i] = vector<ll>(n);
			m[i][i] = 1;
		}
	}

	Matrix operator * (const Matrix& o) {
		int n = int(m.size());
		assert(n == int(o.m.size()));
		Matrix ans(n);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j) {
				ans.m[i][j] = 0;
				for (int k = 0; k < n; ++k)
					ans.m[i][j] += m[i][k] * o.m[k][j];
			}
		return ans;
	}

	ll norm1() const {
		ll ans = 0;
		int n = int(m.size());
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				ans += m[i][j];
		return ans;
	}

	ll trace() const {
		ll ans = 0;
		for (int i = 0; i < int(m.size()); ++i)
			ans += m[i][i];
		return ans;
	}

	void print() const {
		int n = int(m.size());
		for (int i = 0; i < n; ++i, puts(""))
			for (int j = 0; j < n; ++j)
				printf("%2lld ", m[i][j]);
	}

	bool operator < (const Matrix& o) const {
		return llabs(norm1()) > llabs(o.norm1());
	}

	bool operator == (const Matrix& o) const {
		int n = int(m.size());
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				if (m[i][j] != o.m[i][j])
					return false;
		return true;
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

ll applyp(ll x) {
	int ans = 0;
	for (int i = 31; i >= 0; --i)
		ans = (ans << 1) | ((x >> (P[i] - 1)) & 1);
	return ans;
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

int xor_all(int x) {
	int ans = 0;
	while (x > 0) {
		ans ^= (x % 2);
		x /= 2;
	}
	return ans;
}

int cut(int x, int n) {
	return x & ((1 << n) - 1);
}

void print_bits(ll x, int n) {
	for (int i = n - 1; i >= 0; --i) {
		printf("%lld", (x >> i) & 1);
		if ((n - i) % 4 == 0)
			printf(" ");
	}
}

pair<Matrix, int> aone[NB_SBOXES][1 << 8];

void compute_bias_sbox(int sbox, int mask) {
	Matrix& m = aone[sbox][mask].fst;
	aone[sbox][mask].snd = mask;
	m = Matrix(4);
	vector<int> com_pos_bef = get_com_pos(sbox, (sbox - 1 + NB_SBOXES) % NB_SBOXES);
	vector<int> com_pos_aft = get_com_pos(sbox, (sbox + 1) % NB_SBOXES);
	for (int com_val_bef = 0; com_val_bef < 4; ++com_val_bef)
		for (int com_val_aft = 0; com_val_aft < 4; ++com_val_aft) {
			ll& r = m.m[com_val_bef][com_val_aft];
			r = -2;
			for (int rem_val = 0; rem_val < 4; ++rem_val) {
				int in_val = merge_bits(com_pos_aft, rem_val, com_val_aft, 4);
				in_val = merge_bits(com_pos_bef, in_val, com_val_bef, 6);
				r += xor_all(in_val & cut(mask << 2, 6)) == xor_all(do_sbox(sbox, in_val) & (mask >> 4));
			}
//			r = abs(r);
		}
}

ll goal(vector<int> p) {
	Matrix m = aone[0][p[0]].fst;
	for (int i = 1; i < 8; ++i)
		m = m * aone[i][p[i]].fst;
	return m.trace();
}

void generate(vector<pair<Matrix, pll>>& v, Matrix m, ll in, ll out, int i, int imax, int c = 0) {
	if (i >= imax) {
		v.push_back({m, {in, out}});
		return;
	}
	for (int x = 0; x < (1 << 5); ++x) {
		int p = aone[i][x].snd;
		generate(v, m * aone[i][x].fst, (in << 4) | cut(p, 4), (out << 4) | (p >> 4), i + 1, imax, c + 1);
	}
}

map<ll, pair<ll, Frac>> mone;
vector<pair<pll, Frac>> vone;

void compute_best_approx_one() {
	for (int i = 0; i < 8; ++i)
		sort(aone[i], aone[i] + (1 << 8));
	/*
	for (int i = 0; i < 8; ++i)
		for (int j = 1; j < 256; ++j) {
			int k = rand() % j;
			swap(aone[i][j], aone[i][k]);
		}
	*/

	vector<pair<Matrix, pll>> left, right;
	generate(left, Matrix(4), 0, 0, 0, 4);
	generate(right, Matrix(4), 0, 0, 4, 8);
	puts("end of generation");
	sort(left.begin(), left.end());
	sort(right.begin(), right.end());
	priority_queue<pair<ll, pll>, vector<pair<ll, pll>>, greater<pair<ll, pll>>> q;
	for (int i = 0; i < (1 << 12); ++i) {
		auto l = left[i];
		for (int j = 0; j < (1 << 12); ++j) {
			auto r = right[j];
			ll g = (l.fst * r.fst).trace();
			q.push({abs(g), {(l.snd.fst << 16) | r.snd.fst, applyp((l.snd.snd << 16) | r.snd.snd)}});
			if (q.size() > 100000)
				q.pop();
		}
	}
	puts("end of combining");
	while (!q.empty()) {
		ll g = q.top().fst;
		/*
		print_bits(q.top().snd.fst, 32);
		puts("");
		print_bits(q.top().snd.snd, 32);
		puts("");
		printf("Bias = ");
		Frac(g, 1ll << 25).printfloat();
		puts("");
		*/
		mone[q.top().snd.snd] = {q.top().snd.fst, Frac(g, 1ll << 25)}; 
		vone.push_back({q.top().snd, Frac(g, 1ll << 25)});
		q.pop();
	}
	puts("end of push");
}

Frac B[25];
ll gammax[25], gammay[25];

void compute_best_approx(int levelmax, int level = 1, Frac f = Frac(1, 1)) {
	if (level > levelmax)
		return;
	if (level <= 2) {
		for (int i = int(vone.size()) - 1; int(vone.size()) - i < 1e3 && i >= 0; --i) {
			if (level == 1 && i == int(vone.size()) - 1)
				continue;

			gammay[level] = vone[i].fst.snd;
			gammax[level] = vone[i].fst.fst;
			Frac ftot = f * vone[i].snd;
			if (ftot * B[levelmax - level] * (1 << level) >= B[levelmax]) {
				if (level == levelmax) {
					assert(B[levelmax - level] == Frac(1, 1));
//					assert(ftot >= B[levelmax]);
					B[levelmax] = max(B[levelmax], ftot * (1 << (levelmax - 1)));
				}
				else
					compute_best_approx(levelmax, level + 1, ftot);
			}
		}
		return;
	}
	gammay[level] = gammay[level - 2] ^ gammax[level - 1];
	if (mone.find(gammay[level]) == mone.end())
		return;
	gammax[level] = mone[gammay[level]].fst;
	Frac ftot = f * mone[gammay[level]].snd;
	if (ftot * B[levelmax - level] * (1 << level) >= B[levelmax]) {
		if (level == levelmax) {
			assert(B[levelmax - level] == Frac(1, 1));
//			assert(ftot >= B[levelmax]);
			B[levelmax] = max(B[levelmax], ftot * (1 << (levelmax - 1)));
		}
		else
			compute_best_approx(levelmax, level + 1, ftot);
	}
}

int main() {
	srand(time(NULL));
	for (int i = 0; i < 8; ++i)
		for (int mask = 0; mask < (1 << 8); ++mask)
			compute_bias_sbox(i, mask);

	compute_best_approx_one();
	B[0] = Frac(1, 1);
	for (int level = 1; level <= 6; ++level) {
		compute_best_approx(level);
		printf("Best bias for level %d = ", level);
		B[level].printfloat();
		puts("");
	}
	return 0;
}
