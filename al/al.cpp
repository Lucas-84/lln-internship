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

ll gcd(ll a, ll b) {
	if (b == 0) return a;
	return gcd(b, a % b);
}

struct Frac {
	ll num, den;

	Frac() { num = 0; den = 1; }
	Frac(ll num, ll den) : num(num), den(den) {}

	void printfloat() const {
		printf("%.40f", 1. * num / den);
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
				ans += abs(m[i][j]); // here
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

ll applyp(ll x) {
	int ans = 0;
	for (int i = 31; i >= 0; --i) {
		ans = (ans << 1) | ((x >> (32 - P[31 - i])) & 1);
	}
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


pair<Matrix, int> aone[NB_SBOXES][1 << 8];

bool compare_aone(const pair<Matrix, int>& x, const pair<Matrix, int>& y) {
	bool is_x = (x.snd & (3 << 2)) == x.snd;
	bool is_y = (y.snd & (3 << 2)) == y.snd;
	if (is_x && !is_y) return true;
	if (is_y && !is_x) return false;
	return x.fst < y.fst;
}

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

const int NB_VECTORS = 32;
priority_queue<pair<double, ull>, vector<pair<double, ull>>, greater<pair<double, ull>>> bestscal[2][NB_VECTORS];
pthread_mutex_t bestscal_m;
double vect[NB_VECTORS][16];
struct timeval tstart, tnow;

void generate(int side, Matrix m, uint in, uint out, int i, int imax, bool print, int c = 0) {
	if (i >= imax) {
		for (int j = 0; j < NB_VECTORS; ++j) {
			double ps = 0;
			for (int k = 0; k < 4; ++k)
				for (int l = 0; l < 4; ++l)
					ps += m.m[k][l] * vect[j][4 * k + l];
//			m.print();
			pthread_mutex_lock(&bestscal_m);
			bestscal[side][j].push({fabs(ps), (ull(out) << 32) | in});
			if (int(bestscal[side][j].size()) > (1 << 12))
				bestscal[side][j].pop();
			pthread_mutex_unlock(&bestscal_m);
		}
		return;
	}
//	int upp = i == 0 || i == 4 ? 256 : (1 << 3);
//	int upp = i == 0 || i == 4 ? 256 : ((1 << 6) + (1 << 7));
//	int upp = i == 0 || i == 4 ? 256 : (1 << 6);
	int upp = 256;
	for (int x = 0; x < upp; ++x) {
		int p = aone[i][x].snd;
		if (x > 0 && print) {
			gettimeofday(&tnow, NULL);
			ull tnow_usec = tnow.tv_usec + 1000 * 1000ull * tnow.tv_sec;
			ull tstart_usec = tstart.tv_usec + 1000 * 1000ull * tstart.tv_sec;
			ull usec_per_it = (tnow_usec - tstart_usec) / x;
			ull it_rem = upp - x;
			ull usec_rem = usec_per_it * it_rem;
			ull sec_rem = usec_rem / 1000 / 1000;
			ull min_rem = sec_rem / 60;
			ull hour_rem = min_rem / 60;
			printf("%d/%d estimated time remaining: %lluh%llum%llus\n", x, upp, hour_rem, min_rem % 60, sec_rem % 60);
		}
		generate(side, m * aone[i][x].fst, (in << 4) | cut(p, 4), (out << 4) | (p >> 4), i + 1, imax, false, c + 1);
	}
}

void *generate_multi(void *arg) {
	int side = *(int *)arg;
	generate(side, Matrix(4), 0, 0, 4 * side, 4 * (side + 1), side == 0);
	pthread_exit(NULL);
	return NULL;
}

Matrix recompute_matrix(int side, uint in, uint out) {
	int first_box = side == 1 ? 7 : 3;
	int last_box = side == 1 ? 4 : 0;
	Matrix m(4);
	for (int i = first_box; i >= last_box; --i) {
		m = aone[i][cut(in, 4) | (cut(out, 4) << 4)].fst * m;	
		in >>= 4;
		out >>= 4;
	}
	return m;
}

map<uint, vector<pair<uint, Frac>>> mone;
vector<pair<ull, Frac>> vone;

ll shiftrot(ll x) {
	return cut(x << 1, 32) | ((x >> 31) & 1);
}

void compute_best_approx_one() {
	/*
	for (int i = 0; i < 8; ++i)
		sort(aone[i], aone[i] + (1 << 8), compare_aone);
	for (int i = 0; i < (1 << 8); ++i) {
		if (aone[2][i].snd == 0) {
			aone[0][i].fst.print();
			break;
		}
	}
	*/
	/*
	for (int i = 0; i < 8; ++i)
		for (int j = 1; j < 256; ++j) {
			int k = rand() % j;
			swap(aone[i][j], aone[i][k]);
		}
	*/
	for (int i = 0; i < NB_VECTORS; ++i) {
		double vnorm = 0;
		for (int j = 0; j < 16; ++j) {
			vect[i][j] = sqrt(-2 * log((1. * rand() + 1) / (RAND_MAX + 2.))) * cos(2 * M_PI * (1. * rand() + 1) / (RAND_MAX + 2.));
			vnorm += vect[i][j] * vect[i][j];
		}
		assert(vnorm > 0);
		vnorm = sqrt(vnorm);
		for (int j = 0; j < 16; ++j)
			vect[i][j] /= vnorm;
	}
	puts("vectors generated");
	pthread_mutex_init(&bestscal_m, NULL);
	pthread_t tid[2];
	int args[2];
	for (int i = 0; i < 2; ++i) {
		args[i] = i;
		pthread_create(&tid[i], NULL, generate_multi, &args[i]);
//		generate_multi(&args[i]);
	}
	for (int i = 0; i < 2; ++i)
		pthread_join(tid[i], NULL);
	pthread_mutex_destroy(&bestscal_m);
	puts("representants chosen");
	vector<pair<ull, Matrix>> sides[2];
	for (int i = 0; i < NB_VECTORS; ++i) {
		for (int side = 0; side < 2; ++side)
			while (!bestscal[side][i].empty()) {
				ull inout = bestscal[side][i].top().snd;
				sides[side].push_back({inout, recompute_matrix(side, cut(inout, 32), inout >> 32)});
				bestscal[side][i].pop();	
			}
	}
	puts("left and right sides generated");
	for (int i = 0; i < 2; ++i) {
		sort(sides[i].begin(), sides[i].end());
		auto it = unique(sides[i].begin(), sides[i].end());
		sides[i].resize(distance(sides[i].begin(), it));
		printf("Size %d = %d\n", i + 1, int(sides[i].size()));
	}
	priority_queue<pair<ll, ull>, vector<pair<ll, ull>>, greater<pair<ll, ull>>> q;
	for (auto l: sides[0]) {
		uint l_in = cut(l.fst, 32), l_out = l.fst >> 32;
		for (auto m: sides[1]) {
			uint m_in = cut(m.fst, 32), m_out = m.fst >> 32;
//			for (int k: ids[2]) {
//				auto n = sides[2][k];
//				for (int s: ids[3]) {
//					auto r = sides[3][s];
					ll g = (l.snd * m.snd).trace();
					q.push({llabs(g), shiftrot((l_in << 16) | m_in) | (applyp((l_out << 16)| m_out) << 32)});
					if (q.size() > (1 << 26))
						q.pop();
//				}
//			}
		}
//		printf("%f\n", 100. * i / sides[0].size());
	}

	puts("end of combining");
	int c = 0;
	while (!q.empty()) {
		ll g = q.top().fst;
		uint in = cut(q.top().snd, 32);
		uint out = q.top().snd >> 32;
		/*
		print_bits(in, 32);
		puts("");
		print_bits(out, 32);
		puts("");
		printf("Bias = ");
		Frac(g, 1ll << 25).printfloat();
		puts("");
		*/
//		assert(mone.find(out) == mone.end());
		if (in == ((1 << 15)) && out == ((1 << 7) | (1 << 18) | (1 << 24) | (1 << 29)))
			printf("A here: %d\n", c);
		if (in == ((1 << 27) | (1 << 28) | (1 << 30) | (1ll << 31)) && out == ((1 << 15)))
			printf("B here: %d\n", c);
		if (in == ((1 << 29)) && out == ((1 << 15)))
			printf("C here: %d\n", c);
		if (in == ((1 << 15)) && out == ((1 << 7) | (1 << 18) | (1 << 24)))
			printf("D here: %d\n", c);
		if (in == ((1 << 12) | (1 << 16)) && out == ((1 << 7) | (1 << 18) | (1 << 24)))
			printf("E here: %d\n", c);
		if (mone.find(out) == mone.end())
			mone[out] = {{in, Frac(g, 1ll << 25)}}; 
		else
			mone[out].push_back({in, Frac(g, 1ll << 25)});
		vone.push_back({q.top().snd, Frac(g, 1ll << 25)});
		q.pop();
		++c;
	}
	for (auto& it: mone)
		reverse(it.snd.begin(), it.snd.end());
	puts("end of push");
}

//Frac B[25];
priority_queue<double, vector<double>, greater<double>> B[25];
uint gammax[25], gammay[25];
const int NB_APPROX = 20000;

//void compute_best_approx(int levelmax, int level = 1, Frac f = Frac(1, 1), bool nonnul = false) {
void compute_best_approx(int levelmax, int level = 1, double f = 1., bool nonnul = false) {
	if (level > levelmax)
		return;
	if (level <= 2) {
		for (int i = int(vone.size()) - 1; int(vone.size()) - i < 1e5 && i >= 0; --i) {
//			if (level == 1 && i == int(vone.size()) - 1)
//				continue;

			gammay[level] = vone[i].fst >> 32;
			gammax[level] = cut(vone[i].fst, 32);
//			Frac ftot = f * vone[i].snd;
			double ftot = f * vone[i].snd.floatval();
			if (B[levelmax].size() < NB_APPROX || ftot * B[levelmax - level].top() * (1 << level) >= B[levelmax].top()) {
				if (level == levelmax) {
//					assert(B[levelmax - level] == Frac(1, 1));
//					assert(ftot >= B[levelmax]);
					if (nonnul || i != int(vone.size()) - 1) {
						B[levelmax].push(ftot * (1 << (levelmax - 1)));
						if (B[levelmax].size() > NB_APPROX)
							B[levelmax].pop();
					}
				}
				else
					compute_best_approx(levelmax, level + 1, ftot, nonnul || (i != int(vone.size() - 1)));
			}
		}
		return;
	}
	gammay[level] = gammay[level - 2] ^ gammax[level - 1];
	if (mone.find(gammay[level]) == mone.end())
		return;
	for (auto it: mone[gammay[level]]) {
		gammax[level] = it.fst;
//		Frac ftot = f * it.snd;
		double ftot = f * it.snd.floatval();
		if (B[levelmax].size() < NB_APPROX || ftot * B[levelmax - level].top() * (1 << level) >= B[levelmax].top()) {
			if (level == levelmax) {
//				assert(B[levelmax - level] == Frac(1, 1));
//				assert(ftot >= B[levelmax]);
				if (nonnul || gammay[level] != 0) {
					B[levelmax].push(ftot * (1 << (levelmax - 1)));
					if (B[levelmax].size() > NB_APPROX)
						B[levelmax].pop();
				}
			}
			else
				compute_best_approx(levelmax, level + 1, ftot, nonnul || gammay[level] != 0);
		}
//		return;
	}
}

void printdoublepower2(double x) {
	assert(0 <= x && x <= 1);
	if (x == 0) {
		printf("0");
		return;
	}
	ll p = 1;
	while ((1ll << (p + 1)) * x < 2)
		p++;
	printf("%.8f.2^-%lld", x * (1ll << p), p);
}

pair<double, double> lastpq(priority_queue<double, vector<double>, greater<double>> q) {
	double cap = 0;
	double l = q.top();
	while (!q.empty()) {
		l = q.top();
		cap += l * l;
		q.pop();
	}
	return {l, cap};
}

int main() {
	gettimeofday(&tstart, NULL);
	time_t seed = time(NULL);
	srand(seed);
	printf("seed = %ld\n", seed);
	for (int i = 0; i < 8; ++i)
		for (int mask = 0; mask < (1 << 8); ++mask)
			compute_bias_sbox(i, mask);
//	aone[0][(3 << 1) | (1 << 6)].fst.print();
//	aone[1][3 << 2].fst.print();
//	printf("%lld\n", (aone[0][(3 << 1) | (1 << 6)].fst * aone[1][3 << 2].fst).norm1());
//	printf("%f\n", 0.5 + 2 * 1. * (aone[0][(3 << 1) | (1 << 6)].fst * aone[1][3 << 2].fst).norm1() / (4 * 4) / (1 << 6));
//	return 0;
	compute_best_approx_one();
//	B[0] = Frac(1, 1);
	B[0].push(1);
//	for (int level = 1; level <= 16; ++level) {
//		for (int i = 0; i < NB_APPROX; ++i)
//			B[level].push(0);
//	}
	for (int level = 1; level <= 14; ++level) {
		compute_best_approx(level);
		printf("Best bias for level %d = ", level);
		printdoublepower2(B[level].top());
		printf(" ... ");
		auto p = lastpq(B[level]);
		printdoublepower2(p.fst);
		printf(" | Capacity = %.6e\n", 4 * p.snd);
//		B[level].printfloat();
	}
	return 0;
}
