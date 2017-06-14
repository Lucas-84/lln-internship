#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
using namespace std;

const int N = 100;
const int K = 500;

double p[N][4], q[N][4];
double b1[K], b2[K];

double randeps() {
	return 1. * rand() / RAND_MAX / 2;
}

void generate() {
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < 4; ++j) {
			p[i][j] = 0.5 + 1. * rand() / RAND_MAX;
			q[i][j] = 0.5 + 1. * rand() / RAND_MAX;
		}
}

void insert_sorted(double *a, int n, double v) {
	for (int i = 0; i < n; ++i)
		if (v > a[i]) {
			for (int j = n - 1; j > i; --j)
				a[j] = a[j - 1];
			a[i] = v;
			break;
		}
}

void search1() {
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j) {
			double v = 0;
			for (int k = 0; k < 4; ++k)
				v += p[i][k] * q[j][k] + (1 - p[i][k]) * (1 - q[j][k]);
			insert_sorted(b1, K, v);
		}
}

int g_s[4];
vector<double> buf[3 * 3 * 3 * 3][N];

double get_val(const vector<double>& a) {
	double v = 0;
	for (int i = 0; i < 4; ++i) {
		if (g_s[i] == 0) v += a[i];
		else if (g_s[i] == 2) v += 1 - a[i];
	}
	return v;
}

bool compare(const vector<double>& a, const vector<double>& b) {
	return get_val(a) < get_val(b);
}

void search2() {
	for (int i = 0; i < 3 * 3 * 3 * 3; ++i) {
		g_s[0] = i % 3;
		g_s[1] = (i / 3) % 3;
		g_s[2] = (i / 9) % 3;
		g_s[3] = i / 27;
		for (int j = 0; j < N; ++j) {
			buf[i][j] = vector<double>(4);
			for (int k = 0; k < 4; ++k)
				buf[i][j][k] = q[j][k];
		}
		sort(buf[i], buf[i] + N, compare);
	}
	for (int i = 0; i < N; ++i) {
		int power = 1;
		int k = 0;
		for (int j = 0; j < 4; ++j) {
			int s = 0;
			if (fabs(p[i][j] - 0.5) < 1e-6)
				s = 1;
			else if (p[i][j] > 0.5)
				s = 2;
			k += s * power;
			power *= 3;
		}
		for (int j = 0; j < 10; ++j) {
			double v = 0;
			for (int l = 0; l < 4; ++l)
				v += p[i][l] * buf[k][j][l] + (1 - p[i][l]) * (1 - buf[k][j][l]);
			insert_sorted(b2, K, v);
		}
	}
}

void print() {
	for (int i = 0; i < K; ++i)
		printf("[%d]\t%f\t%f\n", i, b1[i], b2[i]);
}

int main() {
//	srand(time(NULL));
	generate();	
	search1();
	search2();
	print();
	return 0;
}
