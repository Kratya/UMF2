#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <ctime>
#include <vector>
#include <set>

using namespace std;

template <class mytype>
class SLAU
{
public:

	int n_x; //узлы областей

	int count_x; //узлы областей + внутренние

	int max_it = 1000;
	mytype eps = 1e-10;

	vector<mytype> ax; //узлы областей
	vector<mytype> x; //узлы c внутренними
	vector<int> nx;
	vector<mytype> qx;

	int m; //количество узлов

	mytype w = 1;//параметр релаксации

	//краевые
	int ncond1;
	vector<int> cond1;

	//матрицы и вектора
	vector<vector<mytype>> G;//локальные матрицы в плотном формате
	vector<vector<mytype>> M;
	vector<vector<mytype>> C;
	vector<mytype> b_loc;

	vector<mytype> Ggg;//глобальные матрицы в разреженном формате
	vector<mytype> Mgg;
	vector<mytype> Gdi;
	vector<mytype> Mdi;
	vector<mytype>Adi;
	vector<mytype>Agg;

	vector<mytype>Adi_new;
	vector<mytype> Agu_new;
	vector<mytype> Agl_new;

	vector<mytype> ig;//портрет
	vector<mytype> jg;//портрет
	vector<mytype> b; //глобальная пр часть
	vector<mytype> d; //глобальная пр часть

	vector<mytype> p_cur; //предыдущее решение по нелинейности
	vector<mytype> p_next; //следующее решение по нелинейности

	vector<mytype> q_cur; //текущее решение по времени
	vector<mytype> q_next; //следующее решение по времени

	vector<mytype> temp1, temp2, temp3;
	vector<mytype> r, z, p;

	int count_t, ntimes;
	vector<mytype> times, at, nt, qt;
	mytype k;

	SLAU<mytype>(string net_file, string cond_file_1, string time_file)
	{
		read_net(net_file);
		init_net();
		read_cond(cond_file_1);
		read_time(time_file);
		init_time();
		resize_func();
	};

	~SLAU<mytype>()
	{
		ax.~vector();
		x.~vector();
		nx.~vector();
		qx.~vector();
		cond1.~vector();

		for (int i = 0; i < 2; i++)
		{
			G[i].~vector();
			M[i].~vector();
			C[i].~vector();
		}

		G.~vector();
		M.~vector();
		C.~vector();
		b_loc.~vector();

		Ggg.~vector();
		Mgg.~vector();
		Gdi.~vector();
		Mdi.~vector();
		Adi.~vector();
		Agg.~vector();

		Adi_new.~vector();
		Agu_new.~vector();
		Agl_new.~vector();

		ig.~vector();
		jg.~vector();
		b.~vector();
		d.~vector();

		p_cur.~vector();
		p_next.~vector();
		q_cur.~vector();
		q_next.~vector();

		temp1.~vector();
		temp2.~vector();
		temp3.~vector();
		r.~vector();
		z.~vector();
		p.~vector();

		times.~vector();
		at.~vector();
		nt.~vector();
		qt.~vector();
	};

	void make_q(vector<mytype>& MV, mytype t);
	void read_net(string net_file);
	void init_net();
	void read_cond(string cond_file_1);
	void read_time(string time_file);
	void init_time();
	void resize_func();
	mytype get_lyam(mytype hx, mytype x, int p, vector<mytype>current);
	mytype get_sigma(mytype hx);
	mytype get_f(mytype x, mytype t);
	mytype get_Ug(mytype x, mytype t);
	void make_portrait();
	void make_local_G(int p, vector<mytype>current);
	void make_local_M(int p);
	void make_local_vec(int p, mytype t);
	void add_local_G(int p);
	void add_local_M(int p);
	void add_local_vec(int p);
	void make_global_G(vector<mytype> current);
	void make_global_M(mytype k);
	void make_global_vec(mytype t);
	void make_A();
	void make_d();
	void add_first(mytype t);
	void LU_sq();
	void mult(vector<mytype>& MV, vector<mytype> vec, vector<mytype>ggl, vector<mytype>ggu, vector<mytype>di);
	mytype skal_mult(vector<mytype> vec1, vector<mytype> vec2);
	mytype norm(vector<mytype> vec);
	void mult_pr(vector<mytype> aa, vector<mytype>& y, vector<mytype> bb);
	void mult_obr(vector<mytype> aa, vector<mytype>& y, vector<mytype> bb);
	void multU(vector<mytype>& MV, vector<mytype>vec);
	void mult_tr(vector<mytype>& MV, vector<mytype> vec);
	void LOS_sq(vector<mytype> sol_pred, vector<mytype>& sol_next);
	void do_smth();
};

template <typename mytype>
void SLAU<mytype>::make_q(vector<mytype>& MV, mytype t)
{
	for (int i = 0; i < m; i++)
	{
		MV[i] = get_Ug(x[i], t);
	}
}

template <typename mytype>
void SLAU<mytype>::read_net(string net_file)
{
	ifstream fin(net_file + ".txt");
	fin >> n_x;

	ax.resize(n_x);
	nx.resize(n_x - 1);
	qx.resize(n_x - 1);

	for (int i = 0; i < n_x; i++)
		fin >> ax[i];
	for (int i = 0; i < n_x - 1; i++)
		fin >> nx[i];
	for (int i = 0; i < n_x - 1; i++)
		fin >> qx[i];

	count_x = 1;
	for (int i = 0; i < n_x - 1; i++)
		count_x += nx[i];

	fin.close();
	m = count_x;
}

template <typename mytype>
void SLAU<mytype>::init_net()
{
	mytype hx, tx;
	for (int j = 0; j < n_x - 1; j++)
	{
		//обработка по х
		if (qx[j] != 1)
			hx = (ax[j + 1] - ax[j]) * (1. - qx[j]) / (1. - pow(qx[j], nx[j]));
		else hx = (ax[j + 1] - ax[j]) / nx[j];

		tx = ax[j];

		for (int k = 0; k < nx[j]; k++)
		{
			x.push_back(tx);
			tx += hx;
			hx *= qx[j];
		}
	}
	x.push_back(tx);
}

template <typename mytype>
void SLAU<mytype>::read_cond(string cond_file_1)
{
	ifstream fin1(cond_file_1 + ".txt");
	fin1 >> ncond1;
	cond1.resize(ncond1);
	for (int i = 0; i < ncond1; i++)
	{
		fin1 >> cond1[i];
	}
	fin1.close();
}

template <typename mytype>
void SLAU<mytype>::read_time(string time_file)
{
	ifstream fin1(time_file + ".txt");
	fin1 >> ntimes;

	at.resize(ntimes);
	nt.resize(ntimes - 1);
	qt.resize(ntimes - 1);

	for (int i = 0; i < ntimes; i++)
		fin1 >> at[i];

	for (int i = 0; i < ntimes - 1; i++)
		fin1 >> nt[i];
	for (int i = 0; i < ntimes - 1; i++)
		fin1 >> qt[i];

	count_t = 1;
	for (int i = 0; i < ntimes - 1; i++)
		count_t += nt[i];

	fin1.close();
}

template <typename mytype>
void SLAU<mytype>::init_time()
{
	mytype ht, tt;
	for (int j = 0; j < ntimes - 1; j++)
	{
		//обработка по t
		if (qt[j] != 1)
			ht = (at[j + 1] - at[j]) * (1. - qt[j]) / (1. - pow(qt[j], nt[j]));
		else ht = (at[j + 1] - at[j]) / nt[j];

		tt = at[j];

		for (int k = 0; k < nt[j]; k++)
		{
			times.push_back(tt);
			tt += ht;
			ht *= qt[j];
		}
	}
	times.push_back(tt);
}

template <typename mytype>
void SLAU<mytype>::resize_func()
{
	G.resize(2);
	M.resize(2);
	C.resize(2);

	for (int i = 0; i < 2; i++)
	{
		G[i].resize(2);
		M[i].resize(2);
		C[i].resize(2);
	}
	b_loc.resize(2);

	Ggg.resize(m - 1);
	Gdi.resize(m);
	Mgg.resize(m - 1);
	Mdi.resize(m);
	Agg.resize(m - 1);
	Adi.resize(m);

	Adi_new.resize(m);
	Agu_new.resize(m - 1);
	Agl_new.resize(m - 1);

	ig.resize(m + 1);
	jg.resize(m - 1);
	b.resize(m);
	d.resize(m);

	p_cur.resize(m);
	p_next.resize(m);
	q_cur.resize(m);
	q_next.resize(m);

	temp1.resize(m);
	temp2.resize(m);
	temp3.resize(m);
	r.resize(m);
	z.resize(m);
	p.resize(m);
}

template <typename mytype>
mytype SLAU<mytype>::get_lyam(mytype hx, mytype x, int p, vector<mytype>current)
{
	mytype arg = (current[p + 1] - current[p]) / hx;
	//return arg * arg * arg;
	return 1;
}

template <typename mytype>
mytype SLAU<mytype>::get_sigma(mytype hx)
{
	return 1;
}

template <typename mytype>
mytype SLAU<mytype>::get_f(mytype x, mytype t)
{
	//return -64 * x * x * x;
	//return -2;
	return 0;
}

template <typename mytype>
mytype SLAU<mytype>::get_Ug(mytype x, mytype t)
{
	//return x * x;
	return x;
}

template <typename mytype>
void SLAU<mytype>::make_portrait()
{
	ig[0] = 0;
	ig[1] = 0;
	int k = 1;
	for (int i = 2; i < m + 1; i++, k++)
	{
		ig[i] = k;
		jg[i - 2] = i - 2;
	}
}

template <typename mytype>
void SLAU<mytype>::make_local_G(int p, vector<mytype>current)
{
	mytype hx = x[p + 1] - x[p];
	mytype lyam1 = get_lyam(hx, x[p], p, current);
	mytype lyam2 = get_lyam(hx, x[p + 1], p, current);
	G[0][0] = (lyam1 + lyam2) / (2 * hx);
	G[0][1] = -G[0][0];
	G[1][0] = -G[0][0];
	G[1][1] = G[0][0];
}

template <typename mytype>
void SLAU<mytype>::make_local_M(int p)
{
	mytype hx = x[p + 1] - x[p];
	mytype sigma = get_sigma(hx);
	M[0][0] = sigma * hx * 2 / 6;
	M[0][1] = sigma * hx / 6;
	M[1][0] = M[0][1];
	M[1][1] = M[0][0];
}

template <typename mytype>
void SLAU<mytype>::make_local_vec(int p, mytype t)
{
	mytype f1, f2, hx;
	hx = x[p + 1] - x[p];
	f1 = get_f(x[p], t);
	f2 = get_f(x[p + 1], t);
	b_loc[0] = hx * (2 * f1 + f2) / 6;
	b_loc[1] = hx * (f1 + 2 * f2) / 6;
}

template <typename mytype>
void SLAU<mytype>::add_local_G(int p)
{
	Gdi[p] += G[0][0];
	Gdi[p + 1] += G[1][1];
	Ggg[p] += G[0][1];
}

template <typename mytype>
void SLAU<mytype>::add_local_M(int p)
{
	Mdi[p] += M[0][0];
	Mdi[p + 1] += M[1][1];
	Mgg[p] += M[0][1];
}

template <typename mytype>
void SLAU<mytype>::add_local_vec(int p)
{
	b[p] += b_loc[0];
	b[p + 1] += b_loc[1];
}

template <typename mytype>
void SLAU<mytype>::make_global_G(vector<mytype> current)
{
	for (int i = 0; i < m - 1; i++)
	{
		Gdi[i] = 0;
		Ggg[i] = 0;
	}
	Gdi[m - 1] = 0;

	for (int p = 0; p < count_x - 1; p++)
	{
		make_local_G(p, current);
		add_local_G(p);
	}
}

template <typename mytype>
void SLAU<mytype>::make_global_M(mytype k)
{
	for (int i = 0; i < m - 1; i++)
	{
		Mdi[i] = 0;
		Mgg[i] = 0;
	}
	Mdi[m - 1] = 0;
	for (int p = 0; p < count_x - 1; p++)
	{
		make_local_M(p);
		add_local_M(p);
	}

	for (int i = 0; i < m; i++)
		Mdi[i] *= k;
	for (int i = 0; i < Mgg.size(); i++)
		Mgg[i] *= k;
}

template <typename mytype>
void SLAU<mytype>::make_global_vec(mytype t)
{
	for (int i = 0; i < m; i++)
		b[i] = 0;
	for (int p = 0; p < count_x - 1; p++)
	{
		make_local_vec(p, t);
		add_local_vec(p);
	}
}

template <typename mytype>
void SLAU<mytype>::make_A()
{
	for (int i = 0; i < m; i++)
		Adi[i] = Gdi[i] + Mdi[i];
	for (int i = 0; i < Mgg.size(); i++)
		Agg[i] = Ggg[i] + Mgg[i];
}

template <typename mytype>
void SLAU<mytype>::make_d()
{
	mult(temp1, q_cur, Mgg, Mgg, Mdi);
	for (int i = 0; i < m; i++)
		d[i] = b[i] + k * temp1[i];
}

template <typename mytype>
void SLAU<mytype>::add_first(mytype t)
{
	long B = 1e+10;
	for (int i = 0; i < ncond1; i++)
	{
		int x_num;
		int ax_num = cond1[i];
		if (ax_num == 0) x_num = 0;
		else x_num = count_x - 1;

		Adi[x_num] = B;
		d[x_num] = B * get_Ug(x[x_num], t);
	}
}

//SLAE
template <typename mytype>
void SLAU<mytype>::LU_sq()
{
	//копирование-инициализация
	for (int i = 0; i < m; i++)
	{
		Adi_new[i] = Adi[i];
	}
	for (int i = 0; i < ig[m]; i++)
	{
		Agl_new[i] = Agg[i];
		Agu_new[i] = Agg[i];
	}

	for (int i = 0; i < m; i++)
	{
		mytype sd = 0; //переменные суммирования

		int i0 = ig[i];
		int i1 = ig[i + 1];

		for (int k = i0; k < i1; k++)
		{
			int j = jg[k];
			mytype sl = 0, su = 0;
			int j0 = ig[j];
			int j1 = ig[j + 1];
			int ki = i0;
			int kj = j0;

			for (; ki < k && kj < j1;)
			{
				int jl = jg[ki];
				int ju = jg[kj];
				if (jl == ju)
				{
					sl += Agu_new[kj] * Agl_new[ki];
					su += Agl_new[kj] * Agu_new[ki];
					ki++; kj++;
				}
				else if (jl < ju) ki++;
				else kj++;
			}
			Agu_new[k] = (Agu_new[k] - su) / Adi_new[j];
			Agl_new[k] = (Agl_new[k] - sl) / Adi_new[j];
			sd += Agu_new[k] * Agl_new[k];
		}

		Adi_new[i] = sqrt(Adi_new[i] - sd);
	}
}

template <typename mytype>
void SLAU<mytype>::mult(vector<mytype>& MV, vector<mytype> vec, vector<mytype>ggl, vector<mytype>ggu, vector<mytype>di)
{
	for (int i = 0; i < m; i++)
	{
		int k0 = ig[i];
		int k1 = ig[i + 1];
		MV[i] = di[i] * vec[i];
		for (int k = k0; k < k1; k++)
		{
			int j = jg[k];
			MV[i] += vec[j] * ggl[k];
			MV[j] += vec[i] * ggu[k];
		}
	}
}

template <typename mytype>
mytype SLAU<mytype>::skal_mult(vector<mytype> vec1, vector<mytype> vec2)
{
	mytype s = 0;
	for (int i = 0; i < m; i++)
	{
		s += vec1[i] * vec2[i];
	}
	return s;
}

template <typename mytype>
mytype SLAU<mytype>::norm(vector<mytype> vec)
{
	mytype sum = 0;
	for (int i = 0; i < m; i++)
		sum += vec[i] * vec[i];
	return sqrt(sum);
}

template <typename mytype>
void SLAU<mytype>::mult_pr(vector<mytype> aa, vector<mytype>& y, vector<mytype> bb)
{

	for (int i = 0; i < m; i++)
	{
		mytype s = 0; //переменные суммирования

		int i0 = ig[i];//индекс 1го элемента в iтой строке
		int i1 = ig[i + 1];

		for (int k = i0; k < i1; k++)
		{
			int j = jg[k];
			s += y[j] * aa[k];
		}
		y[i] = (bb[i] - s) / Adi_new[i];
	}
}

template <typename mytype>
void SLAU<mytype>::mult_obr(vector<mytype> aa, vector<mytype>& y, vector<mytype> bb)
{
	for (int i = 0; i < m; i++)
		y[i] = bb[i];
	for (int i = m - 1; i >= 0; i--)
	{
		int i0 = ig[i];//индекс 1го элемента в iтой строке
		int i1 = ig[i + 1];

		y[i] /= Adi_new[i];

		for (int k = i1 - 1; k >= i0; k--)
		{
			int j = jg[k];
			y[j] -= y[i] * aa[k];
		}
	}
}

template <typename mytype>
void SLAU<mytype>::multU(vector<mytype>& MV, vector<mytype>vec)
{
	for (int i = 0; i < m; i++)
	{
		int k0 = ig[i];
		int k1 = ig[i + 1];
		MV[i] = Adi[i] * vec[i];
		for (int k = k0; k < k1; k++)
		{
			int j = jg[k];
			//MV[i] += vec[j] * al[k];
			MV[j] += vec[i] * Agu_new[k];
		}
	}
}

template <typename mytype>
void SLAU<mytype>::mult_tr(vector<mytype>& MV, vector<mytype> vec)
{
	for (int i = 0; i < m; i++)
	{
		int k0 = ig[i];
		int k1 = ig[i + 1];
		MV[i] = Adi[i] * vec[i];
		for (int k = k0; k < k1; k++)
		{
			int j = jg[k];
			MV[i] += vec[j] * Agg[k];
			MV[j] += vec[i] * Agg[k];
		}
	}
}

template <typename mytype>
void SLAU<mytype>::LOS_sq(vector<mytype> sol_pred, vector<mytype>& sol_next)
{
	mytype skal1, skal2;
	int max_it = 10000;
	mytype err = 1e-20;
	LU_sq();

	for (int i = 0; i < sol_next.size(); i++)
		sol_next[i] = sol_pred[i];

	//      инициализация
	mult(temp1, sol_pred, Agg, Agg, Adi);
	for (int i = 0; i < m; i++)
	{
		temp2[i] = d[i] - temp1[i];
	}
	mult_pr(Agl_new, r, temp2);

	mult_obr(Agu_new, z, r);

	mult(temp1, z, Agg, Agg, Adi);
	mult_pr(Agl_new, p, temp1);

	//iteration
	mytype nev = skal_mult(r, r);
	for (int k = 0; k < max_it && nev > err; k++)
	{
		skal1 = skal_mult(p, r);
		skal2 = skal_mult(p, p);

		mytype alfa = skal1 / skal2;
		for (int i = 0; i < m; i++)
		{
			sol_next[i] += alfa * z[i];
			r[i] -= alfa * p[i];
		}
		mult_obr(Agu_new, temp1, r);
		mult(temp2, temp1, Agg, Agg, Adi);
		mult_pr(Agl_new, temp3, temp2);

		skal1 = skal_mult(p, temp3);

		mytype beta = -skal1 / skal2;

		mult_obr(Agu_new, temp2, r);

		for (int i = 0; i < m; i++)
		{
			z[i] = temp1[i] + beta * z[i];
		}

		for (int i = 0; i < m; i++)
		{
			p[i] = temp3[i] + beta * p[i];
		}

		nev = skal_mult(r, r);
	}
};

template <typename mytype>
void SLAU<mytype>::do_smth()
{
	cout << scientific << setprecision(8);

	ofstream fout("out.txt");
	fout << scientific << setprecision(8);
	//инициализация
	make_q(q_cur, times[0]);
	//for (int i = 0; i < m; i++)
		//q_cur[i] = 0;

	//сборка
	make_portrait();
	make_global_G(q_cur);
	make_global_M(1);

	//итерации по времени
	for (int cur_t = 1; cur_t < times.size(); cur_t++)
	{
		mytype dt = times[cur_t] - times[cur_t - 1];
		k = 1. / dt;
		make_global_vec(times[cur_t]);
		make_d();

		//итерации по нелинейности
		mytype nel_nev = 1;
		mytype dp = 1;

		int iter = 0;
		p_cur = q_cur;

		do
		{
			iter++;
			//сборка А d
			make_global_G(p_cur);
			make_global_M(k);
			make_A();

			add_first(times[cur_t]);
			LOS_sq(p_cur, p_next);

			//релаксация
			for (int i = 0; i < m; i++)
				p_next[i] = w * p_next[i] + (1 - w) * p_cur[i];



			make_global_G(p_next);
			make_global_M(k);
			make_A();
			make_global_vec(cur_t);
			make_d();
			add_first(times[cur_t]);
			mult(temp1, p_next, Agg, Agg, Adi);

			for (int i = 0; i < m; i++)
			{

				temp1[i] -= d[i];
				temp2[i] = p_next[i] - p_cur[i];
			}
			nel_nev = norm(temp1) / norm(d);
			dp = norm(temp2) / norm(p_next);
			p_cur = p_next;
			/*cout << iter << endl;*/
		} while (nel_nev > eps && dp > eps && iter < max_it);

		q_cur = p_next;
		cout << "iter=" << iter << endl;
		for (int i = 0; i < m; i++)
			fout << q_cur[i] << endl;
	}
}