#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>

using namespace std;

const int N = 500000;
const double invN = 1.0/N;
const int c_int = 50;
const double c = (double)c_int;
const double invc = 1.0/c;
const int max_t_step = 5000;
const int epsilon_inv = 50;
const double epsilon = 1.0/epsilon_inv;
const int degmax = 1000;


int d_in, g_in;
int alpha1_in, alpha2_in;
int seed;

double alpha1, alpha2, alpha3;
double d, g;
int N1, N2, N_unnt, N3, N_pn;

int indeg[N];
int outdeg[N];
int j_in[N][degmax];
int cat_in[N][degmax];
int cat_out[N][degmax];

vector<double> m(N, 1.0);

void calc_rhs(const vector<double> &m, vector<double> &dm){
	double m_pene = 0.0;
	for(int i=N1; i<N; i++){
		m_pene += m[i];
	}
	double mu = alpha3 * d * g - d * m_pene * invN;
	for(int i=0; i<N; i++){
		double term_inc = 0.0, term_dec = 0.0;
		for(int j=0; j< indeg[i]; j++){
			term_inc += m[j_in[i][j]]*m[cat_in[i][j]];
		}
		for(int j=0; j< outdeg[i]; j++){
			term_dec += m[cat_out[i][j]];
		}
		term_dec *= m[i];
		double term_elastic;
		if(i < N1){
			term_elastic = 0.0;
		}else if(i < N_unnt){
			term_elastic = - d * m[i];
		}else{
			term_elastic = d * (g - m[i]);
		}
		double term_mu = mu * m[i];
		dm[i] = invc * (term_inc - term_dec) + term_elastic - term_mu;
	}
}

void rk4_step(vector<double> &m){
    vector<double> k1(N), k2(N), k3(N), k4(N), tmp(N);
    calc_rhs(m, k1);
    for (int i = 0; i < N; i++){
        tmp[i] = m[i] + 0.5 * epsilon * k1[i];
        if(tmp[i] < 0.0){
            tmp[i] = 1e-30;
        }
    }
    calc_rhs(tmp, k2);
    for (int i = 0; i < N; i++){
        tmp[i] = m[i] + 0.5 * epsilon * k2[i];
        if(tmp[i] < 0.0){
            tmp[i] = 1e-30;
        }
    }
    calc_rhs(tmp, k3);
    for (int i = 0; i < N; i++){
        tmp[i] = m[i] + epsilon * k3[i];
        if(tmp[i] < 0.0){
            tmp[i] = 1e-30;
        }
    }
    calc_rhs(tmp, k4);
    for (int i = 0; i < N; i++){
        m[i] += (epsilon / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        if(m[i] < 0.0){
            m[i] = 1e-30;
        }
    }
}

int main(int argc, char **argv){
	if(argc > 1){
		g_in = atoi(argv[1]);
		alpha1_in = atoi(argv[2]);
		alpha2_in = atoi(argv[3]);
		seed = atoi(argv[4]);
	}else{
		return 1;
	}
	g = g_in/100.0;
    alpha1 = alpha1_in/10000.0;
    alpha2 = alpha2_in/10000.0;
    alpha3 = 1.0 - alpha1-alpha2;

    N1 = (int)(N*alpha1+0.001);
	N2 = (int)(N*alpha2+0.001);
	N_unnt = N1 + N2;
	N3 = N - N_unnt;
	N_pn = N2 + N3;
	mt19937 gen(seed);
	uniform_real_distribution<double> uniform_real(0.0, 1.0);
	uniform_int_distribution<int> uniform_int_unnt(0, N_unnt-1);

	for(int i=0; i<N; i++){
		indeg[i] = 0;
		for(int j=0; j<degmax; j++){
			j_in[i][j] = 0;
			cat_in[i][j] = 0;
			cat_out[i][j] = 0;
		}
	}
	double norm_c = c/(double)N;
	int temp_k;
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			if(i != j){
				if(uniform_real(gen) < norm_c){
					j_in[i][indeg[i]] = j;
					temp_k = uniform_int_unnt(gen);
					if(temp_k == i || temp_k == j){
						int flg = 1;
						while(flg==1){
							temp_k = uniform_int_unnt(gen);
							if(temp_k != i && temp_k != j){
								flg = 0;
							}
						}
					}
					cat_in[i][indeg[i]] = temp_k;
					cat_out[j][outdeg[j]] = temp_k;
					indeg[i] += 1;
					outdeg[j] += 1;
				}
			}
		}
	}
	double mag, sq_mag;
	double mag1, mag2, mag3, mag_unnt, mag_pn;
	double sq_mag1, sq_mag2, sq_mag3, sq_mag_unnt, sq_mag_pn;
	double dilution_rate;
	for(int d_in = 20; d_in<21; d_in++){
		d = d_in*0.05;
		for(int i=0; i<N; i++){
			m[i] = 1.0;
		}
		fstream f1;
	    stringstream name1;
	    name1 << "result/simulation/" << N << "_" << c_int << "_" << d_in*5 << "_" << g_in << "_" << alpha1_in << "_" << alpha2_in << "_" << seed << "_sparse_long.txt";
		f1.open(name1.str().c_str(), ios::out);
		for(int t = 0; t<max_t_step; t++){
			if(t%5 ==0){
				mag = 0.0; mag1 = 0.0; mag2 = 0.0; mag3 = 0.0; mag_unnt = 0.0; mag_pn = 0.0;
				sq_mag = 0.0; sq_mag1 = 0.0; sq_mag2 = 0.0; sq_mag3 = 0.0; sq_mag_unnt = 0.0; sq_mag_pn = 0.0;
				for(int i=0; i<N; i++){
					double c_tempo = m[i];
					double cc_tempo = c_tempo*c_tempo;
					mag += c_tempo;
					sq_mag += cc_tempo;
					if(i<N_unnt){
						mag_unnt += c_tempo;
						sq_mag_unnt += cc_tempo;
					}else{
						mag3 += c_tempo;
						sq_mag3 += cc_tempo;
					}
					if(i<N1){
						mag1 += c_tempo;
						sq_mag1 += cc_tempo;
					}else{
						mag_pn += c_tempo;
						sq_mag_pn += cc_tempo;
						if(i<N_unnt){
							mag2 += c_tempo;
							sq_mag2 += cc_tempo;
						}
					}
				}
				dilution_rate = alpha3*d*g - d*mag_pn * invN;
				f1 << t*epsilon << " " //1 time
				   << mag/(double)N << " " //2 mean
				   << sq_mag/(double)N << " " //3 second_moment
				   << mag1/(double)N1 << " " //4 mean1
				   << mag2/(double)N2 << " " //5 mean2
				   << mag3/(double)N3 << " " //6 mean3
				   << mag_unnt/(double)N_unnt << " " //7 mean_unnt
				   << mag_pn/(double)N_pn << " " //8 mean_pn
				   << sq_mag1/(double)N1 << " " //9 sq_mean1
				   << sq_mag2/(double)N2 << " " //10 sq_mean2
				   << sq_mag3/(double)N3 << " " //11 sq_mean3
				   << sq_mag_unnt/(double)N_unnt << " " //12 sq_mean_unnt
				   << sq_mag_pn/(double)N_pn << " " //13 sq_mean_pn
				   << dilution_rate << endl; //14 dilution_rate
				if(t%500 ==0){
					fstream f3;
				    stringstream name3;
				    name3 << "result/simulation/" << N << "_" << c_int << "_" << d_in << "_" << g_in << "_" << alpha1_in << "_" << alpha2_in 
				    << "_" << (int)(t*epsilon) << "_" << seed << "_dist.txt";
				    f3.open(name3.str().c_str(), ios::out);
				    for(int i=0; i<N; i++){
				    	if(i < N1){
				    		f3 << i << " " << 1 << " " << indeg[i] << " " << outdeg[i] << " " << m[i] << endl; 
				    	}else if(i < N_unnt){
				    		f3 << i << " " << 2 << " " << indeg[i] << " " << outdeg[i] << " " << m[i] << endl; 
				    	}else{
				    		f3 << i << " " << 3 << " " << indeg[i] << " " << outdeg[i] << " " << m[i] << endl; 
				    	}   
				    }
				    f3.close();
				}
			}
			// heun_step(m);
			rk4_step(m);
		}
		mag = 0.0; mag1 = 0.0; mag2 = 0.0; mag3 = 0.0; mag_unnt = 0.0; mag_pn = 0.0;
		sq_mag = 0.0; sq_mag1 = 0.0; sq_mag2 = 0.0; sq_mag3 = 0.0; sq_mag_unnt = 0.0; sq_mag_pn = 0.0;
		for(int i=0; i<N; i++){
			double c_tempo = m[i];
			double cc_tempo = c_tempo*c_tempo;
			mag += c_tempo;
			sq_mag += cc_tempo;
			if(i<N_unnt){
				mag_unnt += c_tempo;
				sq_mag_unnt += cc_tempo;
			}else{
				mag3 += c_tempo;
				sq_mag3 += cc_tempo;
			}
			if(i<N1){
				mag1 += c_tempo;
				sq_mag1 += cc_tempo;
			}else{
				mag_pn += c_tempo;
				sq_mag_pn += cc_tempo;
				if(i<N_unnt){
					mag2 += c_tempo;
					sq_mag2 += cc_tempo;
				}
			}
		}
		dilution_rate = alpha3*d*g - d*mag_pn * invN;
		f1 << max_t_step*epsilon << " " //1 time
		   << mag/(double)N << " " //2 mean
		   << sq_mag/(double)N << " " //3 second_moment
		   << mag1/(double)N1 << " " //4 mean1
		   << mag2/(double)N2 << " " //5 mean2
		   << mag3/(double)N3 << " " //6 mean3
		   << mag_unnt/(double)N_unnt << " " //7 mean_unnt
		   << mag_pn/(double)N_pn << " " //8 mean_pn
		   << sq_mag1/(double)N1 << " " //9 sq_mean1
		   << sq_mag2/(double)N2 << " " //10 sq_mean2
		   << sq_mag3/(double)N3 << " " //11 sq_mean3
		   << sq_mag_unnt/(double)N_unnt << " " //12 sq_mean_unnt
		   << sq_mag_pn/(double)N_pn << " " //13 sq_mean_pn
		   << dilution_rate << endl; //14 dilution_rate
		f1.close();
	}
}