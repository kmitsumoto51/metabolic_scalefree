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
const double attra_ov_c = 1.0;
const double attra = attra_ov_c*c_int;
const int max_t_step = 5000;
const int epsilon_inv = 50;
const double epsilon = 1.0/epsilon_inv;

int d_in, g_in;
int alpha1_in, alpha2_in;
int seed;

double alpha1, alpha2, alpha3;
double d, g;
int N1=0, N2=0, N_unnt=0, N3=0, N_pn=0;

vector<vector<int>> j_out(N);
vector<vector<int>> j_in(N);
vector<vector<int>> cat_out(N);
vector<vector<int>> cat_in(N);
vector<int> indeg(N, 0);
vector<int> outdeg(N, 0);

vector<int> group_i(N);

vector<double> m(N, 1.0);

void calc_rhs(const vector<double> &m, vector<double> &dm){
    double m_pene = 0.0;
    double g_sum = 0.0;
    for(int i=0; i<N; i++){
        if(group_i[i] == 2){
            m_pene += m[i];
        }else if(group_i[i] == 3){
            m_pene += m[i];
            g_sum += g;
        }
    }
    double mu = d * (g_sum - m_pene) * invN;
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
        if(group_i[i] == 1){
            term_elastic = 0.0;
        }else if(group_i[i] == 2){
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
        d_in = atoi(argv[1]);
        g_in = atoi(argv[2]);
        alpha1_in = atoi(argv[3]);
        alpha2_in = atoi(argv[4]);
		seed = atoi(argv[5]);
	}else{
		return 1;
	}
	mt19937 gen(seed);
	uniform_real_distribution<double> uniform_real(0.0, 1.0);
  int current_nodes = c_int + 1;
  for(int i=0; i<c_int+1; i++){
      for(int j=0; j<c_int+1; j++){
          if(i == j) continue;
          j_out[i].push_back(j);   // i -> j
          j_in[j].push_back(i);    // incoming to j from i
      }
  }
  for(int i=0; i<N; ++i){
      indeg[i]  = (int)j_in[i].size();
      outdeg[i] = (int)j_out[i].size();
  }
  while(current_nodes < N){
      vector<int> pa_list;
      vector<double> pa_weight_list;
      vector<double> pa_weight_cumu_list;

      double pa_weight_sum = 0.0;
      for(int i=0; i<current_nodes; i++){
          pa_list.push_back(i);
          pa_weight_sum += (indeg[i] + attra);
      }

      double prob_sum = 0.0;
      for(int i=0; i<(int)pa_list.size(); i++){
          int v = pa_list[i];
          double p = (indeg[v] + attra) / pa_weight_sum;
          pa_weight_list.push_back(p);
          prob_sum += p;
          pa_weight_cumu_list.push_back(prob_sum);
      }

      for(int nl=0; nl<c_int; nl++){
          double r = uniform_real(gen);

          int pick_idx = -1;
          for(int i=0; i<(int)pa_list.size(); i++){
              if(r < pa_weight_cumu_list[i]){
                  pick_idx = i;
                  break;
              }
          }
          if(pick_idx < 0) pick_idx = (int)pa_list.size() - 1;

          int current_j = pa_list[pick_idx]; 

          j_out[current_nodes].push_back(current_j);
          j_in[current_j].push_back(current_nodes);

          indeg[current_j] += 1;

          double removed_p = pa_weight_list[pick_idx];

          pa_list.erase(pa_list.begin() + pick_idx);
          pa_weight_list.erase(pa_weight_list.begin() + pick_idx);
          pa_weight_cumu_list.erase(pa_weight_cumu_list.begin() + pick_idx);

          if(nl == c_int-1) break;

          double new_sum = 1.0 - removed_p;
          double inv_new_sum = 1.0 / new_sum;

          prob_sum = 0.0;
          for(int i=0; i<(int)pa_list.size(); i++){
              pa_weight_list[i] *= inv_new_sum;
              prob_sum += pa_weight_list[i];
              pa_weight_cumu_list[i] = prob_sum;
          }
      }
      outdeg[current_nodes] = (int)j_out[current_nodes].size();
      indeg[current_nodes]  = (int)j_in[current_nodes].size();

      current_nodes += 1;
  }

  for(int i=0; i<N; ++i){
      indeg[i]  = (int)j_in[i].size();
      outdeg[i] = (int)j_out[i].size();
  }
  d = d_in/100.0;
  g = g_in/100.0;
  alpha1 = alpha1_in/10000.0;
  alpha2 = alpha2_in/10000.0;
  alpha3 = 1.0 - alpha1-alpha2;

  for(int i=0; i<N; i++){
      double r = uniform_real(gen);
      if(r < alpha1){
          group_i[i] = 1;
          N1 += 1; N_unnt += 1;
      }else if(r < alpha1 + alpha2){
          group_i[i] = 2;
          N2 += 1; N_unnt += 1; N_pn += 1;
      }else{
          group_i[i] = 3;
          N3 += 1; N_pn += 1;
      }
  }

  uniform_int_distribution<int> uniform_int(0, N-1);

  cat_in.assign(N, {});
  for(int v=0; v<N; ++v){
      cat_in[v].assign(indeg[v], -1);
  }

  cat_out.assign(N, {});
  for(int u=0; u<N; ++u){
      cat_out[u].clear();
      cat_out[u].reserve(outdeg[u]);
  }

  for(int u=0; u<N; ++u){
      for(int e=0; e<outdeg[u]; ++e){
          int v = j_out[u][e]; 
          int k = uniform_int(gen); 

          if(k == u || k == v || group_i[k] == 3){
              while(true){
                  k = uniform_int(gen);
                  if(k != u && k != v && group_i[k] != 3) break;
              }
          }
          cat_out[u].push_back(k);
          int pos = -1;
          for(int t=0; t<indeg[v]; ++t){
              if(j_in[v][t] == u && cat_in[v][t] == -1){
                  pos = t;
                  break;
              }
          }
          if(pos < 0){
              std::cerr << "ERROR: cannot place catalyst for edge " << u << " -> " << v << "\n";
              std::exit(1);
          }
          cat_in[v][pos] = k;
      }
  }
  for(int v=0; v<N; ++v){
      for(int t=0; t<indeg[v]; ++t){
          if(cat_in[v][t] == -1){
              std::cerr << "ERROR: cat_in not filled at v=" << v << " t=" << t
                        << " source=" << j_in[v][t] << "\n";
              std::exit(1);
          }
      }
  }
  fstream f2;
  stringstream name2;
  name2 << "result/simulation/" << N << "_" << c_int << "_" << (int)(attra_ov_c*1000+0.01) << "_" 
  << d_in << "_" << g_in << "_" << alpha1_in << "_" << alpha2_in << "_" << seed << "_barabasi.txt";
  f2.open(name2.str().c_str(), ios::out);

  double mag, sq_mag;
  double mag1, mag2, mag3, mag_unnt, mag_pn;
  double sq_mag1, sq_mag2, sq_mag3, sq_mag_unnt, sq_mag_pn;
  double dilution_rate;
  for(int t = 0; t<max_t_step; t++){
      if(t%5 ==0){
          mag = 0.0; mag1 = 0.0; mag2 = 0.0; mag3 = 0.0; mag_unnt = 0.0; mag_pn = 0.0;
          sq_mag = 0.0; sq_mag1 = 0.0; sq_mag2 = 0.0; sq_mag3 = 0.0; sq_mag_unnt = 0.0; sq_mag_pn = 0.0;
          for(int i=0; i<N; i++){
              double c_tempo = m[i];
              double cc_tempo = c_tempo*c_tempo;
              mag += c_tempo;
              sq_mag += cc_tempo;
              if(group_i[i] == 1 || group_i[i] == 2){
                  mag_unnt += c_tempo;
                  sq_mag_unnt += cc_tempo;
              }else{
                  mag3 += c_tempo;
                  sq_mag3 += cc_tempo;
              }
              if(group_i[i] == 1){
                  mag1 += c_tempo;
                  sq_mag1 += cc_tempo;
              }else{
                  mag_pn += c_tempo;
                  sq_mag_pn += cc_tempo;
                  if(group_i[i] == 2){
                      mag2 += c_tempo;
                      sq_mag2 += cc_tempo;
                  }
              }
          }
          dilution_rate = alpha3*d*g - d*mag_pn * invN;
          f2 << t*epsilon << " " //1 time
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
      }
      rk4_step(m);
  }
  mag = 0.0; mag1 = 0.0; mag2 = 0.0; mag3 = 0.0; mag_unnt = 0.0; mag_pn = 0.0;
  sq_mag = 0.0; sq_mag1 = 0.0; sq_mag2 = 0.0; sq_mag3 = 0.0; sq_mag_unnt = 0.0; sq_mag_pn = 0.0;
  for(int i=0; i<N; i++){
      double c_tempo = m[i];
      double cc_tempo = c_tempo*c_tempo;
      mag += c_tempo;
      sq_mag += cc_tempo;
      if(group_i[i] == 1 || group_i[i] == 2){
          mag_unnt += c_tempo;
          sq_mag_unnt += cc_tempo;
      }else{
          mag3 += c_tempo;
          sq_mag3 += cc_tempo;
      }
      if(group_i[i] == 1){
          mag1 += c_tempo;
          sq_mag1 += cc_tempo;
      }else{
          mag_pn += c_tempo;
          sq_mag_pn += cc_tempo;
          if(group_i[i] == 2){
              mag2 += c_tempo;
              sq_mag2 += cc_tempo;
          }
      }
  }
  dilution_rate = alpha3*d*g - d*mag_pn * invN;
  f2 << max_t_step*epsilon << " " //1 time
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
  f2.close();

  fstream f3;
  stringstream name3;
  name3 << "result/simulation/" << N << "_" << c_int << "_" << (int)(attra_ov_c*1000+0.01) << "_" 
  << d_in << "_" << g_in << "_" << alpha1_in << "_" << alpha2_in << "_" << seed << "_BA_dist_steady.txt";
  f3.open(name3.str().c_str(), ios::out);
  for(int i=0; i<N; i++){
      f3 << i << " " << group_i[i] << " " << indeg[i] << " " << outdeg[i] << " " << m[i] << endl; 
  }
  f3.close();
}