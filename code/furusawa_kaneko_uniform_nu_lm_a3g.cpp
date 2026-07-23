#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <algorithm>

using namespace std;

int D_in;
int alpha1_in;
int alpha2_in;
int nu_in;
int seed;
const double beta = 0.9;

static constexpr double XMIN = 1e-50;

double X_func(double a, double b){
    return log((a+(1+beta)*b)/(a+(1-beta)*b))/(2*beta*b);
}

double dX_func_a(double a, double b){
    if(beta > 0.99999 && beta < 1.00001){
        return -1.0/(a*(a+2*b));
    }else{
        return -1.0/((a+(1-beta)*b)*(a+(1+beta)*b));
    }
}

double dX_func_b(double a, double b){
    if(beta > 0.99999 && beta < 1.00001){
        return -log((a+2*b)/a)/(2*b*b) + 1.0/(b*(a+2*b));
    }
    return -log((a+(1+beta)*b)/(a+(1-beta)*b))/(2*beta*b*b) + a/(b*(a+(1-beta)*b)*(a+(1+beta)*b));
}

double calc_F1(double mu, double m_in_m_cat, double X2, double alpha1, double D, double D_alpha3_g){
    return D_alpha3_g - D*X2*((1-alpha1)*m_in_m_cat + D_alpha3_g) - mu;
}

double calc_F2(double m_in_m_cat, double X1, double X2, double alpha1, double D_alpha3_g){
    return  m_in_m_cat*(alpha1*X1 + (1-alpha1)*X2) + D_alpha3_g*X2 - 1;
}

double calc_F3(double m_in, double m_cat_nu_m1, double X1, double X2, double alpha1, double alpha2){
    return  m_in*m_cat_nu_m1*(alpha1*X1+alpha2*X2) - (alpha1+alpha2);
}

double calc_dF1_dmu(double m_in_m_cat, double dX2_dmu, double alpha1, double D, double D_alpha3_g){
    return - D*dX2_dmu*((1-alpha1)*m_in_m_cat + D_alpha3_g) - 1;
}

double calc_dF2_dmu(double m_in_m_cat, double dX1_dmu, double dX2_dmu, double alpha1, double D_alpha3_g){
    return m_in_m_cat*(alpha1*dX1_dmu + (1-alpha1)*dX2_dmu) + D_alpha3_g*dX2_dmu;
}

double calc_dF3_dmu(double m_in, double m_cat_nu_m1, double dX1_dmu, double dX2_dmu, double alpha1, double alpha2){
    return m_in*m_cat_nu_m1*(alpha1*dX1_dmu + alpha2*dX2_dmu);
}

double calc_dF1_dmcat(double m_in_m_cat, double nu_m_cat_nu_m1, double m_in, double X2, double dX2_dmcat,
                      double alpha1, double D, double D_alpha3_g){
    return - D*dX2_dmcat*nu_m_cat_nu_m1*((1-alpha1)*m_in_m_cat + D_alpha3_g)
           - D*(1-alpha1)*X2*nu_m_cat_nu_m1*m_in;
}

double calc_dF2_dmcat(double m_in_m_cat, double nu_m_cat_nu_m1, double m_in, double X1, double X2,
                      double dX1_dmcat, double dX2_dmcat, double alpha1, double D_alpha3_g){
    return m_in*nu_m_cat_nu_m1*(alpha1*X1 + (1-alpha1)*X2)
           + m_in_m_cat*nu_m_cat_nu_m1*(alpha1*dX1_dmcat + (1-alpha1)*dX2_dmcat)
           + D_alpha3_g*dX2_dmcat*nu_m_cat_nu_m1;
}

double calc_dF3_dmcat(double m_in, double nu_m1_m_in_m_cat_nu_m2, double nu_m_cat_nu_m1, double m_cat_nu_m1,
                      double X1, double X2, double dX1_dmcat, double dX2_dmcat, double alpha1, double alpha2){
    return nu_m1_m_in_m_cat_nu_m2*(alpha1*X1 + alpha2*X2)
           + m_in*m_cat_nu_m1*nu_m_cat_nu_m1*(alpha1*dX1_dmcat + alpha2*dX2_dmcat);
}

double calc_dF1_dmin(double m_cat_nu, double X2, double D, double alpha1){
    return -D*(1-alpha1)*X2*m_cat_nu;
}

double calc_dF2_dmin(double m_cat_nu, double X1, double X2, double alpha1){
    return m_cat_nu*(alpha1*X1+(1-alpha1)*X2);
}

double calc_dF3_dmin(double m_cat_nu_m1, double X1, double X2, double alpha1, double alpha2){
    return m_cat_nu_m1*(alpha1*X1+alpha2*X2);
}

double calc_energy(double mu, double m_in, double m_cat, double D, double nu_in,
                   double alpha1, double alpha2, double alpha3_g){
    double m_cat_nu = pow(m_cat, nu_in);
    double m_cat_nu_m1 = pow(m_cat, nu_in-1);
    double m_in_m_cat = m_in*m_cat_nu;
    double D_alpha3_g = D*alpha3_g;
    double X1 = X_func(mu, m_cat_nu);
    double X2 = X_func(mu+D, m_cat_nu);
    double F1 = calc_F1(mu, m_in_m_cat, X2, alpha1, D, D_alpha3_g);
    double F2 = calc_F2(m_in_m_cat, X1, X2, alpha1, D_alpha3_g);
    double F3 = calc_F3(m_in, m_cat_nu_m1, X1, X2, alpha1, alpha2);
    return F1*F1+F2*F2+F3*F3;
}


static inline void eval_F_J(double mu, double m_in, double m_cat,
                            double D, int nu_in, double alpha1, double alpha2, double alpha3_g,
                            double F[3], double J[3][3])
{
    double m_cat_nu     = pow(m_cat, nu_in);
    double m_cat_nu_m1  = pow(m_cat, nu_in-1);
    double m_in_m_cat   = m_in*m_cat_nu;

    double nu_m_cat_nu_m1 = nu_in * m_cat_nu_m1;

    double nu_m1_m_in_m_cat_nu_m2;
    if(nu_in == 1){
        nu_m1_m_in_m_cat_nu_m2 = 0.0;
    }else{
        nu_m1_m_in_m_cat_nu_m2 = (nu_in-1)*m_in*pow(m_cat, nu_in-2);
    }

    double D_alpha3_g = D*alpha3_g;

    double X1 = X_func(mu,   m_cat_nu);
    double X2 = X_func(mu+D, m_cat_nu);

    double dX1_dmu   = dX_func_a(mu,   m_cat_nu);
    double dX2_dmu   = dX_func_a(mu+D, m_cat_nu);
    double dX1_dmcat = dX_func_b(mu,   m_cat_nu);
    double dX2_dmcat = dX_func_b(mu+D, m_cat_nu);

    double F1 = calc_F1(mu, m_in_m_cat, X2, alpha1, D, D_alpha3_g);
    double F2 = calc_F2(m_in_m_cat, X1, X2, alpha1, D_alpha3_g);
    double F3 = calc_F3(m_in, m_cat_nu_m1, X1, X2, alpha1, alpha2);

    F[0]=F1; F[1]=F2; F[2]=F3;

    double dF1_dmu = calc_dF1_dmu(m_in_m_cat, dX2_dmu, alpha1, D, D_alpha3_g);
    double dF2_dmu = calc_dF2_dmu(m_in_m_cat, dX1_dmu, dX2_dmu, alpha1, D_alpha3_g);
    double dF3_dmu = calc_dF3_dmu(m_in, m_cat_nu_m1, dX1_dmu, dX2_dmu, alpha1, alpha2);

    double nu_m1_m_cat_nu_m2 = (nu_in==1 ? 0.0 : (nu_in-1)*pow(m_cat, nu_in-2));

    double dF1_dmcat = calc_dF1_dmcat(m_in_m_cat, nu_m_cat_nu_m1, m_in, X2, dX2_dmcat, alpha1, D, D_alpha3_g);
    double dF2_dmcat = calc_dF2_dmcat(m_in_m_cat, nu_m_cat_nu_m1, m_in, X1, X2, dX1_dmcat, dX2_dmcat, alpha1, D_alpha3_g);
    double dF3_dmcat = calc_dF3_dmcat(m_in, nu_m1_m_cat_nu_m2*m_in, nu_m_cat_nu_m1, m_cat_nu_m1,
                                      X1, X2, dX1_dmcat, dX2_dmcat, alpha1, alpha2);

    double dF1_dmin = calc_dF1_dmin(m_cat_nu, X2, D, alpha1);
    double dF2_dmin = calc_dF2_dmin(m_cat_nu, X1, X2, alpha1);
    double dF3_dmin = calc_dF3_dmin(m_cat_nu_m1, X1, X2, alpha1, alpha2);

    J[0][0]=dF1_dmu;   J[0][1]=dF1_dmcat;   J[0][2]=dF1_dmin;
    J[1][0]=dF2_dmu;   J[1][1]=dF2_dmcat;   J[1][2]=dF2_dmin;
    J[2][0]=dF3_dmu;   J[2][1]=dF3_dmcat;   J[2][2]=dF3_dmin;
}

static inline bool solve3(double A[3][3], double b[3], double x[3]){
    int p[3] = {0,1,2};

    for(int k=0;k<3;k++){
        int piv = k;
        double best = fabs(A[p[k]][k]);
        for(int i=k+1;i<3;i++){
            double v = fabs(A[p[i]][k]);
            if(v > best){ best=v; piv=i; }
        }
        if(best < 1e-300) return false;
        swap(p[k], p[piv]);

        for(int i=k+1;i<3;i++){
            double f = A[p[i]][k]/A[p[k]][k];
            for(int j=k;j<3;j++) A[p[i]][j] -= f*A[p[k]][j];
            b[p[i]] -= f*b[p[k]];
        }
    }

    for(int i=2;i>=0;i--){
        double s = b[p[i]];
        for(int j=i+1;j<3;j++) s -= A[p[i]][j]*x[j];
        x[i] = s / A[p[i]][i];
    }
    return true;
}

static inline bool solve_LM(double &mu, double &m_cat, double &m_in,
                            double D, int nu_in, double alpha1, double alpha2, double alpha3_g,
                            int max_iter=200, double tolE=1e-20, double tolStep=1e-10)
{
    double lambda = 1e-3; 
    const double lambda_up = 10.0;
    const double lambda_dn = 0.3;

    double E = calc_energy(mu, m_in, m_cat, D, nu_in, alpha1, alpha2, alpha3_g);

    for(int it=0; it<max_iter; it++){
        double F[3], J[3][3];
        eval_F_J(mu, m_in, m_cat, D, nu_in, alpha1, alpha2, alpha3_g, F, J);

        double A[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        double g[3]    = {0,0,0};

        for(int i=0;i<3;i++){
            for(int a=0;a<3;a++){
                g[a] += J[i][a]*F[i];
                for(int b=0;b<3;b++){
                    A[a][b] += J[i][a]*J[i][b];
                }
            }
        }

        double Ad[3][3];
        for(int a=0;a<3;a++){
            for(int b=0;b<3;b++) Ad[a][b]=A[a][b];
            Ad[a][a] += lambda;
        }

        double bvec[3] = {-g[0], -g[1], -g[2]};
        double dx[3] = {0,0,0};
        if(!solve3(Ad, bvec, dx)){
            return false;
        }

        double mu_new   = mu    + dx[0];
        double mcat_new = m_cat + dx[1];
        double min_new  = m_in  + dx[2];

        if(mu_new   < XMIN) mu_new   = XMIN;
        if(mcat_new < XMIN) mcat_new = XMIN;
        if(min_new  < XMIN) min_new  = XMIN;

        double E_new = calc_energy(mu_new, min_new, mcat_new, D, nu_in, alpha1, alpha2, alpha3_g);

        if(E_new < E){
            double step_norm = fabs(dx[0]) + fabs(dx[1]) + fabs(dx[2]);

            mu = mu_new; m_cat = mcat_new; m_in = min_new;
            E = E_new;

            lambda = max(lambda*lambda_dn, 1e-16);

            if(E < tolE) return true;
            if(step_norm < tolStep) return true;
        }else{
            lambda = min(lambda*lambda_up, 1e30);
            if(lambda > 1e25) return false;
        }
    }
    return (E < tolE);
}

int main(int argc, char **argv){
    if(argc > 1){
        D_in = atoi(argv[1]);
        alpha1_in = atoi(argv[2]);
        alpha2_in = atoi(argv[3]);
        nu_in = atoi(argv[4]);
    }
    double D_div;
    if(nu_in == 1){
        D_div = 0.2;
    }else if(nu_in == 2){
        D_div = 0.02;
    }else{
        D_div = 0.006;
    }
    double D = D_in*D_div;
    double alpha1 = alpha1_in/10000.0;
    double alpha2 = alpha2_in/10000.0;
    double alpha_unnt = alpha1 + alpha2;
    double alpha3 = 1 - alpha_unnt;
    double mu, m_cat, m_in;

    fstream f2;
    stringstream name2;
    name2 << "result/270_" << alpha1_in << "_" << alpha2_in << "_" << (int)(beta*1000) << "_" << nu_in
          << "_mean_uniform_Ddep_LM.txt";
    f2.open(name2.str().c_str(), ios::in);
    if(!f2){
        cerr << "Cannot open init file: " << name2.str() << endl;
        return 1;
    }
    for(int D_out = 1; D_out<101; D_out++){
        double a1,a2,a3_,a4,a5,a6,a7,a8,a9,a10,a11;
        int a12;
        f2 >> a1 >> a2 >> a3_ >> a4 >> a5 >> a6 >> a7 >> a8 >> a9 >> a10 >> a11 >> a12;

        if(a1 > D-0.00001 && a1 < D+0.00001){
            mu = a9;
            m_cat = a10;
            m_in = a11;
        }
    }
    f2.close();

    fstream f1;
    stringstream name1;
    name1 << "result/" << D_in << "_" << alpha1_in << "_" << alpha2_in << "_" << (int)(beta*1000) << "_" << nu_in
          << "_mean_uniform_a3gdep_LM.txt";
    f1.open(name1.str().c_str(), ios::out);

    for(int alpha3_g_in = 5480; alpha3_g_in>0; alpha3_g_in--){
        double alpha3_g = alpha3_g_in/1000.0;
        double g_ratio = alpha3_g/alpha3;

        bool ok = solve_LM(mu, m_cat, m_in, D, nu_in, alpha1, alpha2, alpha3_g,
                           /*max_iter=*/10000, /*tolE=*/1e-20, /*tolStep=*/1e-15);

        double E = calc_energy(mu, m_in, m_cat, D, nu_in, alpha1, alpha2, alpha3_g);

        double m_cat_nu = pow(m_cat, nu_in);
        double m_cat_nu_m1 = pow(m_cat, nu_in-1);
        double m_in_m_cat = m_in*m_cat_nu;
        double D_alpha3_g = D*alpha3_g;
        double X1 = X_func(mu, m_cat_nu);
        double X2 = X_func(mu+D, m_cat_nu);
        double F1v = calc_F1(mu, m_in_m_cat, X2, alpha1, D, D_alpha3_g);
        double F2v = calc_F2(m_in_m_cat, X1, X2, alpha1, D_alpha3_g);
        double F3v = calc_F3(m_in, m_cat_nu_m1, X1, X2, alpha1, alpha2);
        double m1 = X1*m_in_m_cat;
        double m2 = X2*m_in_m_cat;
        double m3 = X2*(D*g_ratio + m_in_m_cat);

        f1 << alpha3_g          //1
           << " " << E          //2
           << " " << F1v        //3
           << " " << F2v        //4
           << " " << F3v        //5
           << " " << m1         //6
           << " " << m2         //7
           << " " << m3         //8
           << " " << mu         //9
           << " " << m_cat      //10
           << " " << m_in       //11
           << " " << (ok?1:0)   //12 
           << endl;

    }

    f1.close();
    return 0;
}