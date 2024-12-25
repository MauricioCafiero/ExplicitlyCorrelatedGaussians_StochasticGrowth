#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include <bits/stdc++.h>
#include <algorithm>

#define M_PI 3.14159265358979323846

using namespace std;
using Eigen::MatrixXf;
using Eigen::VectorXf;
//NECG: He = -2.90338583; Li+ = -7.2796674227; H- =-0.5277157296; Li = -7.47806033; Be = −14.698

string get_Header(int, int, int);

default_random_engine& global_urng(){
    static default_random_engine u{};
    return u;
}

class Energy
{
    public:
        vector<float> vec1;
        int Z, n, m;
        //Nuclear charge, number of electrons, number of basis functions

        void get_vec1 () {
            int i;
            for (i=0; i < vec1.size(); i++){
                cout << "component " << i << " : " << vec1.at(i) << endl;
            }
        }

        Energy(int n_in, int Z_in)
        {
            Z = Z_in;
            n = n_in;
        }

        MatrixXf getMat();

        float calc_E (int m, vector<float> vec1, float* ev) {

            int i,j,k,l,si,bi,ii,jj,low_ind,size_bf = n*(n+1)/2;
            float gs_energy;
            vector<float> energies(m);
            float detAk, detAl, detAkl, T, ER, NA, OV;
            float Z_f = static_cast<float>(Z);
            MatrixXf S(m,m), L(m,m), H(m,m);
            VectorXf Sc(m);
            MatrixXf Ak(n,n),Al(n,n),Jij(n,n),Jii(n,n),Akl(n,n),Akli(n,n),im1(n,n),im2(n,n),ll(n,n),lk(n,n);
            S = MatrixXf::Zero(m,m);
            H = MatrixXf::Zero(m,m);

            //for symmetry: S_total is all matrices for n=4; S_temp hold each matrix for all n; facs = factorial hash table
            //symCoeffs: coefficients for each symmetry term.
            MatrixXf S_temp(n,n);
            int facs[] = {1,1,2,6,24,120};
            float symCoeffs[3][24] = {{1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{4,-2,4,-2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
                                    {16,8,-16,8,-8,-8,8,-16,8,16,16,16,-8,-8,-8,-8,-8,-8,-16,8,8,8,8,-16}};
            MatrixXf S_total(4*facs[4],4);
            S_total = getMat();

            for (si=0; si<facs[n]; si++){

                S_temp = S_total.block(si*4,0,n,n);

                for (i=0; i < m; i++){
                    for (j=0; j < i+1; j++) {
                        //read in lower diagonal and create Ak, Al (makes it positive definite)
                        lk = MatrixXf::Zero(n,n);
                        ll = MatrixXf::Zero(n,n);
                        k = i*size_bf;
                        l = j*size_bf;
                        bi = 0;
                        for (jj=0; jj < n; jj++){
                            for (ii=jj; ii < n; ii++){
                                lk(ii,jj) = vec1.at(k+bi);
                                ll(ii,jj) = vec1.at(l+bi);
                                bi += 1;
                            }
                        }
                        Ak = lk*lk.transpose();
                        Al = ll*ll.transpose();

                        //apply symmetry to the ket
                        im1 = Al*S_temp;
                        im2 = S_temp.transpose()*im1;
                        Al = im2;

                        Akl = Ak + Al;
                        Akli = Akl.inverse();

                        detAk = Ak.determinant();
                        detAk = sqrt(detAk);
                        detAl = Al.determinant();
                        detAl = sqrt(detAl);
                        detAkl = Akl.determinant();

                        OV = sqrt(pow(8.0,n))*pow((detAk*detAl/detAkl),1.5);
                        S(i,j) += symCoeffs[n-2][si]*OV;

                        im1 = Akli*Al;
                        im2 = Ak*im1;
                        T = OV*3.0*im2.trace();

                        ER = 0.0;
                        for (ii=0; ii < n; ii++){
                            for (jj=ii+1; jj < n; jj++){
                                Jij = MatrixXf::Zero(n,n);
                                Jij(ii,ii) = Jij(jj,jj) =  1.0;
                                Jij(ii,jj) = Jij(jj,ii) = -1.0;
                                im1 = Jij*Akli;
                                ER += OV*2.0/sqrt(M_PI*im1.trace());
                            }
                        }

                        NA = 0.0;
                        for (ii=0; ii < n; ii++){
                            Jii = MatrixXf::Zero(n,n);
                            Jii(ii,ii) = 1.0;
                            im1 = Jii*Akli;
                            NA += Z_f*OV*2.0/sqrt(M_PI*im1.trace());
                        }
                        H(i,j) += symCoeffs[n-2][si]*(T - NA + ER);
                    }
                }
            }

            //linear dependence check; saves to separate file
            ofstream linear_dep("Linear_Dependence_info.txt", ios::app);
            if (linear_dep){
                Eigen::SelfAdjointEigenSolver<MatrixXf> ld(S);
                for (i=0; i < m; i++){
                    if (ld.eigenvalues()[i] < 0.001){
                        for (j=0; j < m; j++ ){
                            if (ld.eigenvectors()(j,i) > 0.95){
                                linear_dep << "Detected a possible linear dependence with function: " << i << endl;
                                linear_dep << "Large eigenvector coefficient: " << j << " = " << ld.eigenvectors()(j,i) << endl;
                            }
                        }
                    }
                }
            }
            linear_dep.close();

            Eigen:: GeneralizedSelfAdjointEigenSolver<MatrixXf> es(H,S);

            //get lowest eigenvalue and corresponding eigenvector
            low_ind = 0;
            for (i=1; i < es.eigenvalues().size(); i++){
                if (es.eigenvalues()[i] < es.eigenvalues()[low_ind]){
                    low_ind = i;
                }
            }
            gs_energy = es.eigenvalues()[low_ind];

            Eigen::LLT<MatrixXf> lltOfS(S);
            L = lltOfS.matrixL();
            Sc = L.transpose()*es.eigenvectors().col(low_ind);

            for (i=0; i < m; i++){
                    ev[i] = Sc(i);
                    energies[i] = es.eigenvalues()[i];
            }

            //write all eigenvalues to a separate file
            sort(begin(energies),end(energies));
            ofstream eig_writer("all_eigenvalues_new.csv");
            if (eig_writer){
                for (i=0; i < m; i++){
                    eig_writer << energies[i] << endl;
                }
            }
            eig_writer.close();

            return gs_energy;
        }

	int overlap (int m, vector<float> vec1) {

            int i,j,k,l,si,bi,ii,jj,size_bf = n*(n+1)/2,isOverlap;
            vector<float> overlaps(m);
            float detAk, detAl, detAkl,OV, s_cutoff;
            MatrixXf S(m,m);
            MatrixXf Ak(n,n),Al(n,n),Akl(n,n),Akli(n,n),ll(n,n),lk(n,n);
            S = MatrixXf::Zero(m,m);

            for (i=0; i < m; i++){
                for (j=0; j < i+1; j++) {
                    //read in lower diagonal and create Ak, Al (makes it positive definite)
                    lk = MatrixXf::Zero(n,n);
                    ll = MatrixXf::Zero(n,n);
                    k = i*size_bf;
                    l = j*size_bf;
                    bi = 0;
                    for (jj=0; jj < n; jj++){
                        for (ii=jj; ii < n; ii++){
                            lk(ii,jj) = vec1.at(k+bi);
                            ll(ii,jj) = vec1.at(l+bi);
                            bi += 1;
                        }
                    }
                    Ak = lk*lk.transpose();
                    Al = ll*ll.transpose();

                    Akl = Ak + Al;
                    Akli = Akl.inverse();

                    detAk = Ak.determinant();
                    detAk = sqrt(detAk);
                    detAl = Al.determinant();
                    detAl = sqrt(detAl);
                    detAkl = Akl.determinant();

                    OV = sqrt(pow(8.0,n))*pow((detAk*detAl/detAkl),1.5);
                    S(i,j) = OV;
                }
            }

            //linear dependence check
            s_cutoff = 0.9;
            isOverlap = -1;
            for (i = 0; i < m; i++){
                for (j = 0; j < i+1; j++){
                    OV = S(i,j)/(S(i,i)*S(j,j));
                    if (OV > s_cutoff && i != j){
                        cout << "Large overlap between functions " << i << " and " << j << " with value: " << OV << endl;
                        isOverlap = i;
                    }
                }
            }

            return isOverlap;
        }

};

int main()
{
    vector<float> coeffs(0), good_coeffs(0);
    float gs_k, gs_kp1, rn, sum_c=0.0, diff;
    string st_temp, header;
    stringstream ss_temp;
    int i,j,k,m,n,size_bf,Z,vec_count,basis_i,max_basis, opt_iters, max_iters, isOverlap;

    cout << "Calculate 2-4 electron atomic energies with explicitly correlated Gaussian functions: " << endl;

    ifstream reader("coeffs_in.csv");
    if(!reader) {
        cout << "Could not open file!" << endl;
    }
    else {
        cout << "File opened!" << endl;
        
        ss_temp.clear();
        getline(reader,st_temp,'\n');
	
	ss_temp.clear();
        getline(reader,st_temp,'\n');
        ss_temp << st_temp;
        ss_temp >> Z;

        ss_temp.clear();
        getline(reader,st_temp,'\n');
	
        ss_temp.clear();
        getline(reader,st_temp,'\n');
        ss_temp << st_temp;
        ss_temp >> n;

        ss_temp.clear();
        getline(reader,st_temp,'\n');
	
        ss_temp.clear();
        getline(reader,st_temp,'\n');
        ss_temp << st_temp;
        ss_temp >> max_basis;
        
        ss_temp.clear();
        getline(reader,st_temp,'\n');
	
        ss_temp.clear();
        getline(reader,st_temp,'\n');
        ss_temp << st_temp;
        ss_temp >> max_iters;
	
	size_bf = n*(n+1)/2;
    }

    //cout << "Read in " << m << " basis functions." << endl;
    cout << "Read in atomic number: " << Z << endl;
    cout << "Read in number of electrons: " << n << endl;
    cout << "Read in maximum number of basis functions: " << max_basis << endl;
    cout << "Read in number of iterations per basis function: " << max_iters << endl;
    cout << "==========================================" << endl;

    ofstream linear_dep("Linear_Dependence_info.txt");
    linear_dep.close();
    linear_dep.clear();

    //get filename with atom type, electrons and basis functions
    header = get_Header(Z,n,max_basis);
    header += "_coeffs_new.csv";

    //Grow a Basis Set ======================================================================
    cout << "========== Individual function growth and stochastic optimzation ===========" << endl;
    default_random_engine re(random_device{}());
    normal_distribution<float> gauss_bs_d(1.0,2.5);
    normal_distribution<float> gauss_bs_od(0.01,0.05);
    normal_distribution<float> gauss_d(0.0,0.03);
    normal_distribution<float> gauss_od(0.0,0.015);

    Energy he_atom(n,Z);

    m = 0;
    while (m < max_basis){
        for (j = 0; j < n; j++) {
            for (k = j; k < n; k++){
                if (j == k){
                    rn = gauss_bs_d(re);}
                else {
                    rn = gauss_bs_od(re);}
                good_coeffs.push_back(rn);
            }
        }
        m += 1;
        cout << "Current number of Basis Functions: " << good_coeffs.size()/size_bf << endl;

	isOverlap = he_atom.overlap(m,good_coeffs);
        if (isOverlap !=-1){
            j = isOverlap*size_bf-1;
            good_coeffs.erase(good_coeffs.begin()+j,good_coeffs.begin()+j+size_bf);
            cout << "Removing offending basis function: " << isOverlap << endl;
            m -= 1;
            continue;
        }

        float* eigenvector = new float [m]; 

        gs_k = he_atom.calc_E(m,good_coeffs,eigenvector); 
        cout << "Ground state energy at iteration 1 is: " << gs_k << endl;

        coeffs.assign(good_coeffs.begin(), good_coeffs.end());
        for (opt_iters = 0; opt_iters < max_iters; opt_iters++){
            vec_count = (m-1)*size_bf;
            for (j = 0; j < n; j++) {
                for (k = j; k < n; k++){
                    if (j == k){
                        rn = gauss_d(re);}
                    else {
                        rn = gauss_od(re);}
                    coeffs[vec_count] += rn;
                vec_count += 1;
                }
            }

            gs_kp1 = he_atom.calc_E(m,coeffs,eigenvector);

            if (gs_kp1 < gs_k){
                cout << "Ground state energy at iteration "<< opt_iters <<  " improved: " << gs_kp1 << endl;
                diff = abs(gs_kp1 - gs_k);
                if (diff > 10.0){
                    cout << "Numerical instability: resetting energy!" << endl;
                }
                else{
                    gs_k = gs_kp1;
                    good_coeffs.assign(coeffs.begin(), coeffs.end());
                }
            }
        }
    delete[] eigenvector;
    if (m%10 == 0){
	ofstream writer(header.c_str());
    	if (writer){
	    writer << "Nuclear Charge" << endl;
            writer << Z << endl;
            writer << "Number of electrons" << endl;
            writer << n << endl;
            writer << "Basis set" << endl;
       	    for (i=0; i < good_coeffs.size(); i++){
        	writer << good_coeffs.at(i) << endl;
            }    	
        }
        writer.close();
    }	
    cout << "==========================================" << endl;
    }
    // ======================================================================================

    float* eigenvector = new float [m];
    float good_ev[m];

    coeffs.assign(good_coeffs.begin(), good_coeffs.end());

    gs_k = he_atom.calc_E(m,coeffs, eigenvector);
    for (j=0;j<m;j++){
        good_ev[j] = eigenvector[j];}

    cout << "========== Overall stochastic optimization started ==========" << endl;
    cout << "Ground state energy at iteration 1 is: " << gs_k << endl;

    //add random numbers to each coefficient in the basis set; only keep if energy improves.

    float cutoff;
    int iterations = 1000;

    for (i=2; i <iterations; i++){
        coeffs.assign(good_coeffs.begin(), good_coeffs.end());
        cutoff = exp(-i/iterations/2);
        vec_count = 0;
        for (basis_i = 0; basis_i < m; basis_i++){
            for (j = 0; j < n; j++) {
                for (k = j; k < n; k++){
                    if (good_ev[basis_i] < cutoff){
                        if (j == k){
                            rn = gauss_d(re);}
                        else {
                            rn = gauss_od(re);}
                        coeffs[vec_count] += rn;
                    }
                    vec_count += 1;
                }
            }
        }
        gs_kp1 = he_atom.calc_E(m,coeffs,eigenvector);

        if (gs_kp1 < gs_k){
            cout << "Ground state energy at iteration "<< i <<  " improved: " << gs_kp1 << endl;
            diff = abs(gs_kp1 - gs_k);
            if (diff > 10.0){
                cout << "Numerical instability: resetting energy!" << endl;
            }
            else{
                gs_k = gs_kp1;
                good_coeffs.assign(coeffs.begin(), coeffs.end());
                for (j=0;j<m;j++){
                    good_ev[j] = eigenvector[j];}
            }
        }
    }

    cout << "Basis function weights are: " << endl;
    for (i=0; i < m; i++){
        sum_c += good_ev[i]*good_ev[i];
        cout <<"function "<<setw(3) <<i<<" : " << setw(15) << good_ev[i]*good_ev[i] << endl;
    }
    cout << "Sum of weights: " << sum_c << endl;

    st_temp = "ne" + to_string(n) + "me" +to_string(m) + "_eigenvector.txt";
    ofstream c_out(st_temp.c_str());
        if (c_out){
            for (i=0; i < m; i++){
                c_out << good_ev[i] << endl;
            }
        cout << "Saved best linear coefficients to " << st_temp << endl;
        }
    c_out.close();

    ofstream writer(header.c_str());
        if (writer){
	    writer << "Nuclear Charge" << endl;
            writer << Z << endl;
            writer << "Number of electrons" << endl;
            writer << n << endl;
            writer << "Maximum Number of Basis functions" << endl;
            writer << max_basis << endl;
            writer << "Maximum Number of stochastic iterations per function" << endl;
            writer << max_iters << endl;
            writer << "Basis set" << endl;
            for (i=0; i < good_coeffs.size(); i++){
                writer << good_coeffs.at(i) << endl;
            }
        cout << "Saved best coefficients to " << header << endl;
        }
    writer.close();

    return 0;
}

string get_Header(int Z,int n, int m)
{
    string header;

    switch (Z){
        case 1: {
            string temp = to_string(n) ;
            string temp2 = to_string(m) ;
            header = "H_ne" + temp + "_me" + temp2;
            break;
        }
        case 2: {
            string temp = to_string(n) ;
            string temp2 = to_string(m) ;
            header = "He_ne" + temp + "_me" + temp2;
            break;
        }
        case 3: {
            string temp = to_string(n) ;
            string temp2 = to_string(m) ;
            header = "Li_ne" + temp + "_me" + temp2;
            break;
        }
        case 4: {
            string temp = to_string(n) ;
            string temp2 = to_string(m) ;
            header = "Be_ne" + temp + "_me" + temp2;
            break;
        }
        case 5: {
            string temp = to_string(n) ;
            string temp2 = to_string(m) ;
            header = "B_ne" + temp + "_me" + temp2;
            break;
        }
        default: {
            string temp = to_string(n) ;
            string temp2 = to_string(m) ;
            header = "Unk_atom_ne" + temp + "_me" + temp2;
            break;
        }
    }
    return header;
}

MatrixXf Energy::getMat(){

    MatrixXf S_total(4*24,4);

    S_total <<  1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,1.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,0.0,1.0,0.0,
                0.0,1.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,0.0,1.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,0.0,0.0,1.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                1.0,0.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,0.0,1.0,0.0,
                0.0,1.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,0.0,1.0,0.0,
                0.0,1.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,0.0,1.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,0.0,0.0,1.0,
                1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,0.0,1.0,0.0,
                0.0,1.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,0.0,1.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,1.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,1.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,0.0,0.0,1.0,
                1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,0.0,1.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,1.0,0.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,1.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                0.0,0.0,1.0,0.0,
                1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,0.0,1.0,
                1.0,0.0,0.0,0.0,
                0.0,1.0,0.0,0.0,
                0.0,0.0,1.0,0.0;
    return S_total;
}

