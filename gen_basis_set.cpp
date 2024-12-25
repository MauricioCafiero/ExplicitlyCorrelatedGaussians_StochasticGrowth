#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <random>
#include <bits/stdc++.h>

using namespace std;

default_random_engine& global_urng(){
    static default_random_engine u{};
    return u;
}

string get_Header(int, int);

int main()
{
    vector<float> coeffs(0);
    float rn;
    string header;
    int i,j,k,m,n,basis_i;

    cout << "Generate all coefficients (matrix elements) for n=2-4 explicitly correlated Gaussian functions: " << endl;

    //get filename with atom type, electrons and basis functions
    cout << "Enter number of electrons: " << endl;
    cin >> n;
    cout << "Enter number of basis functions: " << endl;
    cin >> m;

    cout << "Read in " << m << " basis functions." << endl;
    cout << "Read in number of electrons: " << n << endl;

    header = get_Header(n,m);

    default_random_engine re(random_device{}());
//    uniform_real_distribution<float> gauss_d(0.0,0.05);
//    uniform_real_distribution<float> gauss_od(0.0,0.08);
    normal_distribution<float> gauss_d(1.0,2.5);
    normal_distribution<float> gauss_od(0.01,0.05);

    for (basis_i = 0; basis_i < m; basis_i++){
        for (j = 0; j < n; j++) {
            for (k = j; k < n; k++){
                if (j == k){
                    rn = gauss_d(re);}
                else {
                    rn = gauss_od(re);}
                coeffs.push_back(rn);
            }
        }
    }

    header += "_coeffs.txt";
    ofstream writer(header.c_str());
        if (writer){
            writer << n << endl;
            for (i=0; i < coeffs.size(); i++){
                writer << coeffs.at(i) << endl;
            }
        cout << "Saved coefficients to " << header << endl;
        }
    writer.close();

    return 0;
}

string get_Header(int n, int m)
{
    string header;

    string temp = to_string(n) ;
    string temp2 = to_string(m) ;
    header = "New_ne" + temp + "_me" + temp2;

    return header;
}
