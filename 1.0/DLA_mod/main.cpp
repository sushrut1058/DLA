#include "ext.h"
//#include "funcs.cpp"
#include <bits/stdc++.h>
#define PI 3.14159265358979323846
#include <cmath>
#include <ctime>
#include <random>
#include <chrono>
#include <fstream>
//#define seed rand()

//const double pi=3.14159265358979323846;
using namespace std;



int main(){

    ofstream outnew;
	outnew.open ("DLA.csv");

	outnew << "0,1,2"<<endl;
	cout << setprecision(16);
    int ar = 5;
    double ln_ar=log(ar);
    double ln_2ar=log(2*ar);

    double R = 1;
    float s=2;
    int nop=10000;
    float r_max = 80;
    double i_orientations[nop][3];
    gen_ini_dir(nop,i_orientations);

    vector <vector<double>> cluster;
    vector <vector<double>> nl;

    double ir[3]={0,0,0}; double idir[3] = {i_orientations[0][0],i_orientations[0][1],i_orientations[0][2]}; double irod[ar][3];
    get_rod(irod,ir,idir,R,ar);
    for (int i=0;i<ar;i++){
        vector<double> vec {irod[i][0],irod[i][1],irod[i][2]};
        outnew <<vec[0]<<","<<vec[1]<<","<<vec[2]<<endl;
        cluster.push_back(vec);
    }

    double si_r_pe=s*sqrt(2.0/3);
    double si_r_pa=s*sqrt(1.0/3);
    double si_o =s/ln_ar*sqrt(3*(2.0*ln_2ar-1)/(ln_ar));
    double i_coords[3];
    double fac;
    double r_i[3];
    double par_i[3];
    double delr_per;
    double delr_par;
    double del_psi;
    double par[3];
    double mag;
    double r[3];
    double rod[ar][3];
    double overlaps[3];
    bool overlap,success,c_s,c_o;;
    double fin[ar][3];
    double foo;
    double unrot[ar][3];
    double per_i[3],delta[3];
    double tang[3];
    int nl_stepcount = 15; // for neighbourlist
    double limr = 5*R+nl_stepcount*s*R;
    int p=1;


    while (p<nop){
        for (int j=0;j<3;j++) {i_coords[j] = ((double)rand() / RAND_MAX)-0.5;}
        fac = r_max/sqrt(i_coords[0]*i_coords[0]+i_coords[1]*i_coords[1]+i_coords[2]*i_coords[2]);

        for (int j=0;j<3;j++) {
                r_i[j]=i_coords[j]*fac;
        }

        for (int j=0;j<3;j++) {par_i[j] = i_orientations[p][j];}
        int i=0;
        bool totbr=false;
        bool stuck=false;


        unsigned seedd = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seedd);
        normal_distribution<double> dist1(0,si_r_pe); //1
        normal_distribution<double> dist2(0,si_r_pa); //1
        normal_distribution<double> dist3(0,si_o); //1
        while (i<300000){
        //#random walk

            overlap=false;//2
            success=false;//2
            perp_dir_angle(per_i,par_i);//2
            delr_per=dist1(generator);//2
            delr_par=dist2(generator);//2
            del_psi=dist3(generator);//2

            for (int j=0;j<3;j++){
                    delta[j]=delr_par*par_i[j]+delr_per*per_i[j];
                    r[j]= r_i[j]+ delta[j];
            }

            cross(par_i,per_i,tang);

            for (int j=0;j<3;j++) {par[j]=par_i[j]+del_psi*tang[j];}
            mag = sqrt(par[0]*par[0]+par[1]*par[1]+par[2]*par[2]);
            for (int j=0;j<3;j++) par[j]/=mag;

            if (sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])>r_max+200) break;

            if (i%nl_stepcount==0){
                nl.clear();
                neighbour_list(nl,r,cluster,limr);
            }

            if (nl.size() != 0){
                get_rod(rod,r,par,R,ar);
                get_overlaps(nl, overlaps, rod, overlap, success, ar);
            }else{
                overlap = false;
                success = false;
            }

            if (success==true && overlap == false){
                for (int j=0;j<ar;j++){
                    for(int k=0;k<3;k++){
                                fin[j][k] = rod[j][k];
                            }
                }
                stuck=true;
            }


            if (overlap==true){
                get_rod(unrot,r,par_i,R,ar);
                get_overlaps(nl, overlaps, unrot, c_o, c_s, ar);

                while (1){
                    if (c_o==false and c_s==true){
                        success=true;
                        overlap=false;
                        cout << "Perfection after undoing rotation"<<endl;
                        break;
                    }
                    else if (c_o==false and c_s==false){
                        int jj = overlaps[2],ii=overlaps[1];

                        double cl_p[3] = {nl[jj][0],nl[jj][1],nl[jj][2]};
                        double r_p[3] = {unrot[ii][0],unrot[ii][1],unrot[ii][2]};

                        double a=0,b=0,c=2*R;

                        for (int j=0;j<3;j++){
                            a+=(r_p[j]-r[j])*(r_p[j]-r[j]);b+=(cl_p[j]-r[j])*(cl_p[j]-r[j]);
                        }
                        a = sqrt(a);b=sqrt(b);


                        int f1=1,f2=1;
                        double p=0,q=0;

                        for (int j=0;j<3;j++){
                            p+=(rod[ii][j]-r_p[j])*(tang[j]);
                            q+=(r_p[j]-r[j])*(par_i[j]);
                        }

                        if(p<0) f2=-1;
                        if(q<0) f1=-1;

                        double cosine = (a*a+b*b-c*c)/(2*a*b);

                        double A=0,B=0;
                        for(int j=0;j<3;j++){
                            A += f1*(cl_p[j]-r[j])*par_i[j]; B+=f2*(cl_p[j]-r[j])*tang[j];
                        }

                        double phi = asin(A/sqrt(A*A+B*B));
                        double rhs = (cosine*b)/sqrt(A*A+B*B);

                        double theta = asin(rhs) - phi;

                        double coords[3],br=0,axis[3],new_rod[ar][3];

                        for (int j=0;j<3;j++){
                            coords[j]=a*(cos(theta)*f1*par_i[j] + sin(theta)*tang[j]*f2)+r[j];
                            br += (coords[j]-cl_p[j])*(coords[j]-cl_p[j]);
                            axis[j] = coords[j]-r[j];

                        }

                        br=sqrt(br);

                        get_rod(new_rod,r,axis,R,ar);


                        double os[3];
                        bool cs,o;
                        get_overlaps(nl, os, new_rod, o, cs, ar);


                        if (o==false){
                            if (br>2.1*R){
                                break;
                            }
                            stuck=true;

                            for (int j=0;j<ar;j++){
                                    for(int k=0;k<3;k++){
                                        fin[j][k] = new_rod[j][k];
                                    }
                            }
                            c_s = true;
                            c_o = false;
                            break;
                        }
                        else{

                            for (int j=0;j<3;j++) overlaps[j] = os[j];
                            for (int j=0;j<ar;j++){
                                for(int k=0;k<3;k++){
                                    rod[j][k] = new_rod[j][k];
                                }
                            }
                        }
                    }else{

                        double untrans[ar][3];
                        get_rod(untrans,r_i,par_i,R,ar);

                        double os[3]; bool o,c_o,c_s;
                        get_overlaps(nl, os, unrot, c_o, c_s, ar);

                        int rp = os[1], cp = os[2];

                        double x[3],y[3];
                        for (int j=0;j<3;j++){
                            x[j] = untrans[rp][j];
                            y[j] = nl[cp][j];
                        }

                        double diff = frac_distance(delta, x, y ,R);
                        double perf_r[3];
                        double mag =sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
                        for (int j=0;j<3;j++){perf_r[j] = r_i[j] + diff * delta[j]/mag;}
                        double testrod[ar][3];
                        bool suc;

                        get_rod(testrod,perf_r,par_i,R,ar);
                        get_overlaps(nl, os, testrod, o, suc, ar);

                        if (o==false){
                            for (int j=0;j<ar;j++) r[j] = perf_r[j];
                            stuck = true;

                            for (int j=0;j<ar;j++){
                                    for(int k=0;k<3;k++){
                                        fin[j][k] = testrod[j][k];
                                }
                            }

                            break;
                        }
                        else{
                            if (i==0){
                                cout << "Need to increase the firing radius"<<endl;
                                r_max+=10;
                                break;
                            }

                            for (int j=0;j<ar;j++){ for(int k=0;k<3;k++){unrot[j][k] = testrod[j][k];}}
                        }

                    }
                }
            }

            if (stuck==true){
                int l = cluster.size()/5;
                cout << "Appending particle, #id="<< l+1<<" "<<r_max<<endl;
                for (int j=0;j<ar;j++){
                    vector<double> v {fin[j][0],fin[j][1],fin[j][2]};
                    cluster.push_back(v);
                    outnew << setprecision(16);
                    outnew << v[0]<<","<<v[1]<<","<<v[2]<<endl;
                }
                r_max+=.375;
                if(cluster.size()%200==0){
                    r_max+=5;
                }
                break;
            }

            for (int j=0;j<3;j++) {par_i[j]=par[j];}
            for (int j=0;j<3;j++) {r_i[j]=r[j];}
            i+=1;
            if (totbr==true){
                break;
            }

        }

        if (totbr==true){
            break;
        }
        p+=1;
    }
    cout << "done";
    return 0;
}
//2000p 26min21secs
//2000p 34min33secs
//2000p 30min59secs
//2000p 27min47secs
//2000p 29min12secs
//2680p 5230 secs 1hr20min
//3405p 9194 secs 2hr33min

/*
[[47.71964710596988, 2685],
 [45.08129249356216, 2010],
 [47.26558792146913, 2010],
 [45.80278564377659, 2010],
 [45.61608748141891, 1995]]
 */
