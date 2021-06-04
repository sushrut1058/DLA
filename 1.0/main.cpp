#include "ext.h"
#include "funcs.cpp" // please comment this line if you have codeblocks and you opened DLA.cbp to run this code
#include <bits/stdc++.h>
#define PI 3.14159265358979323846
#include <cmath>
#include <ctime>
#include <random>
#include <chrono>
#include <fstream>

using namespace std;



int main(){

    ofstream outnew;
	outnew.open ("DLA.csv");
	outnew << "0,1,2"<<endl;
	cout << setprecision(16);

    int ar = 5; // aspect ratio
    double ln_ar=log(ar);
    double ln_2ar=log(2*ar);
    double R = 1; //radius (remove)
    float s=2; //stepsize
    int nop=10000; //no. of rods to fire (not the number of rods that will form the cluster)
    float r_max = 80; // firing radius
    float rej = 500; // rejection distance ( ||r|| > r_max+rej => reject particle)
    double i_orientations[nop][3]; //initial orientations (remove)

    gen_ini_dir(nop,i_orientations); // assigning

    vector <vector<double>> cluster; // final cluster
    vector <vector<double>> nl; // neighbour list

    double ir[3]={0,0,0}; //seed COM
    double idir[3] = {i_orientations[0][0],i_orientations[0][1],i_orientations[0][2]}; // seed orientation
    double irod[ar][3]; //seed
    get_rod(irod,ir,idir,R,ar); // assigning

    for (int i=0;i<ar;i++){
        vector<double> vec {irod[i][0],irod[i][1],irod[i][2]};
        outnew <<vec[0]<<","<<vec[1]<<","<<vec[2]<<endl;
        cluster.push_back(vec);
    }

    double si_r_pe=s*sqrt(2.0/3); //perpendicular std deviation
    double si_r_pa=s*sqrt(1.0/3); // parallel std deviation
    double si_o =s/ln_ar*sqrt(3*(2.0*ln_2ar-1)/(ln_ar)); //angular std deviation
    double i_coords[3]; //initial rod COM coordinates (on the firing radius)
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
    bool overlap,success,c_s,c_o;
    double fin[ar][3]; //final rod that will be appended to the cluster
    double foo;
    double unrot[ar][3];
    double per_i[3],delta[3];
    double tang[3];
    int nl_stepcount = 30; // for neighbourlist
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

        bool stuck=false;


        unsigned seedd = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seedd);
        normal_distribution<double> dist1(0,si_r_pe); // n0ormal dist.
        normal_distribution<double> dist2(0,si_r_pa); // " "
        normal_distribution<double> dist3(0,si_o); // " "
        while (i<100000){
        //#random walk

            overlap=false;
            success=false;

            perp_dir_angle(per_i,par_i);// generating perpendicular vector

            delr_per=dist1(generator);// random point from the distribution
            delr_par=dist2(generator);// " "
            del_psi=dist3(generator);// " "

            for (int j=0;j<3;j++){
                    delta[j]=delr_par*par_i[j]+delr_per*per_i[j];
                    r[j]= r_i[j]+ delta[j];
            }

            cross(par_i,per_i,tang); // generating a tangent to the rod end

            for (int j=0;j<3;j++) {par[j]=par_i[j]+del_psi*tang[j];} //rotation
            mag = sqrt(par[0]*par[0]+par[1]*par[1]+par[2]*par[2]);
            for (int j=0;j<3;j++) par[j]/=mag; // unit vector

            if (sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])>r_max+rej) break; // reject a particle outside this range

            if (i%nl_stepcount==0){ // update neighbour list after nl_stepcount iterations
                nl.clear(); // clear neighbour list
                neighbour_list(nl,r,cluster,limr);
            }

            if (nl.size() != 0){// if the neighbour list is populated
                get_rod(rod,r,par,R,ar); //get rod
                get_overlaps(nl, overlaps, rod, overlap, success, ar); //check overlaps
            }else{ //no overlaps as no particles nearby
                overlap = false;
                success = false; //success=true means perfect sticking
            }

            if (success==true && overlap == false){
                for (int j=0;j<ar;j++){
                    for(int k=0;k<3;k++){
                                fin[j][k] = rod[j][k]; //final rod
                            }
                }
                stuck=true; //stuck
            }

            /*dealing with overlaps*/
            if (overlap==true){

                get_rod(unrot,r,par_i,R,ar); //undo rotation
                get_overlaps(nl, overlaps, unrot, c_o, c_s, ar); //get overlaps

                while (1){
                    if (c_o==false and c_s==true){ // if rod perfectly sticks after undoing rotation
                        success=true;
                        overlap=false;
                        for (int j=0;j<ar;j++){
                            for(int k=0;k<3;k++){
                                fin[j][k] = unrot[j][k];
                            }
                        }
                        stuck=true; //stuck
                        break;
                    }
                    /*moving a fractional angle*/
                    else if (c_o==false and c_s==false){

                        int jj = overlaps[2],ii=overlaps[1]; //(overlapping) jj = cluster particle, ii=rod particle
                        /*coordinates*/
                        double cl_p[3] = {nl[jj][0],nl[jj][1],nl[jj][2]};
                        double r_p[3] = {unrot[ii][0],unrot[ii][1],unrot[ii][2]};


                        double a=0,b=0,c=2*R;

                        for (int j=0;j<3;j++){
                            a+=(r_p[j]-r[j])*(r_p[j]-r[j]);b+=(cl_p[j]-r[j])*(cl_p[j]-r[j]);
                        }
                        a = sqrt(a);b=sqrt(b);

                        /* moving by a fractional angle*/
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

                        double os[3]; // new overlaps in this position
                        bool cs,o; // booleans for the same
                        get_overlaps(nl, os, new_rod, o, cs, ar);


                        if (o==false){ // stuck
                            if (br>2.01*R){ //for debugging purposes
                                cout << "some error";
                                break;
                            }
                            stuck=true; //stuck

                            for (int j=0;j<ar;j++){
                                    for(int k=0;k<3;k++){
                                        fin[j][k] = new_rod[j][k]; //final rod
                                    }
                            }
                            c_s = true;
                            c_o = false;
                            break;
                        }
                        else{ //another overlap in the new position

                            for (int j=0;j<3;j++) overlaps[j] = os[j]; //overlaps = os(new overlaps)
                            for (int j=0;j<ar;j++){
                                for(int k=0;k<3;k++){
                                    rod[j][k] = new_rod[j][k]; //final rod
                                }
                            }
                        }
                    }
                    /*fractional translation*/
                    else{

                        double untrans[ar][3];
                        get_rod(untrans,r_i,par_i,R,ar); // undoing translation

                        double os[3]; bool o,c_o,c_s;

                        get_overlaps(nl, os, unrot, c_o, c_s, ar); //get overlaps

                        int rp = os[1], cp = os[2];

                        double x[3],y[3];
                        for (int j=0;j<3;j++){
                            x[j] = untrans[rp][j];
                            y[j] = nl[cp][j];
                        }

                        double diff = frac_distance(delta, x, y ,R); // fraction of a step to move

                        double perf_r[3];
                        double mag =sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
                        for (int j=0;j<3;j++){perf_r[j] = r_i[j] + diff * delta[j]/mag;} //corrected position

                        /*test rod again for overlaps in new position*/
                        double testrod[ar][3];
                        bool suc;
                        get_rod(testrod,perf_r,par_i,R,ar);
                        get_overlaps(nl, os, testrod, o, suc, ar);

                        if (o==false){ // no overlaps

                            stuck = true;//stuck
                            for (int j=0;j<ar;j++){
                                    for(int k=0;k<3;k++){
                                        fin[j][k] = testrod[j][k];
                                }
                            }
                            break;
                        }
                        else{ // if overlaps in the new position

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

            /*appending particle to a cluster*/
            if (stuck==true){
                int l = cluster.size()/5;
                cout << "Appending particle, #id="<< l+1<<" "<<r_max<<endl;

                for (int j=0;j<ar;j++){
                    vector<double> v {fin[j][0],fin[j][1],fin[j][2]};
                    cluster.push_back(v);
                    outnew << setprecision(16);
                    outnew << v[0]<<","<<v[1]<<","<<v[2]<<endl;
                }

                r_max+=.21; // 0.21 is just a random number, no good reason for this
                break;
            }

            for (int j=0;j<3;j++) {
                    par_i[j]=par[j];
                    r_i[j]=r[j];
            }

            i+=1;


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

