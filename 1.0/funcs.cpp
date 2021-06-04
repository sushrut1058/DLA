#include <bits/stdc++.h>
#define PI 3.14159265358979323846
#include <cmath>
#include <ctime>
#define seed rand()
const double pi = PI;

using namespace std;

double metropolis(double th[], int N){

    double th_o=pi/2;
    double th_c;
    double al;

    th[0]=th_o;

    for (int i=1;i<N;i++){

            th_c = pi * ((double)rand() / RAND_MAX);
            al = sin(th_c)/sin(th[i-1]);

            if (((double)rand() / RAND_MAX)< min(1.0,al)){
                th[i] = th_c;
            }else{
                th[i] = th[i-1];
            }
            //cout << th[i] << ",";
    }
}

double gen_ini_dir(int np,double f[][3]){
    double pi = PI;
  double th_o[np];
  metropolis(th_o,np);

  double phi_o[np];
  for (int i=0;i<np;i++){
        phi_o[i]=2*pi*((double)rand()/RAND_MAX);
  }

    double c_al[np];
    double c_be[np];
    double c_ga[np];

   for (int p=0;p<np;p++)
    {
        c_ga[p]=cos(th_o[p]);
        if((double)rand()/RAND_MAX<0.5){
        c_be[p]=-sqrt(tan(phi_o[p])*tan(phi_o[p])*(1-cos(th_o[p])*cos(th_o[p]))*(cos(phi_o[p])*cos(phi_o[p])));
      }else{
        c_be[p]=sqrt(tan(phi_o[p])*tan(phi_o[p])*(1-cos(th_o[p])*cos(th_o[p]))*(cos(phi_o[p])*cos(phi_o[p])));
      }
      //Need to resort to the sum of squares as tan(phi) takes bigger values
      if((double)rand()/RAND_MAX<0.5){
        c_al[p]=sqrt(abs(1-c_be[p]*c_be[p]-c_ga[p]*c_ga[p]));
      }else{
        c_al[p]=-sqrt(abs(1-c_be[p]*c_be[p]-c_ga[p]*c_ga[p]));
      }
    }


    for (int i=0;i<np;i++){
        f[i][0]=c_al[i];
        f[i][1]=c_be[i];
        f[i][2]=c_ga[i];
    }
}


double perp_dir_angle(double perp[],double par[]){
    double pi = M_PI;
  double th=acos(par[2]);
  double phi=atan(par[1]/par[0]);
  int n=0;
  if(par[0]<0 && par[1]>=0){
      n=1;
    }
  else if (par[0]<0 && par[1]<0){
    n=1;
  }
  else if (par[0]>=0 && par[1]<=0){
    n=2;
  }
  phi=phi+n*pi;

  if (abs(par[0])<pow(10,-8)){
      if (abs(par[1])<pow(10,-8)){
        phi=(double)rand()/RAND_MAX*2*pi;
      }
    }
    double al=(double)rand()/RAND_MAX*2*pi;

  perp[0]=cos(phi)*cos(th)*cos(al)-sin(phi)*sin(al);
  perp[1]=sin(phi)*cos(th)*cos(al)+cos(phi)*sin(al);
  perp[2]=-sin(th)*cos(al);
    //cout <<"[" <<perp[0] <<","<<perp[1] <<","<<perp[2] <<"],"<<endl;
}


double get_rod(double rod[][3],double r[],double direction[],double R,int ar){
    int j=0;
    double mag = sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
    for (int i=0;i<3;i++){
        direction[i] = direction[i]/mag;
    }
    for (int k=-int(ar/2);k<=(ar/2);k++){
        if (ar%2!=0){
            for (int i=0;i<3;i++){
                rod[j][i] = r[i]-2*R*(k)*direction[i];
            }
            j+=1;
        }else{
            if (k==ar/2) {break;}
            for (int i=0;i<3;i++){
                rod[j][i] = r[i]+2*(k/2)*R*(pow(-1,k))*direction[i];
            }
            j+=1;
        }
    }
}


double frac_distance(double direction[], double c1[], double c2[],double R){
    //dist = ((c1-c2)[0]**2+(c1-c2)[1]**2+(c1-c2)[2]**2)**.5
    double vec[3] = {c1[0]-c2[0],c1[1]-c2[1],c1[2]-c2[2]};
    double dist = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    double mag = sqrt(direction[0]*direction[0]+direction[1]*direction[1]+direction[2]*direction[2]);
    double foo = -(vec[0]*direction[0] + vec[1]*direction[1] + vec[2]*direction[2])/mag;
    double perp = dist*sin(acos(foo/dist));
    double theta = asin(perp/2*R);
    double phi = asin(perp/dist);
    double diff = foo - 2.00000000*R*cos(theta) ;
    return diff;
}

double cross(double a[], double b[], double c[]){
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}


double get_overlaps(vector<vector<double>>cluster,double overlaps[] ,double rod[][3],bool &overlap,bool &success,int ar){
    double R=1.0;
    double d;
    double min_d = 1000;
    double v1[3],v2[3],d1,d2,dot1=0,dot2=0;
    overlap = false;
    success = false;
    for (int j=0;j<cluster.size();j++){
        for (int p=0;p<3;p++){
            v1[p] = cluster[j][p]-rod[0][p];
            v2[p] = cluster[j][p]-rod[ar-1][p];
        }

        dot1 = v1[0]*(v1[0]-v2[0])+v1[1]*(v1[1]-v2[1])+v1[2]*(v1[2]-v2[2]);
        dot2 = v2[0]*(v1[0]-v2[0])+v2[1]*(v1[1]-v2[1])+v2[2]*(v1[2]-v2[2]);
        d1 = sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
        d2 = sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);

        if(dot1*dot2>0){
            if (min(d1,d2)<=2.01*R){
                for(int i=0;i<ar;i++){
                    d = sqrt((cluster[j][0]-rod[i][0])*(cluster[j][0]-rod[i][0])+(cluster[j][1]-rod[i][1])*(cluster[j][1]-rod[i][1])+(cluster[j][2]-rod[i][2])*(cluster[j][2]-rod[i][2]));
                    if (d<2.01*R){
                        success = true;
                        if (d<1.99999999*R){
                            overlap = true;
                            if (min_d>d){
                                overlaps[0] = d;
                                overlaps[1] = i;
                                overlaps[2] = j;
                                min_d = d;
                            }
                        }
                    }
                }

            }
        }else{
            if (min(d1,d2)>2.01*R){
                double theta = acos(dot1/(d1*sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2]))));
                if (d1*sin(theta)<=2.01*R){
                    for(int i=0;i<ar;i++){
                        d = sqrt((cluster[j][0]-rod[i][0])*(cluster[j][0]-rod[i][0])+(cluster[j][1]-rod[i][1])*(cluster[j][1]-rod[i][1])+(cluster[j][2]-rod[i][2])*(cluster[j][2]-rod[i][2]));
                        if (d<2.01*R){
                            success = true;
                            if (d<1.99999999*R){
                                overlap = true;
                                if (min_d>d){
                                    overlaps[0] = d;
                                    overlaps[1] = i;
                                    overlaps[2] = j;

                                    min_d = d;
                                }
                            }
                        }
                    }
                }
            }else{// #the minimum of the distances of the ends is less than 2.01R
              //  cout << 32;
                for(int i=0;i<ar;i++){
                        d = sqrt((cluster[j][0]-rod[i][0])*(cluster[j][0]-rod[i][0])+(cluster[j][1]-rod[i][1])*(cluster[j][1]-rod[i][1])+(cluster[j][2]-rod[i][2])*(cluster[j][2]-rod[i][2]));
                        if (d<2.01*R){
                            success = true;
                            if (d<1.99999999*R){
                                overlap = true;
                                if (min_d>d){
                                    overlaps[0] = d;
                                    overlaps[1] = i;
                                    overlaps[2] = j;

                                    min_d = d;
                                }
                            }
                        }
                    }
                }
            }
        }
}


double neighbour_list(vector <vector<double>> &nl, double r[3], vector <vector<double>> cluster, double limr){
    double d;
    int i=0;

    for (int j=0;j<cluster.size();j++){ // iterate through the cluster
        //d = distance between the COM and a cluster particle
        d = sqrt((cluster[j][0]-r[0])*(cluster[j][0]-r[0]) + (cluster[j][1]-r[1])*(cluster[j][1]-r[1]) + (cluster[j][2]-r[2])*(cluster[j][2]-r[2]));
        //m = min(m,d);

        if (d<=limr){ // distance within the range
            //append it to neighbour list
            vector<double> v {cluster[j][0],cluster[j][1],cluster[j][2]};
            nl.push_back(v);
            i++;
        }
    }
}


/*
double get_overlaps(vector<vector<double>>cluster,double overlaps[] ,double rod[][3],bool &overlap,bool &success,int ar){

    int R=1;
    int k=0;
    double d;
    double min_d = 1000;
    overlap = false;
    success = false;
    for (int j=0;j<cluster.size();j++){
        for(int i=0;i<ar;i++){
            d = sqrt((cluster[j][0]-rod[i][0])*(cluster[j][0]-rod[i][0])+(cluster[j][1]-rod[i][1])*(cluster[j][1]-rod[i][1])+(cluster[j][2]-rod[i][2])*(cluster[j][2]-rod[i][2]));
            //cout << d <<endl;
            //overlaps[0] = d;
            if (d<2.01*R){
                success = true;
                if (d<1.99999999*R){
                    overlap = true;
                    if (min_d>d){
                        overlaps[0] = d;
                        overlaps[1] = i;
                        overlaps[2] = j;

                        min_d = d;
                    }
                    k++;
                }
            }
        }
    }
}*/
