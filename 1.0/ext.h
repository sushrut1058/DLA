#pragma once
#include <vector>
using namespace std;

double metropolis(double[], int);
double gen_ini_dir(int,double[][3]);
double perp_dir_angle(double[] ,double[] );
double get_rod(double[][3],double[],double[],double,int);
double get_overlaps(vector <vector<double>>,double[],double[][3],bool &,bool &,int);
double frac_distance(double[],double[],double[],double);
double cross(double[], double[], double[]);
double neighbour_list(vector <vector<double>> &, double[], vector <vector<double>>, double);
