//
// Created by Michael Patterson on 27/06/2020.
//

#include<cmath>
#include <vector>
#include "helpers.h"

using namespace std;

double interpolate(vector<vector<double>> *Data, double x, bool extrapolate,int collumn)
{
    int size = (*Data).size();
    double starting_time = (*Data)[0][0];
    int i = 0;                                                                  // find left end of interval for interpolation
    if ( x >= (*Data)[size - 2][0]-starting_time)                                                 // special case: beyond right end
    {
        i = size - 2;
    }
    else
    {
        while ( x > (*Data)[i+1][0]-starting_time) i++;
    }

    double xL = (*Data)[i][0]-starting_time, yL = (*Data)[i][collumn], xR = (*Data)[i+1][0]-starting_time, yR = (*Data)[i+1][collumn];
//    cout << xL << " " << yL << " "<< xR << " " << yR << endl; // points on either side (unless beyond ends)
    if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }

    double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

    return yL + dydx * ( x - xL );                                              // linear interpolation
}

double interpolate(vector<vector<double>> *Data,double y[], double x, bool extrapolate)
{
    int size = (*Data).size();
    double starting_time = (*Data)[0][0];

    int i = 0;                                                                  // find left end of interval for interpolation
    if ( x >= (*Data)[size - 2][0]-starting_time)                                                 // special case: beyond right end
    {
        i = size - 2;
    }
    else
    {
        while ( x > (*Data)[i+1][0]-starting_time) i++;
    }
    double xL = (*Data)[i][0]-starting_time, yL = y[i], xR = (*Data)[i+1][0]-starting_time, yR = y[i+1]; // points on either side (unless beyond ends)
    if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }

    double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

    return yL + dydx * ( x - xL );                                              // linear interpolation
}

double interpolate(double xdata[],double ydata[], double x,int size, bool extrapolate)
{
    int i = 0;                                                                  // find left end of interval for interpolation
    if ( x >= xdata[size - 2] )                                                 // special case: beyond right end
    {
        i = size - 2;
    }
    else
    {
        while ( x > xdata[i+1] ) i++;
    }
    double xL = xdata[i], yL = ydata[i], xR = xdata[i+1], yR = ydata[i+1]; // points on either side (unless beyond ends)
    if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
    {
        if ( x < xL ) yR = yL;
        if ( x > xR ) yL = yR;
    }

    double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

    return yL + dydx * ( x - xL );                                              // linear interpolation
}

double polyfit(double point,const double coeff[]){
    return pow(point,8)* coeff[0] + pow(point,7)* coeff[1] + pow(point,6)* coeff[2]
           + pow(point,5)* coeff[3] + pow(point,4)* coeff[4] + pow(point,3)* coeff[5]
           + pow(point,2)* coeff[6] + pow(point,1)* coeff[7] +coeff[8];
}

void polyfit(double points[],double y[],int length,const double coeff[]){
    for(int i=0;i<length;i++){
        y[i] = pow(points[i],8)* coeff[0] + pow(points[i],7)* coeff[1] + pow(points[i],6)* coeff[2]
               + pow(points[i],5)* coeff[3] + pow(points[i],4)* coeff[4] + pow(points[i],3)* coeff[5]
               + pow(points[i],2)* coeff[6] + pow(points[i],1)*coeff[7] +coeff[8];
    }
}

void polyfit(vector<vector<double>> *data,double y[],int length,const double coeff[],int x_col){
    for(int i=0;i<length;i++){
        y[i] = pow((*data)[i][x_col],8)* coeff[0] + pow((*data)[i][x_col],7)* coeff[1] + pow((*data)[i][x_col],6)* coeff[2]
               + pow((*data)[i][x_col],5)* coeff[3] + pow((*data)[i][x_col],4)* coeff[4] + pow((*data)[i][x_col],3)* coeff[5]
               + pow((*data)[i][x_col],2)* coeff[6] + pow((*data)[i][x_col],1)*coeff[7] +coeff[8];
    }

}



//vector<vector<double>> **dataSplitter(vector<vector<double>> data,double split_points[],int charge_points_length){
//
//    double minimum = 0,maximum =0;
//
//    vector<vector<double>> **datapointers = new vector<vector<double>>*[charge_points_length];
//
//
//    for(int i= 0;i<charge_points_length;i++){
//        if(i==0){
//            minimum = split_points[i];
//        }
//        else{
//            minimum = split_points[i-1];
//        }
//        if(i ==(charge_points_length-1)){
//            maximum = split_points[i];
//        }
//        else{
//            maximum = split_points[i+1];
//        }
//
//        vector<vector<double>> *vect = new vector<vector<double>>;
//
//        for(int j=0;j<data.size();j++){
//
//            if (data[j][3] > minimum && data[j][3] <= maximum){
//                vector<double> row{data[j][0],data[j][1],data[j][2],data[j][3]};
//                vect->push_back(row);
//            }
//        }
//        datapointers[i] = vect;
//    }
//
//    return datapointers;
//}


double L2_norm_distance(double calculated_voltage[],vector<vector<double>> *measured_data,int length){
    double sum = 0;
    for(int i=0 ;i<length;i++){
        sum += pow( (*measured_data)[i][2]-calculated_voltage[i],2);
    }
    return sqrt(sum);
}

double Trapz(vector<vector<double>> *measured_data,double calculated_voltage[]){

    int length = measured_data->size();
    double dt;
    double ans=0;
    for(int i=0;i<length-1;i++){
        dt = (*measured_data)[i+1][0] - (*measured_data)[i][0];
        double diff_0 = pow( (*measured_data)[i][2]-calculated_voltage[i],2);
        double diff_1 = pow( (*measured_data)[i+1][2]-calculated_voltage[i+1],2);
        ans += ((diff_0+diff_1)/2.0) * dt;
    }
    return sqrt(ans);
}


