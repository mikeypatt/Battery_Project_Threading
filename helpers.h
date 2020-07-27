//
// Created by Michael Patterson on 27/06/2020.
//
#ifndef IMPERIAL_BATTERY_HELPERS_H
#define IMPERIAL_BATTERY_HELPERS_H

double interpolate(std::vector<std::vector<double>> *Data, double x, bool extrapolate,int collumn);
double interpolate(std::vector<std::vector<double>> *Data,double y[], double x, bool extrapolate);
double interpolate(double xdata[],double ydata[], double x,int size, bool extrapolate);

double polyfit(double point,const double coeff[]);
void polyfit(double points[],double y[],int length,const double coeff[]);
void polyfit(std::vector<std::vector<double>> *data,double y[],int length,const double coeff[],int x_col);

std::vector<std::vector<double>> **dataSplitter(std::vector<std::vector<double>> data, double split_points[],int charge_points_length);
double* splitting_points_ref();
double* splitting_points_opt();

double L2_norm_distance(double calculated_voltage[],std::vector<std::vector<double>> *measured_data,int length);
double Trapz(std::vector<std::vector<double>> *measured_data,double calculated_voltage[]);

#endif //IMPERIAL_BATTERY_HELPERS_H
