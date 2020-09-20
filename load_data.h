#ifndef IMPERIAL_BATTERY_LOAD_DATA_H
#define IMPERIAL_BATTERY_LOAD_DATA_H

void readInData(const char* filename,std::vector<std::vector<double>>& data);
void cumTrapz(double batteryCap,std::vector<std::vector<double>>& data);

std::vector<std::vector<double>> **dataSplitter(std::vector<std::vector<double>> data, double split_points[],int charge_points_length);
double* splitting_points_ref();
double* splitting_points_opt();
double* splitting_points_new();

#endif //IMPERIAL_BATTERY_LOAD_DATA_H
