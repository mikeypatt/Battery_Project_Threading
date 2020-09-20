//
// Created by Michael Patterson on 20/09/2020.
//

#ifndef BATTERY_PROJECT_THREADING_DYNAMIC_OPTIMISATION_H
#define BATTERY_PROJECT_THREADING_DYNAMIC_OPTIMISATION_H
void *dynamic_sim (void *id);
void *dynamic_sim_post (void *id);

double dynamic_objective(unsigned n, const double *x, double *grad, void *my_func_data);

void dynamic_fit(std::vector<std::vector<double>> ***data,starting_Params* carryOvers[][no_datasets],
                 double* outputvoltage[][no_datasets],int index,double Reff_table[],
                 double Rct_table[],double C_table[],double R0_table[],double spoints_ref[]);


void dynamic_simulate(std::vector<std::vector<double>> *data,double output_voltage[],starting_Params* starting_data,double R0,double Reff,double Rct,
                      double C,double* Reff_table,double* Rct_table,
                      double* C_table,double* R0_table,double * spoints_ref,int index,starting_Params* ending_data = nullptr,bool final = false);

void dynamic_interpolation(std::vector<std::vector<double>> *data,double* spoints_ref,int index,double* points,double point,double* table);




#endif //BATTERY_PROJECT_THREADING_DYNAMIC_OPTIMISATION_H
