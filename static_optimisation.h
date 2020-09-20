//
// Created by Michael Patterson on 30/06/2020.
//
#ifndef IMPERIAL_BATTERY_STATIC_SIMULATE_H
#define IMPERIAL_BATTERY_STATIC_SIMULATE_H

void static_fit(std::vector<std::vector<double>> ***data,double carry_over_parameters[],
                starting_Params* carryOvers[][no_datasets],double* outputvoltage[][no_datasets],int index,
                double Reff_table[],double Rct_table[],double C_table[],double R0_table[]);

void *static_sim (void *id);
void *static_sim_post (void *id);
double static_objective(unsigned n, const double *x, double *grad, void *my_func_data);

void static_simulate(std::vector<std::vector<double>> *data,double output_voltage[],starting_Params* starting_data,double R0,double Reff,
                     double Rct,double C,starting_Params* ending_data = nullptr,bool final = false);


#endif //IMPERIAL_BATTERY_STATIC_SIMULATE_H
