//
// Created by Michael Patterson on 30/06/2020.
//
#ifndef IMPERIAL_BATTERY_OPTIMISATIONS_H
#define IMPERIAL_BATTERY_OPTIMISATIONS_H

const int no_datasets = 16;
const double OCV_Coef[9] = {162.819296346443,-626.424040821280,994.123599474504,-838.370905010509,395.859371472140,-94.4306297230054,4.31313232297881,3.37833790202489,2.92273089870673};
const double R0_Coef[9] = {0.428373375339626,-3.09155625471716, 8.49955501406088,-12.1797208311519,10.1834668524870,-5.16594525582102,1.57736163925228,-0.271762771100002,0.0460650738087307};


struct starting_Params{
    double Vc;
    double starting_charge1;
    double starting_charge2;
    double I2;
    double I2_bar;
    starting_Params(double Vc,double starting_charge1,double starting_charge2,double I2,double I2_bar):
            Vc(Vc),starting_charge1(starting_charge1),starting_charge2(starting_charge2),I2(I2),I2_bar(I2_bar){}
};


double static_sim_mini(std::vector<std::vector<double>> *data,starting_Params* starting_data,double R0,double Reff,double Rct,double C);



struct Thread_Static{

    int id;
    std::vector<std::vector<double>> *data_set;
    double* output_voltage;
    starting_Params* starting_data;
    double R0;
    double Reff;
    double Rct;
    double C;
    double* answer;
    int length;

    starting_Params* ending_data;
    bool final;

    Thread_Static(int id,std::vector<std::vector<double>> *data_set,double* output_voltage,starting_Params* starting_data,
                  double R0,double Reff,double Rct,double C,int length,double* answer = nullptr,starting_Params* ending_data = nullptr, bool final = false):
            id(id),data_set(data_set),output_voltage(output_voltage),starting_data(starting_data),R0(R0),Reff(Reff),Rct(Rct),C(C),answer(answer),length(length),
            ending_data(ending_data),final(final) {}

};

struct Thread_Dynamic{

    int id;
    std::vector<std::vector<double>> *data_set;
    double* output_voltage;
    starting_Params* starting_data;
    double R0;
    double Reff;
    double Rct;
    double C;
    double* answer;
    int length;
    double* Reff_table;
    double* Rct_table;
    double* C_table;
    double* R0_table;
    double* s_points_opt;
    int index;

    starting_Params* ending_data;
    bool final;

    Thread_Dynamic(int id,std::vector<std::vector<double>> *data_set,double* output_voltage,starting_Params* starting_data,
                   double R0,double Reff,double Rct,double C,int length,double* Reff_table,double* Rct_table,
                   double* C_table,double* R0_table,double* s_points_opt,int index,
                   double* answer = nullptr,starting_Params* ending_data = nullptr, bool final = false):
            id(id),data_set(data_set),output_voltage(output_voltage),starting_data(starting_data),R0(R0),Reff(Reff),Rct(Rct),C(C),answer(answer),
            length(length),Reff_table(Reff_table),Rct_table(Rct_table),C_table(C_table),R0_table(R0_table),s_points_opt(s_points_opt),
            index(index),ending_data(ending_data),final(final) {}

};

struct optData{
    std::vector<std::vector<double>> *** data;
    int index;
    starting_Params** carryOvers;
    double* Reff_table = nullptr;
    double* Rct_table = nullptr;
    double* C_table = nullptr;
    double * R0_table = nullptr;
    double* s_points_opt = nullptr;

    optData(std::vector<std::vector<double>>*** data,int index,starting_Params** carryOvers,double *Reff_table=nullptr,double *Rct_table=nullptr,double *C_table=nullptr,double *R0_table=nullptr,double *s_points_opt=nullptr):
            data(data),index(index),carryOvers(carryOvers),Reff_table(Reff_table),Rct_table(Rct_table),
            C_table(C_table),R0_table(R0_table),s_points_opt(s_points_opt){}
};

double startingV_con_1(unsigned n, const double *x, double *grad, void *my_func_data);
double OCV_Con(unsigned n, const double *x, double *grad, void *my_func_data);


void static_fit(std::vector<std::vector<double>> ***data,double carry_over_parameters[],
                starting_Params* carryOvers[][no_datasets],double* outputvoltage[][no_datasets],int index,
                double Reff_table[],double Rct_table[],double C_table[],double R0_table[]);

void *static_sim (void *id);
void *static_sim_post (void *id);
double static_objective(unsigned n, const double *x, double *grad, void *my_func_data);

void static_simulate(std::vector<std::vector<double>> *data,double output_voltage[],starting_Params* starting_data,double R0,double Reff,
                     double Rct,double C,starting_Params* ending_data = nullptr,bool final = false);


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


#endif IMPERIAL_BATTERY_OPTIMISATIONS_H
