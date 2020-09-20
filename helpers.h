//
// Created by Michael Patterson on 27/06/2020.
//
#ifndef IMPERIAL_BATTERY_HELPERS_H
#define IMPERIAL_BATTERY_HELPERS_H
const int no_datasets = 16;
const double OCV_Coef[9] = {162.819296346443,-626.424040821280,994.123599474504,-838.370905010509,395.859371472140,-94.4306297230054,4.31313232297881,3.37833790202489,2.92273089870673};
const double R0_Coef[9] = {0.428373375339626,-3.09155625471716, 8.49955501406088,-12.1797208311519,10.1834668524870,-5.16594525582102,1.57736163925228,-0.271762771100002,0.0460650738087307};

double interpolate(std::vector<std::vector<double>> *Data, double x, bool extrapolate,int collumn);
double interpolate(std::vector<std::vector<double>> *Data,double y[], double x, bool extrapolate);
double interpolate(double xdata[],double ydata[], double x,int size, bool extrapolate);

double polyfit(double point,const double coeff[]);
void polyfit(double points[],double y[],int length,const double coeff[]);
void polyfit(std::vector<std::vector<double>> *data,double y[],int length,const double coeff[],int x_col);

std::vector<std::vector<double>> **dataSplitter(std::vector<std::vector<double>> data, double split_points[],int charge_points_length);
double* splitting_points_ref();
double* splitting_points_opt();
double* splitting_points_new();

double L2_norm_distance(double calculated_voltage[],std::vector<std::vector<double>> *measured_data,int length);
double Trapz(std::vector<std::vector<double>> *measured_data,double calculated_voltage[]);


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



double OCV_Con(unsigned n, const double *x, double *grad, void *my_func_data){

    optData *d = (optData*) my_func_data;
    double ans = 0;
    double inter;
    for(int i=0;i<no_datasets;i++) {
        starting_Params *param = d->carryOvers[i];
        std::vector<std::vector<double>> ***relevantDataPtr = d->data;
        double voltage = (*(relevantDataPtr[i])[d->index])[1][2];
        inter = voltage - static_sim_mini(relevantDataPtr[i][d->index], param,x[3], x[0], x[1], x[2]);
        ans += pow(inter, 2.0);
    }

    return sqrt(ans);
}

double startingV_con_1(unsigned n, const double *x, double *grad, void *my_func_data){

    optData *d = (optData*) my_func_data;
    double ans=0;
    for(int i=0;i<no_datasets;i++) {
        starting_Params *param = d->carryOvers[i];
        std::vector<std::vector<double>> ***relevantDataPtr = d->data;
        double current = (*(relevantDataPtr[i])[d->index])[0][1];
        double voltage = (*(relevantDataPtr[i])[d->index])[0][2];
        ans+= pow(voltage - (x[4+i] - current * x[3]),2.0);
    }

    return sqrt(ans);

}

#endif //IMPERIAL_BATTERY_HELPERS_H
