#include <iostream>
#include <cstring>
#include <vector>
#include<cmath>
#include <nlopt.hpp>
#include <fstream>
#include <random>
#include "load_data.h"
#include "helpers.h"
#include "implicitRK4.h"

using namespace std;

const int no_datasets = 16;
const double OCV_Coef[9] = {162.819296346443,-626.424040821280,994.123599474504,-838.370905010509,395.859371472140,-94.4306297230054,4.31313232297881,3.37833790202489,2.92273089870673};

struct starting_Params{
    double Vc;
    double starting_charge1;
    double starting_charge2;
    double I2;
    double I2_bar;
    starting_Params(double Vc,double starting_charge1,double starting_charge2,double I2,double I2_bar):
    Vc(Vc),starting_charge1(starting_charge1),starting_charge2(starting_charge2),I2(I2),I2_bar(I2_bar){}
};


struct Thread_Static{

    int id;
    vector<vector<double>> *data_set;
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

    Thread_Static(int id,vector<vector<double>> *data_set,double* output_voltage,starting_Params* starting_data,
           double R0,double Reff,double Rct,double C,int length,double* answer = nullptr,starting_Params* ending_data = nullptr, bool final = false):
            id(id),data_set(data_set),output_voltage(output_voltage),starting_data(starting_data),R0(R0),Reff(Reff),Rct(Rct),C(C),answer(answer),length(length),
            ending_data(ending_data),final(final) {}

};

struct Thread_Dynamic{

    int id;
    vector<vector<double>> *data_set;
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

    Thread_Dynamic(int id,vector<vector<double>> *data_set,double* output_voltage,starting_Params* starting_data,
                  double R0,double Reff,double Rct,double C,int length,double* Reff_table,double* Rct_table,
                  double* C_table,double* R0_table,double* s_points_opt,int index,
                  double* answer = nullptr,starting_Params* ending_data = nullptr, bool final = false):
            id(id),data_set(data_set),output_voltage(output_voltage),starting_data(starting_data),R0(R0),Reff(Reff),Rct(Rct),C(C),answer(answer),
            length(length),Reff_table(Reff_table),Rct_table(Rct_table),C_table(C_table),R0_table(R0_table),s_points_opt(s_points_opt),
            index(index),ending_data(ending_data),final(final) {}

};


void *static_sim (void *id);
void *static_sim_post (void *id);

void *dynamic_sim (void *id);
void *dynamic_sim_post (void *id);

typedef pair<double,double> pairs;
double static_objective(unsigned n, const double *x, double *grad, void *my_func_data);


void semiStatic_fit(vector<vector<double>> ***data,double carry_over_parameters[],
        starting_Params* carryOvers[][no_datasets],
        double* outputvoltage[][no_datasets],int index,double Reff_table[],double Rct_table[],
        double C_table[],double R0_table[]);

void static_simulate(std::vector<std::vector<double>> *data,double output_voltage[],starting_Params* starting_data,double R0,double Reff,
double Rct,double C,starting_Params* ending_data = nullptr,bool final = false);


double dynamic_objective(unsigned n, const double *x, double *grad, void *my_func_data);

void dynamic_fit(vector<vector<double>> ***data,starting_Params* carryOvers[][no_datasets],
        double* outputvoltage[][no_datasets],int index,double Reff_table[],
        double Rct_table[],double C_table[],double R0_table[],double spoints_ref[]);


void dynamic_simulate(vector<vector<double>> *data,double output_voltage[],starting_Params* starting_data,double R0,double Reff,double Rct,
                      double C,double* Reff_table,double* Rct_table,
                      double* C_table,double* R0_table,double * spoints_ref,int index,starting_Params* ending_data = nullptr,bool final = false);

void dynamic_interpolation(vector<vector<double>> *data,double* spoints_ref,int index,double* points,double point,double* table);



struct optData{
    vector<vector<double>> *** data;
    int index;
    starting_Params** carryOvers;
    double* Reff_table = nullptr;
    double* Rct_table = nullptr;
    double* C_table = nullptr;
    double * R0_table = nullptr;
    double* s_points_opt = nullptr;

    optData(vector<vector<double>>*** data,int index,starting_Params** carryOvers,double *Reff_table=nullptr,double *Rct_table=nullptr,double *C_table=nullptr,double *R0_table=nullptr,double *s_points_opt=nullptr):
            data(data),index(index),carryOvers(carryOvers),Reff_table(Reff_table),Rct_table(Rct_table),
            C_table(C_table),R0_table(R0_table),s_points_opt(s_points_opt){}
};

int main() {

    vector<vector<double>> data02ConstCurrent;
    vector<vector<double>> data02LongPulse;
    vector<vector<double>> data02ShortPulse;
    vector<vector<double>> data03ConstCurrent;
    vector<vector<double>> data03LongPulse;
    vector<vector<double>> data03ShortPulse;
    vector<vector<double>> data04ConstCurrent;
    vector<vector<double>> data04LongPulse;
    vector<vector<double>> data04ShortPulse;
    vector<vector<double>> data05ConstCurrent;
    vector<vector<double>> data05LongPulse;
    vector<vector<double>> data05ShortPulse;

    vector<vector<double>> data1CPulse;
    vector<vector<double>> dataDriveCycleLong;
    vector<vector<double>> dataDriveCycleShort;
    vector<vector<double>> dataPulses;

    const char* CC_02_CC = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.2_const_current.txt" ;
    readInData(CC_02_CC,data02ConstCurrent);
    data02ConstCurrent.pop_back();
    const char* CC_02_Long = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.2_long_pulse.txt" ;
    readInData(CC_02_Long,data02LongPulse);
    data02LongPulse.pop_back();
    const char* CC_02_Short = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.2_short_pulse.txt" ;
    readInData(CC_02_Short,data02ShortPulse);
    data02ShortPulse.pop_back();

    const char* CC_03_CC = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.3CC.txt" ;
    readInData(CC_03_CC,data03ConstCurrent);
    data03ConstCurrent.pop_back();
    const char* CC_03_Long = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.3_long_pulse.txt" ;
    readInData(CC_03_Long,data03LongPulse);
    data03LongPulse.pop_back();
    const char* CC_03_Short = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.3_short_pulse.txt" ;
    readInData(CC_03_Short,data03ShortPulse);
    data03ShortPulse.pop_back();

    const char* CC_04_CC = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.4CC.txt" ;
    readInData(CC_04_CC,data04ConstCurrent);
    data04ConstCurrent.pop_back();
    const char* CC_04_Long = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.4_long_pulse.txt" ;
    readInData(CC_04_Long,data04LongPulse);
    data04LongPulse.pop_back();
    const char* CC_04_Short = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.4_short_pulse.txt" ;
    readInData(CC_04_Short,data04ShortPulse);
    data04ShortPulse.pop_back();

    const char* CC_05_CC = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.5CC.txt" ;
    readInData(CC_05_CC,data05ConstCurrent);
    data05ConstCurrent.pop_back();
    const char* CC_05_Long = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.5_long_pulse.txt" ;
    readInData(CC_05_Long,data05LongPulse);
    data05LongPulse.pop_back();
    const char* CC_05_Short = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_0.5_short_pulse.txt" ;
    readInData(CC_05_Short,data05ShortPulse);
    data05ShortPulse.pop_back();

    const char* CC_1_pulses = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_1C_pulse.txt" ;
    readInData(CC_1_pulses,data1CPulse);
    data1CPulse.pop_back();
    const char* CC_drive_long = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_drive_cycles_long.txt" ;
    readInData(CC_drive_long,dataDriveCycleLong);
    dataDriveCycleLong.pop_back();
    const char* CC_drive_short = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_drive_cycles_short.txt" ;
    readInData(CC_drive_short,dataDriveCycleShort);
    dataDriveCycleShort.pop_back();
    const char* CC_pulses = "/Users/michaelpatterson/CLionProjects/Imperial_Battery/raw_data/data_pulses.txt" ;
    readInData(CC_pulses,dataPulses);
    dataPulses.pop_back();

    double* spoints_ref = splitting_points_ref();

    int charge_point_length = 41;
    //split the data into the different opt bits
    vector<vector<double>> ** CC_02_CC_data;
    CC_02_CC_data = dataSplitter(data02ConstCurrent,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_02_Long_data;
    CC_02_Long_data = dataSplitter(data02LongPulse,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_02_short_data;
    CC_02_short_data = dataSplitter(data02ShortPulse,spoints_ref,charge_point_length);

    //split the data into the different opt bits
    vector<vector<double>> ** CC_03_CC_data;
    CC_03_CC_data = dataSplitter(data03ConstCurrent,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_03_Long_data;
    CC_03_Long_data = dataSplitter(data03LongPulse,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_03_short_data;
    CC_03_short_data = dataSplitter(data03ShortPulse,spoints_ref,charge_point_length);

    //split the data into the different opt bits
    vector<vector<double>> ** CC_04_CC_data;
    CC_04_CC_data = dataSplitter(data04ConstCurrent,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_04_Long_data;
    CC_04_Long_data = dataSplitter(data04LongPulse,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_04_short_data;
    CC_04_short_data = dataSplitter(data05ShortPulse,spoints_ref,charge_point_length);

    //split the data into the different opt bits
    vector<vector<double>> ** CC_05_CC_data;
    CC_05_CC_data = dataSplitter(data05ConstCurrent,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_05_Long_data;
    CC_05_Long_data = dataSplitter(data05LongPulse,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_05_short_data;
    CC_05_short_data = dataSplitter(data05ShortPulse,spoints_ref,charge_point_length);

    vector<vector<double>> ** CC_1C_pulse_data;
    CC_1C_pulse_data = dataSplitter(data1CPulse,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_Drive_Long_data;
    CC_Drive_Long_data = dataSplitter(dataDriveCycleLong,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_Drive_Short_data;
    CC_Drive_Short_data = dataSplitter(dataDriveCycleShort,spoints_ref,charge_point_length);
    vector<vector<double>> ** CC_Pulses_data;
    CC_Pulses_data = dataSplitter(dataPulses,spoints_ref,charge_point_length);

    //clears the data to save some memory!
    data02ConstCurrent.clear();
    data02LongPulse.clear();
    data02ShortPulse.clear();
    data03ConstCurrent.clear();
    data03LongPulse.clear();
    data03ShortPulse.clear();
    data04ConstCurrent.clear();
    data04LongPulse.clear();
    data04ShortPulse.clear();
    data05ConstCurrent.clear();
    data05LongPulse.clear();
    data05ShortPulse.clear();
    data1CPulse.clear();
    dataDriveCycleLong.clear();
    dataDriveCycleShort.clear();
    dataPulses.clear();

    // I think this is a vector pointer that points to an array of vector pointer pointers of size number of datasets
    vector<vector<double>> ***fitting_data = new vector<vector<double>>**[no_datasets];
    fitting_data[0] = CC_02_CC_data; // TODO change thissss back
    fitting_data[1] = CC_02_Long_data;
    fitting_data[2] = CC_02_short_data;

    fitting_data[3] = CC_03_CC_data;
    fitting_data[4] = CC_03_Long_data;
    fitting_data[5] = CC_03_short_data;
    fitting_data[6] = CC_04_CC_data;
    fitting_data[7] = CC_04_Long_data;
    fitting_data[8] = CC_04_short_data;
    fitting_data[9] = CC_05_CC_data;
    fitting_data[10] = CC_05_Long_data;
    fitting_data[11] = CC_05_short_data;

    fitting_data[12] = CC_1C_pulse_data;
    fitting_data[13] = CC_Drive_Long_data;
    fitting_data[14] = CC_Drive_Short_data;
    fitting_data[15] = CC_Pulses_data;

    double Reff_table[charge_point_length];
    double Rct_table[charge_point_length];
    double R0_table[charge_point_length];
    double C_table[charge_point_length];
    double alpha_table[charge_point_length];

    //coefficeicnts of OCV and RO need to work out how to fit this
    double OCV_Coef[9] = {162.819296346443,-626.424040821280,994.123599474504,-838.370905010509,395.859371472140,-94.4306297230054,4.31313232297881,3.37833790202489,2.92273089870673};

    // creates an struct for the starting values
    starting_Params* starting_params[41][no_datasets];

    double carry_over_parameters[4] = {0.02,0.02,0.02,0.02};
    double* outputVoltage[charge_point_length][no_datasets];

    //reverse iterates from a high soc to a lower one
    for(int i=40;i>=0;i--){
        semiStatic_fit(fitting_data,carry_over_parameters,starting_params,outputVoltage,
                       i,Reff_table, Rct_table, C_table, R0_table);
    }

    // delete all the heap dataaaaaa
    ofstream out_stream;

    out_stream.open("semistatic_tables.txt");
    for(int i=0;i<41;i++){
        out_stream << R0_table[i] << " " << Reff_table[i] << " " << Rct_table[i] << " " << C_table[i] << endl;
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_02C_short.txt");
    for(int i=0;i<41;i++){
        int length = (*(fitting_data[0])[i]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[0])[i])[j][2] << " " << outputVoltage[i][0][j] << endl;
        }
    }

    cout << "finished semistatic" << endl;

    //dynamic bit
    double* outputVoltagedynamic[41][no_datasets];
    int order[41];
    for(int i=0;i<41;i++)
        order[i] = i;
    int dynamic_index = 0;

    //reverse iterates from a high soc to a lower one
    auto rng = std::default_random_engine {};

    for(int i=0;i<10;i++) {

        cout << "iter: " << i << endl;
        std::shuffle(std::begin(order), std::end(order), rng);

        for(int j=0;j<41;j++){
            dynamic_index = order[j];
            dynamic_fit(fitting_data,starting_params,outputVoltagedynamic,dynamic_index,
                    Reff_table, Rct_table, C_table, R0_table,spoints_ref);
        }

        out_stream.open("dynamic_tables.txt");
        for(int j=0;j<41;j++){
            out_stream << Reff_table[j] << " " << Rct_table[j] << " " << C_table[j] << endl;
        }
        out_stream.close();


        out_stream.open("dynamic_voltage_02C_CC.txt");
        for(int ii=0;i<41;i++){
            int length = (*(fitting_data[0])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[0])[ii])[j][2] << " " << outputVoltagedynamic[ii][0][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_02C_long.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[1])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[1])[ii])[j][2] << " " << outputVoltagedynamic[ii][1][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_02C_short.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[2])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[2])[ii])[j][2] << " " << outputVoltagedynamic[ii][2][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_03C_CC.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[3])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[3])[ii])[j][2] << " " << outputVoltagedynamic[ii][3][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_03C_long.txt");
        for(int ii=0;i<41;i++){
            int length = (*(fitting_data[4])[i]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[4])[i])[j][2] << " " << outputVoltagedynamic[i][4][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_03C_short.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[5])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[5])[ii])[j][2] << " " << outputVoltagedynamic[ii][5][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_04C_CC.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[6])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[6])[ii])[j][2] << " " << outputVoltagedynamic[ii][6][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_04C_long.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[7])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[7])[ii])[j][2] << " " << outputVoltagedynamic[ii][7][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_04C_short.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[8])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[8])[ii])[j][2] << " " << outputVoltagedynamic[ii][8][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_05C_CC.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[9])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[9])[ii])[j][2] << " " << outputVoltagedynamic[ii][9][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_05C_long.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[10])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[10])[ii])[j][2] << " " << outputVoltagedynamic[ii][10][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_05C_short.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[11])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[11])[ii])[j][2] << " " << outputVoltagedynamic[ii][11][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_1C_pulse.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[12])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[12])[ii])[j][2] << " " << outputVoltagedynamic[ii][12][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_drive_long.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[13])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[13])[ii])[j][2] << " " << outputVoltagedynamic[ii][13][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_drive_short.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[14])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[14])[ii])[j][2] << " " << outputVoltagedynamic[ii][14][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_pulses.txt");
        for(int ii=0;ii<41;ii++){
            int length = (*(fitting_data[15])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[15])[ii])[j][2] << " " << outputVoltagedynamic[ii][15][j] << endl;
            }
        }
        out_stream.close();

        //delete output voltages
        for(int j =0;j<41;j++){
            for(int ii=0;ii<no_datasets;ii++){
                delete[] outputVoltage[j][ii];
            }
        }
        
    }


    delete spoints_ref;

    return 0;
}

void dynamic_fit(vector<vector<double>> ***data,starting_Params* carryOvers[][no_datasets],
        double* outputvoltage[][no_datasets],int index,double Reff_table[],double Rct_table[],
        double C_table[],double R0_table[],double spoints_ref[]){

    const int number_of_parameters = 4;

    double lb[number_of_parameters];
    double ub[number_of_parameters];
    lb[0] = 0.001; // Reff lb
    lb[1] = 0.001;// Rct lb
    lb[2] = 0.001;// C lb
    lb[3] = 0.001; //R0


    ub[0] = 1; // Reff ub
    ub[1] = 1; // Rct ub
    ub[2] = 1;  // C ub
    ub[3] = 1;  // C ub

    double x[number_of_parameters];
    x[0] = Reff_table[index];
    x[1] = Rct_table[index];
    x[2] = C_table[index];
    x[3] = R0_table[index];

    //create addtional datastruct
    optData addData(data,index,carryOvers[index],Reff_table,Rct_table,C_table,R0_table,spoints_ref);

    //creates optimsation object
    nlopt_opt opt; // NLOPT_LD_SLSQP
    opt = nlopt_create(NLOPT_LN_COBYLA,number_of_parameters);

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    nlopt_set_stopval(opt,-HUGE_VAL);

    nlopt_set_ftol_abs(opt, 1e-30);
    nlopt_set_xtol_rel(opt, 1e-30);
    nlopt_set_maxeval(opt,1);

    nlopt_set_min_objective(opt,dynamic_objective,&addData);
    double minf;
    if (nlopt_optimize(opt, x, &minf) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found minimum at f(%g,%g,%g) = %0.10g\n", x[0],x[1],x[2],minf);
    }
    nlopt_destroy(opt);

    //filling up lookup tables with values
    Reff_table[index] = x[0];
    Rct_table[index] = x[1];
    C_table[index] = x[2];
    R0_table[index] = x[3];

    vector<vector<double>>* relevant_data_set = nullptr;

    double S;
    if (index == 40 ||index == 39) {
        for (int i = 0; i<no_datasets;i++) {
            relevant_data_set = (data[i])[index];
            carryOvers[index][i]->Vc = (*relevant_data_set)[0][2] + x[3] * (*relevant_data_set)[0][1];
            S = (*relevant_data_set)[0][3];
            initial0CV(carryOvers[index][i]->Vc,OCV_Coef,S, x[3], x[1], (*relevant_data_set)[0][1]);
            carryOvers[index][i]->starting_charge1 = S;
            carryOvers[index][i]->starting_charge2 = S;
        }
    }

    //Creating the simulation Threads
    Thread_Dynamic* simulating_threads[no_datasets];
    pthread_t simulationid[no_datasets];
    int sim_id;

    // filling up the output voltage data and starting charges for the dataset (index-2)number_of_datasets
    for(int i = 0;i<no_datasets;i++){

        relevant_data_set = (data[i])[index];

        if(index>1) {
            carryOvers[index - 2][i] = new starting_Params(0, 0, 0, 0, 0);
            sim_id = i+1;
            simulating_threads[i] = new Thread_Dynamic(sim_id,relevant_data_set,outputvoltage[index][i],carryOvers[index][i],x[3], x[0], x[1], x[2],
                                                      relevant_data_set->size(),Reff_table,Rct_table,C_table,R0_table,spoints_ref,index,nullptr,carryOvers[index - 2][i],true);
            pthread_create (&simulationid[i],NULL,static_sim_post,(void*) simulating_threads[i]);

        }
        else {
            sim_id = i + 1;
            simulating_threads[i] = new Thread_Dynamic(sim_id,relevant_data_set,outputvoltage[index][i],carryOvers[index][i],x[3], x[0], x[1], x[2],
                                                       relevant_data_set->size(),Reff_table,Rct_table,C_table,R0_table,spoints_ref,index);
            pthread_create(&simulationid[i], NULL, static_sim_post, (void *) simulating_threads[i]);
        }
    }

    for(auto & i : simulationid)
    {
        pthread_join (i, NULL) ;
    }

    for(int i=0;i<no_datasets;i++)
    {
        delete simulating_threads[i];
    }

}

void *dynamic_sim_post (void *producer_parameter)
{

    //Intialise the structure passed in to thread function
    Thread_Dynamic* d = (Thread_Dynamic*) producer_parameter;
    vector<vector<double>> *data = d->data_set;
    double* output_voltage = d->output_voltage;
    starting_Params* starting_data = d->starting_data;
    starting_Params* ending_data = d->ending_data;
    bool final = d->final;
    double R0 = d->R0;
    double Reff = d->Reff;
    double Rct = d->Rct;
    double C = d->C;
    int length = d->length;

    output_voltage = new double[length];

    dynamic_simulate(data,output_voltage,starting_data,R0,Reff,Rct,C,d->Reff_table,d->Rct_table,d->C_table,d->R0_table,
            d->s_points_opt,d->index,ending_data,final);

    pthread_exit(0);

}

double dynamic_objective(unsigned n, const double *x, double *grad, void *my_func_data){

    optData *d = (optData*) my_func_data;

    double Reff = x[0];
    double Rct = x[1];
    double C = x[2];
    double R0 = x[3];

    vector<vector<double>> ***relevantDataPtr = d->data;
    vector<vector<double>>* relevant_data_set = nullptr;

    double S;
    if (d->index == 40 ||d->index == 39) {
        for (int i = 0; i < no_datasets; i++) {
            relevant_data_set = (relevantDataPtr[i])[d->index];
            d->carryOvers[i]->Vc = (*relevant_data_set)[0][2] + R0 * (*relevant_data_set)[0][1];
            S = (*relevant_data_set)[0][3];
            initial0CV(d->carryOvers[i]->Vc, OCV_Coef, S, R0, Rct, (*relevant_data_set)[0][1]);
            d->carryOvers[i]->starting_charge1 = S;
            d->carryOvers[i]->starting_charge2 = S;
        }
    }

    double *output_voltage  = nullptr;
    double answer = 0;
    size_t length=0;

    //Creating the simulation Threads
    Thread_Dynamic* simulating_threads[no_datasets];
    pthread_t simulationid[no_datasets];
    int sim_id;

    //to put each threads answer in;
    double answers[no_datasets];

    for (int i=0;i<no_datasets;i++){

        relevant_data_set = (relevantDataPtr[i])[d->index];
        length = relevant_data_set->size();

        sim_id = i+1;
        simulating_threads[i] = new Thread_Dynamic(sim_id,relevant_data_set,output_voltage,d->carryOvers[i],R0,Reff,Rct,
                C,length,d->Reff_table,d->Rct_table,d->C_table,d->R0_table,d->s_points_opt,d->index,&answers[i]);

        pthread_create (&simulationid[i],NULL,dynamic_sim,(void*) simulating_threads[i]);
    }

    for(auto & i : simulationid)
    {
        pthread_join (i, NULL) ;
    }

    for(int i=0;i<no_datasets;i++)
    {
        delete simulating_threads[i];
        answer += answers[i];
    }

    return answer;

}

void *dynamic_sim (void *producer_parameter)
{

    //Intialise the structure passed in to thread function
    Thread_Dynamic* d = (Thread_Dynamic *) producer_parameter;

    vector<vector<double>> *data = d->data_set;
    double* output_voltage = d->output_voltage;
    starting_Params* starting_data = d->starting_data;
    double R0 = d->R0;
    double Reff = d->Reff;
    double Rct = d->Rct;
    double C = d->C;
    int length = d->length;

    output_voltage = new double[length];

    double* answer = d->answer;

    dynamic_simulate(data,output_voltage,starting_data,R0,Reff,Rct,C,d->Reff_table,d->Rct_table,d->C_table,d->R0_table,d->s_points_opt,
            d->index);

    *answer = L2_norm_distance(output_voltage,data,length);

    //stops the memory from leaking for the output voltage array
    delete[] output_voltage;

    pthread_exit(0);

}


void dynamic_interpolation(vector<vector<double>> *data,double* spoints_ref,int index,double* points,double point,double* table){
    table[index] = point;
    for(int i=0;i<data->size();i++){
        points[i] = interpolate(spoints_ref,table,(*data)[i][3],41,true);
    }
}


void dynamic_simulate(vector<vector<double>> *data,double output_voltage[],starting_Params* starting_data,double R0,double Reff,double Rct,
                      double C,double* Reff_table,double* Rct_table,
                      double* C_table,double* R0_table,double * spoints_ref,int index,starting_Params* ending_data,bool final){

    int size = (int)(*data).size();

    double starting_time = (*data)[0][0];

    // interpolation functions
    double* _Reff = new double[size];
    double* _Rct= new double[size];
    double* _C= new double[size];
    double* _R0= new double[size];

    dynamic_interpolation(data,spoints_ref,index,_Reff,Reff,Reff_table);
    dynamic_interpolation(data,spoints_ref,index,_Rct,Rct,Rct_table);
    dynamic_interpolation(data,spoints_ref,index,_C,C,C_table);
    dynamic_interpolation(data,spoints_ref,index,_R0,R0,R0_table);

    vector<double>* new_data = new vector<double>;
    new_data->push_back((*data)[0][0]);
    double ittime=0.0;
    for(size_t i=1;i<data->size();i++) {
        if (((*data)[i][0] - (*data)[i-1][0])  > 0.1) {
            ittime = (*data)[i - 1][0];
            //while time is smaller than the next one
            while (ittime < (*data)[i][0]) {
                ittime += 0.1;
                new_data->push_back(ittime);
            }
        }else{
            new_data->push_back((*data)[i][0]);
        }
    }

    int length = new_data->size();
    double dt[length-1];

    double* time = new double[length];
    double* I = new double[length];
    double* dt_= new double[length];

    double* Reff_ = new double[length];
    double* Rct_= new double[length];
    double* C_= new double[length];
    double* R0_= new double[length];

    for(int i=0;i<length;i++){
        time[i] = (*new_data)[i] -(*new_data)[0];
        I[i] = interpolate(data, time[i], true, 1);
        Reff_[i] = interpolate(data,_Reff,time[i] , true);
        Rct_[i] = interpolate(data,_Rct,time[i] , true);
        C_[i] = interpolate(data,_C,time[i] , true);
        R0_[i] = interpolate(data,_R0,time[i] , true);
        if(i>0){
            dt_[i-1] = time[i] - time[i-1];
        }
    }

    //save some memory
    delete new_data;
    delete[] _Reff;
    delete[] _Rct;
    delete[] _C;
    delete[] _R0;

    for(int i=1;i<length;i++){
        dt[i-1] = time[i] - time[i-1];
    }

    // now for the simulation
    double* V = new double[length];
    double* Vc = new double[length];
    double* I2_bar=new double[length];
    double* I2=new double[length];
    double* OCV1=new double[length];
    double* OCV2=new double[length];
    double* S1=new double[length];
    double* S2=new double[length];

    double Q = 5.0 * 3600;

    Vc[0] = starting_data->Vc;
    V[0] = Vc[0] - I[0] * R0_[0];

    S1[0] = starting_data->starting_charge1;
    S2[0] =  starting_data->starting_charge2;

    OCV1[0] = polyfit(S1[0],OCV_Coef);
    OCV2[0] = polyfit(S2[0],OCV_Coef);

    I2[0] = starting_data->I2;
    I2_bar[0] = starting_data->I2_bar;

    double odeI[2];
    for(size_t i=1;i<length;i++){

        //calculating stuff from the previous timestep
        S2[i] = S2[i - 1] - dt[i-1] / (Q/2.0) * I2_bar[i - 1];
        S1[i] = S1[i - 1] - dt[i-1]/(Q/2.0) * I2[i - 1];

        OCV1[i] = polyfit(S1[i],OCV_Coef);
        OCV2[i] = polyfit(S2[i],OCV_Coef);

        odeI[0] = I[i - 1];
        odeI[1] = I[i];
        Vc[i] = implicitRK4(time[i-1], Vc[i-1], dt[i-1],odeI, &OCV1[i - 1], &C_[i-1], &Rct_[i-1],false);

        I2_bar[i] = (OCV2[i] - OCV1[i]) / Reff;

        if(OCV1[i] == OCV2[i])
            I2[i] = I[i];
        else
            I2[i] = ((OCV1[i] - Vc[i]) / Rct) - ((OCV2[i] - OCV1[i]) / Reff);

        V[i] = Vc[i] - I[i] * R0;
    }

    for(int i =0;i<data->size();i++)
        output_voltage[i] = interpolate(time,V,(*data)[i][0]-starting_time,length,true);

    if(final){
        ending_data->Vc = Vc[length-1];
        ending_data->starting_charge1 = S1[length-1];
        ending_data->starting_charge2 = S2[length-1];
        ending_data->I2 = I2[length-1];
        ending_data->I2_bar = I2_bar[length-1];
    }

    // now for the simulation
    delete[] V;
    delete[] Vc;
    delete[] I2_bar;
    delete[] I2;
    delete[] OCV1;
    delete[] OCV2;
    delete[] S1;
    delete[] S2;

    delete[] Reff_;
    delete[] Rct_;
    delete[] C_;
    delete[] R0_;

    delete[] I;
    delete[] time;
    delete[] dt_;

}


void semiStatic_fit(vector<vector<double>> ***data,double carry_over_parameters[],
        starting_Params* carryOvers[][no_datasets],double* outputvoltage[][no_datasets],int index,
        double Reff_table[],double Rct_table[],double C_table[],double R0_table[]){

    const int number_of_parameters = 4;

    double lb[number_of_parameters];
    double ub[number_of_parameters];

    lb[0] = 0.001; // Reff lb
    lb[1] = 0.001;// Rct lb
    lb[2] = 0.001;// C lb
    lb[3] = 0.001;

    ub[0] = 1.0; // Reff ub
    ub[1] = 1.0; // Rct ub
    ub[2] = 1.0;  //C ub
    ub[3] = 1.0; //R0 ub

    double x [number_of_parameters];
    x[0] = carry_over_parameters[0];
    x[1] = carry_over_parameters[1];
    x[2] = carry_over_parameters[2];
    x[3] = carry_over_parameters[3];

    //create addtional datastruct
    optData addData(data,index,carryOvers[index]);

    //creates optimsation object
    nlopt_opt opt; // NLOPT_LD_SLSQP
    opt = nlopt_create(NLOPT_LN_COBYLA , number_of_parameters);

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    nlopt_set_stopval(opt,-HUGE_VAL);

    nlopt_set_ftol_abs(opt, 1e-30);
    nlopt_set_xtol_rel(opt, 1e-30);

    nlopt_set_maxeval(opt,1);

    nlopt_set_min_objective(opt,static_objective,&addData);
    double minf;
    if (nlopt_optimize(opt, x, &minf) < 0) {
        printf("nlopt failed!\n");
    }
    else {
        printf("found minimum at f(%g,%g,%g) = %0.10g\n", x[0],x[1],x[2],minf);
    }
    nlopt_destroy(opt);

    //filling up lookup tables with values
    Reff_table[index] = x[0];
    Rct_table[index] = x[1];
    C_table[index] = x[2];
    R0_table[index] = x[3];

    carry_over_parameters[0] = x[0];
    carry_over_parameters[1] = x[1];
    carry_over_parameters[2] = x[2];
    carry_over_parameters[3] = x[3];


    vector<vector<double>>* relevant_data_set = nullptr;

    double S;
    if (index == 40 ||index == 39) {
        for (int i = 0; i<no_datasets;i++) {
            relevant_data_set = (data[i])[index];
            carryOvers[index][i]->Vc = (*relevant_data_set)[0][2] + x[3] * (*relevant_data_set)[0][1];
            S = (*relevant_data_set)[0][3];
            initial0CV(carryOvers[index][i]->Vc,OCV_Coef,S, x[3], x[1], (*relevant_data_set)[0][1]);
            carryOvers[index][i]->starting_charge1 = S;
            carryOvers[index][i]->starting_charge2 = S;
        }
    }

    //Creating the simulation Threads
    Thread_Static* simulating_threads[no_datasets];
    pthread_t simulationid[no_datasets];
    int sim_id;

    // filling up the output voltage data and starting charges for the dataset (index-2)number_of_datasets
    for(int i = 0;i<no_datasets;i++){

        relevant_data_set = (data[i])[index];
        outputvoltage[index][i] = new double[relevant_data_set->size()];

        if(index>1) {
            carryOvers[index - 2][i] = new starting_Params(0, 0, 0, 0, 0);
            sim_id = i+1;

            simulating_threads[i] = new Thread_Static(sim_id,relevant_data_set,outputvoltage[index][i],carryOvers[index][i],x[3], x[0], x[1], x[2],
                    relevant_data_set->size(), nullptr,carryOvers[index - 2][i],true);
            pthread_create (&simulationid[i],NULL,static_sim_post,(void*) simulating_threads[i]);

        }
        else {
            sim_id = i + 1;
            simulating_threads[i] = new Thread_Static(sim_id, relevant_data_set, outputvoltage[index][i], carryOvers[index][i],
                                               x[3], x[0], x[1], x[2], relevant_data_set->size());
            pthread_create(&simulationid[i], NULL, static_sim_post, (void *) simulating_threads[i]);
        }
    }

    for(auto & i : simulationid)
    {
        pthread_join (i, NULL) ;
    }

    for(int i=0;i<no_datasets;i++)
    {
        delete simulating_threads[i];
    }

}

void *static_sim_post (void *producer_parameter)
{

    //Intialise the structure passed in to thread function
    Thread_Static* thread_struct = (Thread_Static *) producer_parameter;
    vector<vector<double>> *data = thread_struct->data_set;
    double* output_voltage = thread_struct->output_voltage;
    starting_Params* starting_data = thread_struct->starting_data;
    starting_Params* ending_data = thread_struct->ending_data;
    bool final = thread_struct->final;
    double R0 = thread_struct->R0;
    double Reff = thread_struct->Reff;
    double Rct = thread_struct->Rct;
    double C = thread_struct->C;
    int length = thread_struct->length;


    static_simulate(data,output_voltage,starting_data,R0,Reff,Rct,C,ending_data,final);

    pthread_exit(0);

}

double static_objective(unsigned n, const double *x, double *grad, void *my_func_data){

    optData *d = (optData*) my_func_data;
    double Reff = x[0];
    double Rct = x[1];
    double C = x[2];
    double R0 = x[3];

    vector<vector<double>> ***relevantDataPtr = d->data;
    vector<vector<double>>* relevant_data_set = nullptr;


    double S;
    double Vc;
    if (d->index == 40 ||d->index == 39) {
        for (int i = 0; i < no_datasets; i++) {
            relevant_data_set = (relevantDataPtr[i])[d->index];

            // get the carryover struct initiated
            Vc = (*relevant_data_set)[0][2] + R0 * (*relevant_data_set)[0][1];
            S = (*relevant_data_set)[0][3];
            initial0CV(Vc,OCV_Coef,S,R0,Rct,(*relevant_data_set)[0][1]);
            d->carryOvers[i] = new starting_Params(Vc,S,S,0,0);
        }
    }


    double *output_voltage  = nullptr;
    double answer = 0;
    int length=0;

    //Creating the simulation Threads
    Thread_Static* simulating_threads[no_datasets];
    pthread_t simulationid[no_datasets];
    int sim_id;

    //to put each threads answer in;
    double answers[no_datasets];

    for (int i=0;i<no_datasets;i++){

        relevant_data_set = (relevantDataPtr[i])[d->index];
        length = relevant_data_set->size();

        sim_id = i+1;
        simulating_threads[i] = new Thread_Static(sim_id,relevant_data_set,output_voltage,d->carryOvers[i],R0,Reff,Rct,C,length,&answers[i]);
        pthread_create (&simulationid[i],NULL,static_sim,(void*) simulating_threads[i]);
    }

    for(auto & i : simulationid)
    {
        pthread_join (i, NULL) ;
    }

    for(int i=0;i<no_datasets;i++)
    {
        delete simulating_threads[i];
        answer += answers[i];
    }

    return answer;
}

void *static_sim (void *producer_parameter)
{

    //Intialise the structure passed in to thread function
    Thread_Static* thread_struct = (Thread_Static *) producer_parameter;

    vector<vector<double>> *data = thread_struct->data_set;
    double* output_voltage = thread_struct->output_voltage;
    starting_Params* starting_data = thread_struct->starting_data;
    double R0 = thread_struct->R0;
    double Reff = thread_struct->Reff;
    double Rct = thread_struct->Rct;
    double C = thread_struct->C;
    int length = thread_struct->length;

    output_voltage = new double[length];

    double* answer = thread_struct->answer;

    static_simulate(data,output_voltage,starting_data,R0,Reff,Rct,C);

    *answer = L2_norm_distance(output_voltage,data,length);

    //stops the memory from leaking for the output voltage array
    delete[] output_voltage;

    pthread_exit(0);

}

void static_simulate(vector<vector<double>> *data,double output_voltage[],starting_Params* starting_data,double R0,double Reff,
        double Rct,double C,starting_Params* ending_data,bool final){

    double starting_time = (*data).at(0).at(0);

    vector<double>* new_data = new vector<double>;
    new_data->push_back((*data)[0][0]);

    double ittime=0.0;
    for(size_t i=1;i<data->size();i++) {
        if (((*data)[i][0] - (*data)[i-1][0])  > 0.1) {
            ittime = (*data)[i - 1][0];
            //while time is smaller than the next one
            while (ittime < (*data)[i][0]) {
                ittime += 0.1;
                new_data->push_back(ittime);
            }
        }else{
            new_data->push_back((*data)[i][0]);
        }
    }

    size_t length = new_data->size();
    double dt[length-1];//,time[length];

    double* time = new double[length];
    double* I = new double[length];
    double* dt_= new double[length];

    for(size_t i=0;i<length;i++){
        time[i] = (*new_data)[i] -(*new_data)[0];
        I[i] = interpolate(data, time[i], true, 1);
        if(i>0){
            dt_[i-1] = time[i] - time[i-1];
        }
    }
    delete new_data;

    for(size_t i=1;i<length;i++){
        dt[i-1] = time[i] - time[i-1];
    }

    // now for the simulation
    double* V = new double[length];
    double* Vc = new double[length];
    double* I2_bar=new double[length];
    double* I2=new double[length];
    double* OCV1=new double[length];
    double* OCV2=new double[length];
    double* S1=new double[length];
    double* S2=new double[length];

    double Q = 5.0 * 3600;

    Vc[0] =starting_data->Vc;
    V[0] = Vc[0] - I[0] * R0;

    S1[0] = starting_data->starting_charge1;
    S2[0] = starting_data->starting_charge2;

    OCV1[0] = polyfit(S1[0],OCV_Coef);
    OCV2[0] = polyfit(S2[0],OCV_Coef);

    I2[0] = starting_data->I2;
    I2_bar[0] = starting_data->I2_bar;

    double I_[2];

    for(size_t i=1;i<length;i++){

        //calculating stuff from the previous timestep
        S2[i] = S2[i - 1] - dt[i-1] / (Q/2.0) * I2_bar[i - 1];
        S1[i] = S1[i - 1] - dt[i-1]/(Q/2.0) * I2[i - 1];

        OCV1[i] = polyfit(S1[i],OCV_Coef);
        OCV2[i] = polyfit(S2[i],OCV_Coef);

        I_[0] = I[i - 1];
        I_[1] = I[i];
        Vc[i] = implicitRK4(time[i-1], Vc[i-1], dt[i-1],I_, &OCV1[i - 1], &C, &Rct,true);

        I2_bar[i] = (OCV2[i] - OCV1[i]) / Reff;

        if(OCV1[i] == OCV2[i])
            I2[i] = I[i];
        else
            I2[i] = ((OCV1[i] - Vc[i]) / Rct) - ((OCV2[i] - OCV1[i]) / Reff);

        V[i] = Vc[i] - I[i] * R0;
    }

    cout << V[length-1] << endl;

    for(size_t i =0;i<data->size();i++)
        output_voltage[i] = interpolate(time, V, (*data)[i][0] - starting_time, length, true);

    if(final){
        ending_data->Vc = Vc[length-1];
        ending_data->starting_charge1 = S1[length-1];
        ending_data->starting_charge2 = S2[length-1];
        ending_data->I2 = I2[length-1];
        ending_data->I2_bar = I2_bar[length-1];
    }

    // now for the simulation
    delete[] V;
    delete[] Vc;
    delete[] I2_bar;
    delete[] I2;
    delete[] OCV1;
    delete[] OCV2;
    delete[] S1;
    delete[] S2;

    delete[] I;
    delete[] time;
    delete[] dt_;

}



