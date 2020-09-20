//
// Created by Michael Patterson on 30/06/2020.
//
#include <vector>
#include <nlopt.hpp>
#include <random>
#include <pthread.h>
#include <cctype>
#include <algorithm>


#include "optimisations.h"
#include "implicitRK4.h"
#include "helpers.h"

using namespace std;


void static_fit(vector<vector<double>> ***data,double carry_over_parameters[],
                starting_Params* carryOvers[][no_datasets],double* outputvoltage[][no_datasets],int index,
                double Reff_table[],double Rct_table[],double C_table[],double R0_table[]){

    const int number_of_parameters = 4 + no_datasets*3;

    double lb[number_of_parameters];
    double ub[number_of_parameters];
    double x [number_of_parameters];

    if(index==40 || index==39){
        x[0] = carry_over_parameters[0];
        x[1] = carry_over_parameters[1];
        x[2] = 0.347126;//carry_over_parameters[2];
        x[3] =  (polyfit((*data[0][index])[0][3], R0_Coef));

        for(int i =0;i<no_datasets;i++){
            //starting Vcs
            x[4+i] = (*data[i][index])[0][2] + x[3] * (*data[i][index])[0][1];
            //starting S1
            x[4+no_datasets+i] = (*data[i][index])[0][3];
            initial0CV(x[4+i],OCV_Coef,x[4+no_datasets+i],x[3],x[1],(*data[i][index])[0][1]);
            //setting S2 = S1
            x[4+(2*no_datasets)+i] = x[4+no_datasets+i];
        }
    }else{
        x[0] = Reff_table[index+2];
        x[1] = Rct_table[index+2];
        x[2] = 0.347126;
        x[3] = R0_table[index+2];
        for(int i =0;i<no_datasets;i++){
            //starting Vcs
            x[4+i] = carryOvers[index][i]->Vc;
            //starting S1
            x[4+no_datasets+i] =carryOvers[index][i]->starting_charge1;
            //setting S2
            x[4+(2*no_datasets)+i] = carryOvers[index][i]->starting_charge2;
        }
    }

    lb[0] = 0.003;
    lb[1] = 0.001;
    ub[0] = 0.7;
    ub[1] = 0.1;

    lb[2] = 0.001;
    ub[2]  = 1;

    lb[3] = 0.001;
    ub[3]  = 2;

    for(int i =0;i<no_datasets;i++)
    {
        if(index==40 || index==39) {
            ub[4+i] = x[4+i];
            lb[4+i] = x[4+i];
            ub[4+no_datasets+i] = x[4+no_datasets+i];
            lb[4+no_datasets+i] =  x[4+no_datasets+i];
            ub[4+(2*no_datasets)+i] = x[4+(2*no_datasets)+i];
            lb[4+(2*no_datasets)+i] = x[4+(2*no_datasets)+i];
        }else{
            ub[4+i] = x[4+i] + 0.03*x[4+i];
            lb[4+i] = x[4+i] - 0.03*x[4+i];
            ub[4+no_datasets+i] = x[4+no_datasets+i] + 0.03*x[4+no_datasets+i];
            lb[4+no_datasets+i] = x[4+no_datasets+i] - 0.03*x[4+no_datasets+i];
            ub[4 + (2 * no_datasets) + i] = x[4 + (2 * no_datasets) + i] + 0.03*x[4 + (2 * no_datasets) + i];
            lb[4 + (2 * no_datasets) + i] = x[4 + (2 * no_datasets) + i] - 0.03*x[4 + (2 * no_datasets) + i];
        }
    }

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

    nlopt_add_equality_constraint(opt,startingV_con_1, &addData, 1e-30);
    nlopt_add_equality_constraint(opt,OCV_Con, &addData, 1e-30);
    //nlopt_add_inequality_constraint(opt,charges_con, &addData, 1e-8);
    nlopt_set_min_objective(opt,static_objective,&addData);

    nlopt_set_maxeval(opt,1);
    double minf;
    int res = nlopt_optimize(opt, x, &minf);

    if (res <0) {
        printf("nlopt failed!\n");
    }
    else{
        printf("index %d found minimum at f(%g,%g,%g,%g,%g,%g,%g) = %0.10g\n",index, x[0],x[1],x[2], x[3],x[4],x[5],x[6],minf);
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



    vector<vector<double>>* relevant_data_set = nullptr;

    for (int i = 0; i<no_datasets;i++) {
        carryOvers[index][i]->Vc = x[4+i];
        carryOvers[index][i]->starting_charge1 = x[4+no_datasets+i];
        carryOvers[index][i]->starting_charge2 = x[4+(2*no_datasets)+i];
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

    if (d->index == 40 ||d->index == 39) {
        for (int i = 0; i < no_datasets; i++) {
            d->carryOvers[i] = new starting_Params(x[4+i],x[4+no_datasets+i],x[4+(2*no_datasets)+i],0,0);
        }
    }else{
        for (int i = 0; i < no_datasets; i++) {
            d->carryOvers[i]->Vc = x[4+i];
            d->carryOvers[i]->starting_charge1 = x[4+no_datasets+i];
            d->carryOvers[i]->starting_charge2 = x[4+(2*no_datasets)+i];
        }
    }

    double *output_voltage  = nullptr;
    double answer = 0;
    int length=0;

    struct timespec start, end, elapsed;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);



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

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    elapsed.tv_sec = end.tv_sec - start.tv_sec;
    elapsed.tv_nsec = end.tv_nsec - start.tv_nsec;
    if (elapsed.tv_nsec < 0) {
        elapsed.tv_nsec += 1e9;
        elapsed.tv_sec -= 1;
    }
    printf("Elapsed time = %f seconds.\n", (double)elapsed.tv_sec + (double)elapsed.tv_nsec / 1.0e9);

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

    *answer = Trapz(data,output_voltage);

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
        if (((*data)[i][0] - (*data)[i-1][0]+0.1)  > 0.1) {
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

    double* dt = new double[length-1];
    double* time = new double[length];
    double* I = new double[length];

    for(size_t i=0;i<length;i++){
        time[i] = (*new_data)[i] -(*new_data)[0];
        I[i] = interpolate(data, time[i], true, 1);
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

    //calculating stuff from the previous timestep
    S2[0] = starting_data->starting_charge2;
    S1[0] = starting_data->starting_charge1;
    Vc[0] = starting_data->Vc;


    OCV1[0] = polyfit(S1[0],OCV_Coef);
    OCV2[0] = polyfit(S2[0],OCV_Coef);

    I2[0] = ((OCV1[0] - Vc[0]) / Rct) - ((OCV2[0] - OCV1[0]) / Reff);
    I2_bar[0] = (OCV2[0] - OCV1[0]) / Reff;

    V[0] = Vc[0] - I[0] * R0;


    double newI,newOCV,Fn;

    for(size_t i=1;i<length;i++){

        //calculating stuff from the previous timestep
        S2[i] = S2[i - 1] - dt[i-1] / (Q/2.0) * I2_bar[i - 1];
        S1[i] = S1[i - 1] - dt[i-1] / (Q/2.0) * I2[i - 1];

        OCV1[i] = polyfit(S1[i],OCV_Coef);
        OCV2[i] = polyfit(S2[i],OCV_Coef);

        newI =  (I[i-1] + I[i])/2.0;
        newOCV =  (OCV1[i-1] + OCV1[i])/2.0;
        Fn = (newOCV/(C*Rct) - newI/C);
        Vc[i] = Vc[i-1] * exp(-dt[i-1]/(C*Rct)) + C*Rct * Fn * (1 - exp(-dt[i-1]/(C*Rct)));

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

    delete[] I;
    delete[] time;
    delete[] dt;

}


void dynamic_fit(vector<vector<double>> ***data,starting_Params* carryOvers[][no_datasets],
                 double* outputvoltage[][no_datasets],int index,double Reff_table[],double Rct_table[],
                 double C_table[],double R0_table[],double spoints_ref[]){

    const int number_of_parameters = 4 + no_datasets*3;

    double lb[number_of_parameters];
    double ub[number_of_parameters];
    double x [number_of_parameters];

    x[0] = Reff_table[index];
    x[1] = Rct_table[index];
    x[2] = C_table[index];
    x[3] = R0_table[index];

    if(index==40 || index==39){
        for(int i =0;i<no_datasets;i++){
            //starting Vcs
            x[4+i] = (*data[i][index])[0][2] + x[3] * (*data[i][index])[0][1];
            //starting S1
            x[4+no_datasets+i] = (*data[i][index])[0][3];
            initial0CV(x[3+i],OCV_Coef,x[4+no_datasets+i],x[3],x[1],(*data[i][index])[0][1]);
            //setting S2 = S1
            x[4+(2*no_datasets)+i] = x[4+no_datasets+i];
        }
    }else{
        for(int i =0;i<no_datasets;i++){
            //starting Vcs
            x[4+i] = carryOvers[index][i]->Vc;
            //starting S1
            x[4+no_datasets+i] =carryOvers[index][i]->starting_charge1;
            //setting S2
            x[4+(2*no_datasets)+i] = carryOvers[index][i]->starting_charge2;
        }
    }

    lb[0] = 0.003;
    lb[1] = 0.001;
    ub[0] = 0.7;
    ub[1] = 0.1;

    lb[2] = x[2] - 0.03*x[2];
    ub[2]  = x[2] + 0.03*x[2];

    lb[3] = x[3] - 0.03*x[3];
    ub[3]  = x[3] + 0.03*x[3];

    for(int i =0;i<no_datasets;i++)
    {
        if(index==40 || index==39) {
            ub[4+i] = x[4+i];
            lb[4+i] = x[4+i];
            ub[4+no_datasets+i] = x[4+no_datasets+i];
            lb[4+no_datasets+i] =  x[4+no_datasets+i];
            ub[4+(2*no_datasets)+i] = x[4+(2*no_datasets)+i];
            lb[4+(2*no_datasets)+i] = x[4+(2*no_datasets)+i];
        }else{
            ub[4+i] = x[4+i] + 0.03*x[4+i];
            lb[4+i] = x[4+i] - 0.03*x[4+i];
            ub[4+no_datasets+i] = x[4+no_datasets+i] + 0.03*x[4+no_datasets+i];
            lb[4+no_datasets+i] = x[4+no_datasets+i] - 0.03*x[4+no_datasets+i];
            ub[4 + (2 * no_datasets) + i] = x[4 + (2 * no_datasets) + i] + 0.03*x[4 + (2 * no_datasets) + i];
            lb[4 + (2 * no_datasets) + i] = x[4 + (2 * no_datasets) + i] - 0.03*x[4 + (2 * no_datasets) + i];
        }
    }

    //create addtional datastruct
    optData addData(data,index,carryOvers[index],Reff_table,Rct_table,C_table,R0_table,spoints_ref);

    //creates optimsation object
    nlopt_opt opt; // NLOPT_LD_SLSQP
    opt = nlopt_create(NLOPT_LN_COBYLA , number_of_parameters);

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    nlopt_set_stopval(opt,-HUGE_VAL);

    nlopt_set_ftol_abs(opt, 1e-30);
    nlopt_set_xtol_rel(opt, 1e-30);

    nlopt_add_equality_constraint(opt,startingV_con_1, &addData, 1e-30);
    nlopt_add_equality_constraint(opt,OCV_Con, &addData, 1e-30);
    //nlopt_add_inequality_constraint(opt,charges_con, &addData, 1e-8);
    nlopt_set_min_objective(opt,dynamic_objective,&addData);

    nlopt_set_maxeval(opt,10);
    double minf;
    int res = nlopt_optimize(opt, x, &minf);

    if (res <0) {
        printf("nlopt failed!\n");
    }
    else{
        printf("index %d found minimum at f(%g,%g,%g,%g,%g,%g,%g) = %0.10g\n",index, x[0],x[1],x[2], x[3],x[4],x[5],x[6],minf);
    }

    nlopt_destroy(opt);

    //filling up lookup tables with values
    Reff_table[index] = x[0];
    Rct_table[index] = x[1];
    C_table[index] = x[2];
    R0_table[index] = x[3];

    vector<vector<double>>* relevant_data_set = nullptr;

    for (int i = 0; i<no_datasets;i++) {
        carryOvers[index][i]->Vc = x[4+i];
        carryOvers[index][i]->starting_charge1 = x[4+no_datasets+i];
        carryOvers[index][i]->starting_charge2 = x[4+(2*no_datasets)+i];
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
            pthread_create (&simulationid[i],NULL,dynamic_sim_post,(void*) simulating_threads[i]);

        }
        else {
            sim_id = i + 1;
            simulating_threads[i] = new Thread_Dynamic(sim_id,relevant_data_set,outputvoltage[index][i],carryOvers[index][i],x[3], x[0], x[1], x[2],
                                                       relevant_data_set->size(),Reff_table,Rct_table,C_table,R0_table,spoints_ref,index);
            pthread_create(&simulationid[i], NULL, dynamic_sim_post, (void *) simulating_threads[i]);
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

    for (int i = 0; i < no_datasets; i++) {
        d->carryOvers[i]->Vc = x[4+i];
        d->carryOvers[i]->starting_charge1 = x[4+no_datasets+i];
        d->carryOvers[i]->starting_charge2 = x[4+(2*no_datasets)+i];
    }


    double *output_voltage  = nullptr;
    double answer = 0;
    size_t length=0;

    clock_t start, end;
    start = clock();

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

    *answer = Trapz(data,output_voltage);

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
        if (((*data)[i][0] - (*data)[i-1][0]+0.1)  > 0.1) {
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

    I2[0] = ((OCV1[0] - Vc[0]) / Rct_[0]) - ((OCV2[0] - OCV1[0]) / Reff_[0]);
    I2_bar[0] = (OCV2[0] - OCV1[0]) / Reff_[0];

    double newI,newOCV,Fn,newC,newRct;
    for(size_t i=1;i<length;i++){

        //calculating stuff from the previous timestep
        S2[i] = S2[i - 1] - dt[i-1] / (Q/2.0) * I2_bar[i - 1];
        S1[i] = S1[i - 1] - dt[i-1]/(Q/2.0) * I2[i - 1];

        OCV1[i] = polyfit(S1[i],OCV_Coef);
        OCV2[i] = polyfit(S2[i],OCV_Coef);

        newI =  (I[i - 1] + I[i])/2.0;
        newOCV =  (OCV1[i-1] + OCV1[i])/2.0;
        newC =  (C_[i] + C_[i-1])/2.0;
        newRct =  (Rct_[i] + Rct_[i-1])/2.0;
        Fn = (newOCV/(newC*newRct) - newI/newC);
        Vc[i] = Vc[i-1] * exp(-dt[i-1]/(newC*newRct)) + newC*newRct * Fn * (1 - exp(-dt[i-1]/(newC*newRct)));

        I2_bar[i] = (OCV2[i] - OCV1[i]) / Reff_[i];

        if(OCV1[i] == OCV2[i])
            I2[i] = I[i];
        else
            I2[i] = ((OCV1[i] - Vc[i]) / Rct_[i]) - ((OCV2[i] - OCV1[i]) / Reff_[i]);

        V[i] = Vc[i] - I[i] * R0_[i];
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

double static_sim_mini(vector<vector<double>> *data,starting_Params* starting_data,double R0,double Reff,double Rct,double C)
{

    double Q = 5.0 * 3600;

    //now for the simulation
    double V[2];
    double Vc[2];
    double I2_bar[2];
    double I2[2];
    double OCV1[2];
    double OCV2[2];
    double S1[2];
    double S2[2];

    //calculating stuff from the previous timestep
    S2[0] = starting_data->starting_charge2;
    S1[0] = starting_data->starting_charge1;
    Vc[0] = starting_data->Vc;

    double I_[2];

    OCV1[0] = polyfit(S1[0],OCV_Coef);
    OCV2[0] = polyfit(S2[0],OCV_Coef);

    I2[0] = ((OCV1[0] - Vc[0]) / Rct) - ((OCV2[0] - OCV1[0]) / Reff);
    I2_bar[0] = (OCV2[0] - OCV1[0]) / Reff;

    V[0] = Vc[0] - (*data)[0][1] * R0;
    double dt = (*data)[1][0] - (*data)[0][0];

    double newI,newOCV,Fn;
    for(size_t i=1;i<2;i++){

        //calculating stuff from the previous timestep
        S2[i] = S2[i - 1] - dt / (Q/2.0) * I2_bar[i - 1];
        S1[i] = S1[i - 1] - dt/ (Q/2.0) * I2[i - 1];

        OCV1[i] = polyfit(S1[i],OCV_Coef);
        OCV2[i] = polyfit(S2[i],OCV_Coef);

        newI =  ((*data)[i-1][1] + (*data)[i][1])/2.0;
        newOCV =  (OCV1[i-1] + OCV1[i])/2.0;
        Fn = (newOCV/(C*Rct) - newI/C);
        Vc[i] = Vc[i-1] * exp(-dt/(C*Rct)) + C*Rct * Fn * (1 - exp(-dt/(C*Rct)));

        I2_bar[i] = (OCV2[i] - OCV1[i]) / Reff;

        if(OCV1[i] == OCV2[i])
            I2[i] = (*data)[i][1];
        else
            I2[i] = ((OCV1[i] - Vc[i]) / Rct) - ((OCV2[i] - OCV1[i]) / Reff);

        V[i] = Vc[i] - (*data)[i][1] * R0;
    }

    return V[1];

}

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

