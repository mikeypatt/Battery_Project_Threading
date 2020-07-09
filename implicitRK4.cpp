//
// Created by Michael Patterson on 29/06/2020.
//
#include <iostream>
#include <cmath>
#include <nlopt.hpp>
#include <vector>
#include "implicitRK4.h"
#include "helpers.h"

double solveKs(unsigned n, const double *x, double *grad, void *my_func_data){

    kOptData *d = (kOptData*) my_func_data;
    double t1,y1,t2,y2;
    double f1,f2,newdt,In1,OCV1_ ,C_,Rct_;
    t1 = d->time + ((1/2.0)-(1/6.0)*sqrt(3))*d->dt;
    y1 = d->Vc + (1/4.0)*d->dt*x[0] + ((1/4.0)-(1/6.0)*sqrt(3))*d->dt*x[1];

    newdt = ((1/2.0)-(1/6.0)*sqrt(3))*d->dt;
    In1 =  (d->I)[0] + ((d->I)[1] - (d->I)[0]) * (newdt/d->dt);
    OCV1_ =  (d->OCV1)[0] + ((d->OCV1)[1] - (d->OCV1)[0]) *(newdt/d->dt);

    if(!d->static_){
        C_ =  (d->C)[0] + ((d->C)[1] - (d->C)[0]) * (newdt/d->dt);
        Rct_ =  (d->Rct)[0] + ((d->Rct)[1] - (d->Rct)[0]) *(newdt/d->dt);
    }

    t2 = d->time + ((1/2.0)+(1/6.0)*sqrt(3))*d->dt;
    y2 = d->Vc + ((1/4.0)+(1/6.0)*sqrt(3))*d->dt*x[0] + (1/4.0)*d->dt*x[1];

    f1 = 0 - (x[0] - stateDeriv_implicit(t1,y1,newdt,d->I,d->OCV1,d->C,d->Rct,d->static_,In1,OCV1_,C_,Rct_));
    f2 = 0 - (x[1] - stateDeriv_implicit(t2,y2,newdt,d->I,d->OCV1,d->C,d->Rct,d->static_,In1,OCV1_,C_,Rct_));

    return sqrt(pow(f1,2) + pow(f2,2));
}


double stateDeriv_implicit(double time, double Vc, double dt, double* I, double* OCV1, double* C, double* Rct,bool static_,double newI,double newOCV1,double newC, double newRct) {
    //double newVc = (OCV1 - Vc)/(C*Rct) - (I/C);
    double newVc;
    if (static_)
        newVc = ((dt/2) * (-1 /(Rct[0] * C[0]) - 1/(Rct[0] * C[0])))*Vc + (dt/2) * (OCV1[0]/(Rct[0] * C[0]) + newOCV1/(Rct[0] * C[0])) - (dt/2) * (I[0]/C[0] + newI/C[0]);
//        newVc = ((dt/2) * (-1 /(Rct[0] * C[0]) - 1/(Rct[0] * C[0])))*Vc + (dt/2) * (OCV1[0]/(Rct[0] * C[0]) + OCV1[1]/(Rct[0] * C[0])) - (dt/2) * (I[0]/C[0] + I[1]/C[0]);
    else
        newVc = ((dt/2) * (-1 /(Rct[0] * C[0]) - 1/(newRct * newC)))*Vc + (dt/2) * (OCV1[0]/(Rct[0] * C[0]) + newOCV1/(newRct * newC)) - (dt/2) * (I[0]/C[0] + newI/newC);
//        newVc = ((dt/2) * (-1 /(Rct[0] * C[0]) - 1/(Rct[1] * C[1])))*Vc + (dt/2) * (OCV1[0]/(Rct[0] * C[0]) + OCV1[1]/(Rct[1] * C[1])) - (dt/2) * (I[0]/C[0] + I[1]/C[1]);
    return newVc;
}


double implicitRK4(double time, double Vc, double dt, double* I, double* OCV1, double* C, double* Rct,bool static_){

    const int number_of_parameters = 2;
    double lb[number_of_parameters];
    double ub[number_of_parameters];
    lb[0] = -2; // k1 lb
    lb[1] = -1;//  k2 lb
    ub[0] = -0.01; // k1 ub
    ub[1] = 2; // k2 ub

    kOptData addData(time, Vc, dt, I, OCV1, C, Rct, static_);

    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_COBYLA, number_of_parameters);

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    nlopt_set_stopval(opt, -HUGE_VAL);
    nlopt_set_ftol_abs(opt, 1e-30);

    double k[2];
    k[1] = 0.01;
    k[0] = -0.01;

    nlopt_set_xtol_rel(opt, 1e-30);

    nlopt_set_min_objective(opt, solveKs, &addData);
    double minf;
    nlopt_optimize(opt, k, &minf);

    nlopt_destroy(opt);

    return  Vc + (1/2.0)*dt*k[0] + (1/2.0)*dt*k[1];
}


double constaint1(unsigned n, const double *x, double *grad, void *my_func_data){
    SOptData *d = (SOptData*) my_func_data;
    double OCV1 = polyfit(x[0],d->OCV_Coef);
    return -(d->Vc - (OCV1 - d->I*d->Rct - d->R0*d->I));
}

void initial0CV(double Vc,const double OCV_Coef[],double& starting_charge1,double R0,double Rct,double I){

    const int number_of_parameters = 1;
    double lb[number_of_parameters];
    double ub[number_of_parameters];
    lb[0] = starting_charge1 - 0.05*starting_charge1 ; // s1 lb
    ub[0] = starting_charge1 + 0.05*starting_charge1 ;
    SOptData addData(Vc,OCV_Coef,Rct,R0,I);

    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_COBYLA, number_of_parameters);

    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);

    nlopt_set_stopval(opt,-HUGE_VAL);
    nlopt_set_ftol_abs(opt, 1e-30);

    double s[number_of_parameters];
    s[0] = starting_charge1;

    nlopt_add_inequality_constraint(opt, constaint1,&addData, 1e-8);

    nlopt_set_xtol_rel(opt, 1e-30);

    nlopt_set_min_objective(opt,solveS,&addData);
    double minf;
    nlopt_optimize(opt, s, &minf);

    nlopt_destroy(opt);

    starting_charge1 = s[0];

}

double solveS(unsigned n, const double *x, double *grad, void *my_func_data){
    SOptData *d = (SOptData*) my_func_data;
    double OCV1 = polyfit(x[0],d->OCV_Coef);
    return  d->Vc - (OCV1 - d->I*d->Rct - d->R0*d->I);
}







