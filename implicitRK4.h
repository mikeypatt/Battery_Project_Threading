//
// Created by Michael Patterson on 29/06/2020.
//
#ifndef UNTITLED1_IMPLICITRK4_H
#define UNTITLED1_IMPLICITRK4_H

double implicitRK4(double time, double Vc, double dt, double* I, double* OCV1, double* C, double* Rct,bool static_);
double solveKs(double time, double Vc, double dt, double I, double OCV1, double C, double Rct);
double stateDeriv_implicit(double time, double Vc, double dt, double* I, double* OCV1, double* C, double* Rct,bool static_,double newI,double newOCV1,double newC, double newRct);

struct kOptData{
    double time;
    double Vc;
    double dt;
    double* I;
    double* OCV1;
    double* C;
    double* Rct;
    bool static_;
    kOptData(double time,double Vc,double dt, double* I,double* OCV1,double* C,double* Rct,bool static_):
            time(time),Vc(Vc),dt(dt),I(I),OCV1(OCV1),C(C),Rct(Rct),static_(static_){}
};


struct SOptData{
    double Vc;
    const double* OCV_Coef;
    double Rct;
    double R0;
    double I;
    SOptData(double Vc,const double OCV_Coef[],double Rct,double R0,double I):Vc(Vc),OCV_Coef(OCV_Coef),Rct(Rct),R0(R0),I(I){}
};
double constaint1(unsigned n, const double *x, double *grad, void *my_func_data);
void initial0CV(double Vc, const double *OCV_Coef, double& starting_charge1, double R0, double Rct, double I);
double solveS(unsigned n, const double *x, double *grad, void *my_func_data);
#endif //UNTITLED1_IMPLICITRK4_H
