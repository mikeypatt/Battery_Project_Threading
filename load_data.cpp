#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>

#include "load_data.h"

using namespace std;

void readInData(const char* filename, vector<vector<double>>& data)
{
    /* function has the file name and a 2D vector as input and fills up the vector
     * with the files inputs of time current and voltage
     */
    int i;
    double time,current,voltage,currentNumber,charge;
    char line[300];

    //opens the file
    ifstream input_stream;
    input_stream.open(filename);
    if(input_stream.fail()){
        cout << "could not open file" << endl;
        return;
    }

//    vector<vector<double>> **datapointers = new vector<vector<double>>*[4];

    //gets a line of the file
    while(!input_stream.eof()){
        input_stream.getline(line,300);
        stringstream stream(line);
        i = 0;

        //iterates through the entire line
        while(!stream.eof()){
            stream >> currentNumber;
            if(i==0)
                time = currentNumber;
            if(i==1)
                current = currentNumber;
            if(i==2)
                voltage = currentNumber;
            if(i==3)
                charge = currentNumber;
            i++;
        }

        vector<double> row{time,current,voltage,charge};
        data.push_back(row);
    }


    //closes the input stream
    input_stream.close();

}


void cumTrapz(double batteryCap,vector<vector<double>>& data){

    int length = data.size();
    double dt,max=-9999;
    int width = data[0].size();

    for(int i=1;i<length-1;i++){
        dt = data[i][0] - data[i-1][0];
        if(i>0) {
            data[i][width - 1] = (data[i - 1][width - 1] + (((-1*data[i][1]/batteryCap) + (-1*data[i + 1][1]/batteryCap)) / 2) * dt);
        }else {
            data[i][width - 1] = (-1 * ((data[i][1]/batteryCap + data[i + 1][1]/batteryCap) / 2) * dt);
        }
        if(data[i][width-1]>max)
            max = data[i][width-1];
    }

    for(int i=0;i<length;i++){
        data[i][width-1] = data[i][width-1] + (1-max);
    }
}

vector<vector<double>> **dataSplitter(vector<vector<double>> data,double split_points[],int charge_points_length){

    double minimum = 0,maximum =0;

    vector<vector<double>> **datapointers = new vector<vector<double>>*[charge_points_length];

    for(int i= 0;i<charge_points_length;i++){
        if(i==0){
            minimum = split_points[i];
        }
        else{
            minimum = split_points[i-1];
        }
        if(i ==(charge_points_length-1)){
            maximum = split_points[i];
        }
        else{
            maximum = split_points[i+1];
        }

        vector<vector<double>> *vect = new vector<vector<double>>;

        for(auto & j : data){
            if (j[3] > minimum && j[3] <= maximum){
                vector<double> row{j[0],j[1],j[2],j[3]};
                vect->push_back(row);
            }
        }
        datapointers[i] = vect;
    }

    return datapointers;
}

double* splitting_points_opt(){

    int length = 11;
    double * splitting_points = new double[length];
    splitting_points[0] = 0.15;
    splitting_points[1] = 0.25;
    splitting_points[2] = 0.3;
    splitting_points[3] =0.34;
    splitting_points[4] =0.4;
    splitting_points[5] =0.59;
    splitting_points[6] =0.66;
    splitting_points[7] = 0.80;
    splitting_points[8] =0.87;
    splitting_points[9] =0.94;
    splitting_points[10]=1.0;
    return splitting_points;
}

double* splitting_points_ref(){

    int length = 41;
    double* split_points = new double[length];

    // segments of the break points
    split_points[0] = 0.15;
    split_points[5] = 0.225;
    split_points[33] = 0.9125;

    for(int i = 0;i<=length;i++){
        if(0<i && i<5){
            split_points[i] =  split_points[i-1] + 0.0125;
        }
        if(i>5 && i < 33){
            split_points[i] =  split_points[i-1] + 0.025;
        }
        if(i>33 && i < 41){
            split_points[i] =  split_points[i-1] + 0.0125;
        }
    }
    return split_points;
}

double* splitting_points_new(){

    int length = 86;
    double* split_points = new double[length];

    // segments of the break points
    split_points[0] = 0.15;
    for(int i = 1;i<length;i++){
        split_points[i] =  split_points[i-1] + 0.01;
    }
    return split_points;
}

