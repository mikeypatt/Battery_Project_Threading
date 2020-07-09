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




//void cumTrapz(double batteryCap,vector<vector<double>>& data){
//
//    int length = data.size();
//    double dt,max=-9999;
//    int width = data[0].size();
//
//    for(int i=1;i<length-1;i++){
//        dt = data[i][0] - data[i-1][0];
//        if(i>0) {
//            data[i][width - 1] = (data[i - 1][width - 1] + (((-1*data[i][1]/batteryCap) + (-1*data[i + 1][1]/batteryCap)) / 2) * dt);
//        }else {
//            data[i][width - 1] = (-1 * ((data[i][1]/batteryCap + data[i + 1][1]/batteryCap) / 2) * dt);
//        }
//        if(data[i][width-1]>max)
//            max = data[i][width-1];
//    }
//
//    for(int i=0;i<length;i++){
//        data[i][width-1] = data[i][width-1] + (1-max);
//    }
//}