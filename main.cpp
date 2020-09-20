#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>

#include "load_data.h"
#include "optimisations.h"

using namespace std;


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
    fitting_data[0] = CC_02_CC_data;
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

    // creates an struct for the starting values
    starting_Params* starting_params[41][no_datasets];

    double carry_over_parameters[3] = {0.123611,0.0136325,0.32427};
    double* outputVoltage[charge_point_length][no_datasets];

    //reverse iterates from a high soc to a lower one
    for(int i=40;i>=0;i--){
        static_fit(fitting_data,carry_over_parameters,starting_params,outputVoltage,
                       i,Reff_table, Rct_table, C_table, R0_table);
    }

    // delete all the heap dataaaaaa
    ofstream out_stream;

    out_stream.open("semistatic_tables.txt");
    for(int i=0;i<charge_point_length;i++){
        out_stream << R0_table[i] << " " << Reff_table[i] << " " << Rct_table[i] << " " << C_table[i] << endl;
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_02C_CC.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[0])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[0])[ii])[j][2] << " " << outputVoltage[ii][0][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_02C_long.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[1])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[1])[ii])[j][2] << " " << outputVoltage[ii][1][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_02C_short.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[2])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[2])[ii])[j][2] << " " << outputVoltage[ii][2][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_03C_CC.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[3])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[3])[ii])[j][2] << " " << outputVoltage[ii][3][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_03C_long.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[4])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[4])[ii])[j][2] << " " << outputVoltage[ii][4][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_03C_short.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[5])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[5])[ii])[j][2] << " " << outputVoltage[ii][5][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_04C_CC.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[6])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[6])[ii])[j][2] << " " << outputVoltage[ii][6][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_04C_long.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[7])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[7])[ii])[j][2] << " " << outputVoltage[ii][7][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_04C_short.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[8])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[8])[ii])[j][2] << " " << outputVoltage[ii][8][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_05C_CC.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[9])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[9])[ii])[j][2] << " " << outputVoltage[ii][9][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_05C_long.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[10])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[10])[ii])[j][2] << " " << outputVoltage[ii][10][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_05C_short.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[11])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[11])[ii])[j][2] << " " << outputVoltage[ii][11][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_1C_pulse.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[12])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[12])[ii])[j][2] << " " << outputVoltage[ii][12][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_drive_long.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[13])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[13])[ii])[j][2] << " " << outputVoltage[ii][13][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_drive_short.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[14])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[14])[ii])[j][2] << " " << outputVoltage[ii][14][j] << endl;
        }
    }
    out_stream.close();

    out_stream.open("semistatic_voltage_pulses.txt");
    for(int ii=0;ii<charge_point_length;ii++){
        int length = (*(fitting_data[15])[ii]).size();
        for(int j=0;j<length;j++){
            out_stream << (*(fitting_data[15])[ii])[j][2] << " " << outputVoltage[ii][15][j] << endl;
        }
    }
    out_stream.close();

    //delete output voltages
    for(int j =0;j<charge_point_length;j++){
        for(int ii=0;ii<no_datasets;ii++){
            delete[] outputVoltage[j][ii];
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

        for(int j=0;j<charge_point_length;j++){
            dynamic_index = order[j];
            dynamic_fit(fitting_data,starting_params,outputVoltagedynamic,dynamic_index,
                    Reff_table, Rct_table, C_table, R0_table,spoints_ref);
        }

        out_stream.open("dynamic_tables.txt");
        for(int j=0;j<charge_point_length;j++){
            out_stream << Reff_table[j] << " " << Rct_table[j] << " " << C_table[j] << endl;
        }
        out_stream.close();


        out_stream.open("dynamic_voltage_02C_CC.txt");
        for(int ii=0;i<charge_point_length;i++){
            int length = (*(fitting_data[0])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[0])[ii])[j][2] << " " << outputVoltagedynamic[ii][0][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_02C_long.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[1])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[1])[ii])[j][2] << " " << outputVoltagedynamic[ii][1][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_02C_short.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[2])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[2])[ii])[j][2] << " " << outputVoltagedynamic[ii][2][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_03C_CC.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[3])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[3])[ii])[j][2] << " " << outputVoltagedynamic[ii][3][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_03C_long.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[4])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[4])[ii])[j][2] << " " << outputVoltagedynamic[ii][4][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_03C_short.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[5])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[5])[ii])[j][2] << " " << outputVoltagedynamic[ii][5][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_04C_CC.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[6])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[6])[ii])[j][2] << " " << outputVoltagedynamic[ii][6][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_04C_long.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[7])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[7])[ii])[j][2] << " " << outputVoltagedynamic[ii][7][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_04C_short.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[8])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[8])[ii])[j][2] << " " << outputVoltagedynamic[ii][8][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_05C_CC.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[9])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[9])[ii])[j][2] << " " << outputVoltagedynamic[ii][9][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_05C_long.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[10])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[10])[ii])[j][2] << " " << outputVoltagedynamic[ii][10][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_05C_short.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[11])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[11])[ii])[j][2] << " " << outputVoltagedynamic[ii][11][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_1C_pulse.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[12])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[12])[ii])[j][2] << " " << outputVoltagedynamic[ii][12][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_drive_long.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[13])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[13])[ii])[j][2] << " " << outputVoltagedynamic[ii][13][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_drive_short.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[14])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[14])[ii])[j][2] << " " << outputVoltagedynamic[ii][14][j] << endl;
            }
        }
        out_stream.close();

        out_stream.open("dynamic_voltage_pulses.txt");
        for(int ii=0;ii<charge_point_length;ii++){
            int length = (*(fitting_data[15])[ii]).size();
            for(int j=0;j<length;j++){
                out_stream << (*(fitting_data[15])[ii])[j][2] << " " << outputVoltagedynamic[ii][15][j] << endl;
            }
        }
        out_stream.close();

        //delete output voltages
        for(int j =0;j<charge_point_length;j++){
            for(int ii=0;ii<no_datasets;ii++){
                delete[] outputVoltagedynamic[j][ii];
            }
        }
    }

    delete spoints_ref;

    return 0;
}





