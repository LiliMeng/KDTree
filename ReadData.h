/**
 * File: ReadData.h
 * Author: Lili Meng (lilimeng1103@gmail.com)
 * read the line and column data from .txt file
 */

#ifndef READDATA_H
#define READDATA_H

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace std;

class readData{
public:

    readData(string filename);
    vector <vector<double>> allDataPointsVec;
    int get_num_of_elements();
    int get_num_of_dimensions();
    void printDataID();
};

readData::readData( string filename )
{
    std::ifstream fin(filename.c_str(),std::ios::in);
    if (!fin.is_open())
    {
        cout<<"cannot open the file"<<endl;
    }

    istringstream istr;
    string str;
    vector <double> dataPointVec;
    double oneDimension;
    while(getline(fin,str))
    {
        istr.str(str);
        while(istr>>oneDimension)
        {
            dataPointVec.push_back(oneDimension);
        }
        allDataPointsVec.push_back(dataPointVec);
        dataPointVec.clear();
        istr.clear();
        str.clear();
    }
    fin.close();

}

int readData::get_num_of_elements(){

    return allDataPointsVec.size();
}

int readData::get_num_of_dimensions(){

return allDataPointsVec[0].size();
}


void readData::printDataID(){

   int numOfDimensions=allDataPointsVec[0].size();
   int numOfElements=allDataPointsVec.size();


   double dataPoints[numOfElements][numOfDimensions];
   for(int i=0; i<numOfElements; i++)
   {
        for(int j=0; j<numOfDimensions; j++)
        {
             dataPoints[i][j]=allDataPointsVec[i][j];
             cout<<"A"<<i<<j<<endl;
        }
   }

 }


#endif // READDATA_H
