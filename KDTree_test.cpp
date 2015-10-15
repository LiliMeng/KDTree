#include <iostream>
#include <assert.h>
#include <unordered_map>
#include "KDTree.h"
#include "ReadData.h"
#include "Point.h"

using namespace std;


//A brute force method to check the distances
double distance_sq(const vector<double> & data0, const vector<double> & data1)
{
    assert(data0.size() == data1.size());
    double sq = 0.0;
    for (int i=0; i<data0.size(); i++)
    {
        double dif = data0[i]-data1[i];
        sq +=dif*dif;
    }
    return sq;
}



int main(int argc, const char * argv[])
{

    int K = 3;
    vector<vector<double> > dataset;
    ReadData rd1("sample_data.txt");
    dataset=rd1.allDataPointsVec;
    int N=dataset.size();
    //query_point
    vector<double> query_point;
    vector<vector<double> > query_point_dataset;
    ReadData rd2("query_points.txt");
    int N2 = rd2.get_num_of_elements();
    int dim2 = rd2.get_num_of_dimensions();
    query_point_dataset=rd2.allDataPointsVec;
    query_point=query_point_dataset[1];

    KDTree<128, size_t> kd;

    Point<128> key;
    for(int i=0; i<128; i++)
    {
        key[i]=query_point[i];
    }


    vector<size_t> pointIndices;
    for(int i=0; i<N; i++)
    {
        pointIndices.push_back(i);
    }

    kd.buildTree(dataset, pointIndices);



    vector<size_t> indices = kd.kNNValues(key, 3);


    for (int i = 0; i<indices.size(); i++)
    {
        cout<<"Using KD Tree Search K Nearest Neigbour : The number "<<i+1<<" nearest neighbor index is  "<<indices[i]<<endl;
    }



    // brute force
    unordered_map<double, int> brute_force_htable;
    vector<double> brute_force_vec;

    for (int i = 0; i< N; i++)
    {
        double dist = distance_sq(query_point, dataset[i]);

        brute_force_htable.insert({dist, i});
        brute_force_vec.push_back(dist);
    }


    std::sort(brute_force_vec.begin(), brute_force_vec.end());


    for(int i = 0; i<K; i++)
    {
        cout<<"Using Brute-Force method Search: The number "<<i+1<<" nearest neighbor index is  "<< brute_force_htable[brute_force_vec[i]]<<"\t"<<"The brute-force distance is "<<brute_force_vec[i]<<endl;
    }

    /**Compare the KD-Tree with the Brute-Force Method**/
    for (int i = 0; i<indices.size(); i++)
    {
        if(indices[i]==brute_force_htable[brute_force_vec[i]])
        {
            cout<<"Comparing with the Brute-force method, the Exact K-Nearest Neighbour search by KD-Tree program is correct"<<endl;
        }
    }



    return 0;
}

