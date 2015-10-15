/*************************************************
 * File: KDTree_test.cpp
 * Author: Lili Meng (lilimeng1103@gmail.com)
 * File containing several test cases that can be
 * used to verify the correctness of the KDTree
 * implementation.
 */

#include <iostream>
#include <assert.h>
#include <unordered_map>
#include "KDTree.h"
#include "ReadData.h"
#include "Point.h"

using namespace std;

/* A utility function to construct a Point from a range of iterators. */
template <size_t N, typename IteratorType>
Point<N> PointFromRange(IteratorType begin, IteratorType end) {
    Point<N> result;
    copy(begin, end, result.begin());
    return result;
}

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
    int N = rd1.get_num_of_elements();
    int dim = rd1.get_num_of_dimensions();
    dataset=rd1.allDataPointsVec;

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

  double dataPoints[N][dim];
   for(int i=0; i<N; i++)
   {
        for(int j=0; j<dim; j++)
        {
             dataPoints[i][j]=dataset[i][j];
        }
   }


    for (size_t i = 0; i < N; ++i)
    {
        kd.insert(PointFromRange<128>(dataPoints[i], dataPoints[i] + 128), i);
    }

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

    /** Compare the KD-Tree with the Brute-Force Method**/
    for (int i = 0; i<indices.size(); i++)
    {
        if(indices[i]==brute_force_htable[brute_force_vec[i]])
        {
            cout<<"Comparing with the Brute-force method, the Exact K-Nearest Neighbour search by KD-Tree program is correct"<<endl;
        }
    }



    return 0;
}

