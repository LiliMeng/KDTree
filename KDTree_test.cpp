/*************************************************
 * File: KDTree_test.cpp
 * Author: Lili Meng (lilimeng1103@gmail.com)
 * Based on  test-harness.cpp by Keith Schwarz (htiek@cs.stanford.edu)
 * File containing several test cases that can be
 * used to verify the correctness of the KDTree
 * implementation.
 */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cstdarg>
#include <set>
#include <fstream>

#include "KDTree.h"
#include "ReadData.h"

using namespace std;

#define BasicKDTreeTestEnabled          1 // Step one checks
#define NearestNeighborTestEnabled      1 // Step two checks

/* A utility function to construct a Point from a range of iterators. */
template <size_t N, typename IteratorType>
Point<N> PointFromRange(IteratorType begin, IteratorType end) {
    Point<N> result;
    copy(begin, end, result.begin());
    return result;
}


/* Utility function that pauses until the user hits ENTER. */
void PressEnterToContinue() {
  /* Use getline to stall until receiving input. */
  string line;
  getline(cin, line);
}

/* This function is what the test suite uses to ensure that the KDTree works
 * correctly.  It takes as parameters an expression and description, along
 * with a file and line number, then checks whether the condition is true.
 * If so, it prints that the test passed.  Otherwise, it reports that the
 * test fails and points the caller to the proper file and line.
 */
void DoCheckCondition(bool expr, const string& rationale, const string& file, int line) {
  /* It worked!  Congrats. */
  if (expr) {
    cout << "PASS: " << rationale << endl;
    return;
  }

  /* Uh oh!  The test failed! */
  cout << "FAIL: " << rationale << endl;
  cout << "  Error at " << file << ", line " << line << endl;
  cout << "  (ENTER to continue)" << endl;

  /* Pause so that the test fail stands out. */
  PressEnterToContinue();
}

/* Reports that an unexpected error occurred that caused a test to fail. */
void FailTest(const exception& e) {
  cerr << "TEST FAILED: Unexpected exception: " << e.what() << endl;
  PressEnterToContinue();
}

/* This macro takes in an expression and a string, then invokes
 * DoCheckCondition passing in the arguments along with the file
 * and line number on which the macro was called.  This makes it
 * easier to track down the source of bugs if a test case should
 * fail.
 */
#define CheckCondition(expr, rationale) DoCheckCondition(expr, rationale, __FILE__, __LINE__)

/* Utility function to delimit the start and end of test cases. */
void PrintBanner(const string& header) {
  cout << "\nBeginning test: " << header << endl;
  cout << setw(40) << setfill('-') << "" << setfill(' ') << endl;
}

/* Utility function to signal that a test isn't begin run. */
void TestDisabled(const string& header) {
  cout << "== Test " << header << " NOT RUN: press ENTER to continue ==" << endl;

  /* Pause for the user to hit enter. */
  PressEnterToContinue();
}

/* Utility function to signal the end of a test. */
void EndTest() {
  cout << "== end of test: press ENTER to continue ==" << endl;
  PressEnterToContinue();
}

/* Basic test: Can we build a small tree and look up the elements it contains? */
void BasicKDTreeTest() try {
#if BasicKDTreeTestEnabled
  PrintBanner("Basic KDTree Test");

  /* Construct the KDTree. */
  KDTree<128, size_t> kd;
  CheckCondition(true, "KDTree construction completed.");

  /* Check basic properties of the KDTree. */
  CheckCondition(kd.dimension() == 128, "Dimension is 128.");
  CheckCondition(kd.size() == 0,      "New KD tree has no elements.");
  CheckCondition(kd.empty(),          "New KD tree is empty.");

  /* Add some elements. */

  /* Read Data */
  ReadData rd1("sample_data.txt");
  int row1=rd1.get_num_of_elements();
  int N1=rd1.get_num_of_dimensions();
  vector <vector<double>> dataPointsVec1=rd1.allDataPointsVec;

   double dataPoints[row1][N1];
   for(int i=0; i<row1; i++)
   {
        for(int j=0; j<N1; j++)
        {
             dataPoints[i][j]=dataPointsVec1[i][j];
        }
   }


  for (size_t i = 0; i < row1; ++i)
    kd.insert(PointFromRange<128>(dataPoints[i], dataPoints[i] + 128), i);

  /* Check basic properties again. */
  CheckCondition(kd.size() == 1000, "After adding 1000 elements, KDTree has size 1000.");
  CheckCondition(!kd.empty(),    "After adding  elements, KDTree is not empty.");

  /* Make sure that the elements we built the tree out of are still there. */
  CheckCondition(kd.contains(PointFromRange<128>(dataPoints[0], dataPoints[0] + 128)), "New KD tree has element zero.");
  CheckCondition(kd.contains(PointFromRange<128>(dataPoints[1], dataPoints[1] + 128)), "New KD tree has element one.");
  CheckCondition(kd.contains(PointFromRange<128>(dataPoints[2], dataPoints[2] + 128)), "New KD tree has element two.");

  /* Make sure that the values of these points are correct. */
  for (size_t i = 0; i < row1; ++i)
    CheckCondition(kd.at(PointFromRange<128>(dataPoints[i], dataPoints[i] + 128)) == i, "New KD tree has correct values.");

  EndTest();
#else
  TestDisabled("BasicKDTreeTest");
#endif
} catch (const exception& e) {
  FailTest(e);
}

void NearestNeighborTest() try {
#if NearestNeighborTestEnabled
  PrintBanner("Nearest Neighbor Test");

   /* Add some elements. */

	 /* Read Data */
  ReadData rd1("sample_data.txt");
  int row1=rd1.get_num_of_elements();
  int N1=rd1.get_num_of_dimensions();
  vector <vector<double>> dataPointsVec1=rd1.allDataPointsVec;

   double dataPoints[row1][N1];
   for(int i=0; i<row1; i++)
   {
        for(int j=0; j<N1; j++)
        {
             dataPoints[i][j]=dataPointsVec1[i][j];
        }
   }


  /* Build two data sets - a set of tree points and a set of test points. */

  /* Each test point should be right next to the corresponding point in the
   * tree point set.
   */

   /* Read Data */
  ReadData rd2("query_points.txt");
  int row2=rd1.get_num_of_elements();
  int N2=rd1.get_num_of_dimensions();
  vector <vector<double>> dataPointsVec2=rd1.allDataPointsVec;

   double testPoints[row2][N2];
   for(int i=0; i<row2; i++)
   {
        for(int j=0; j<N2; j++)
        {
             testPoints[i][j]=dataPointsVec2[i][j];
        }
   }

  /* Build up a data set from the first set of points. */
  KDTree<128, size_t> kd;
  for (size_t i = 0; i < row1; ++i)
    kd.insert(PointFromRange<128>(dataPoints[i], dataPoints[i] + 128), i);

  /* Check that calling nearest neighbor on a point in the tree yields that point. */
  for (size_t i = 0; i < row1; ++i)
    CheckCondition(kd.kNNValue(PointFromRange<128>(dataPoints[i], dataPoints[i] + 128), 1) == i, "Nearest neighbor of element is that element.");

  /* Check that calling nearest neighbor on the test points yields the correct tree points. */
  for (size_t i = 0; i < row2; ++i)
    CheckCondition(kd.kNNValue(PointFromRange<128>(testPoints[i], testPoints[i] + 128), 1) == i, "Test point yielded correct nearest neighbor.");

  EndTest();
#else
  TestDisabled("NearestNeighborTest");
#endif
} catch (const exception& e) {
  FailTest(e);
}

int main() {
  /* Step one Tests */
  BasicKDTreeTest();

  /* Step Two Tests */
  NearestNeighborTest();



#if (BasicKDTreeTestEnabled && \
     NearestNeighborTestEnabled)
  cout << "All tests completed!  If they passed, you should be good to go!" << endl << endl;
#else
  cout << "Not all tests were run.  Enable the rest of the tests, then run again." << endl << endl;
#endif

  PressEnterToContinue();
}
