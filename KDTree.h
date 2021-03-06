/**
 * File: KDTree.h
 * Author Lili Meng (lilimeng1103@gmail.com) based on the code by Keith Schwarz (htiek@cs.stanford.edu) and JohnYangSam
 * Thanks a lot for the discussion with Jimmy Chen, Victor Gan, Keith Schwarz.
 * ------------------------
 * Perform constructing trees, efficient exact query for k-nearest neighbors based on Bounded priority queue kd-tree,
 * Best-Bin-First(BBF) query for approximate k-nearest neighbors search.
 * For more BBF query, please refer to
 * Beis, J. S. and Lowe, D. G.  Shape indexing using approximate nearest-neighbor search in high-dimensional spaces.
 *
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 **/

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include <stdexcept>
#include <cmath>
#include <queue>
#include <vector>
#include <assert.h>
#include <stack>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <stdio.h>

#include "Point.h"
#include "BoundedPQueue.h"
#include "ReadData.h"

using namespace std;


/* A utility function to construct a Point from a range of iterators. */
template <size_t N, typename IteratorType>
Point<N> PointFromRange(IteratorType begin, IteratorType end) {
    Point<N> result;
    copy(begin, end, result.begin());
    return result;
}

class distance_index
{
public:
    int index_;
    double distance_;

    distance_index(int index, double distance)
    {
        index_ = index;
        distance_ = distance;
    }

    bool operator < (const distance_index & other) const
    {
        return this->distance_ < other.distance_;
    }

    bool operator > (const distance_index & other) const
    {
        return this->distance_ > other.distance_;
    }

};

template<typename T>
class TreeNode_distance
{
public:
    T node_;
    double distance_;

    TreeNode_distance(T node, double distance)
    {
        node_=node;
        distance_ = distance;
    }

    bool operator < (const TreeNode_distance & other) const
    {
        return this->distance_ < other.distance_;
    }

     bool operator > (const TreeNode_distance & other) const
    {
        return this->distance_ > other.distance_;
    }

};




template <size_t N, typename ElemType>
class KDTree {
public:
    /** Constructor: KDTree();
     * Usage: KDTree<3, int> myTree;
     * ----------------------------------------------------
     * Constructs an empty KDTree.
     **/
    KDTree();

    /** Destructor: ~KDTree()
     * Usage: (implicit)
     * ----------------------------------------------------
     * Cleans up all resources used by the KDTree.
     **/
    ~KDTree();

    /**
     * KDTree(const KDTree& rhs);
     * KDTree& operator=(const KDTree& rhs);
     * Usage: KDTree<3, int> one = two;
     * Usage: one = two;
     * -----------------------------------------------------
     * Deep-copies the contents of another KDTree into this one.
     **/
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);

    /** size_t dimension() const;
     * Usage: size_t dim = kd.dimension();
     * ----------------------------------------------------
     * Returns the dimension of the points stored in this KDTree.
     **/
    size_t dimension() const;

    /**size_t size() const;
     * bool empty() const;
     * Usage: if (kd.empty())
     * ----------------------------------------------------
     * Returns the number of elements in the kd-tree and whether the tree is
     * empty.
     **/
    size_t size() const;
    bool empty() const;

    /** bool contains(const Point<N>& pt) const;
     * Usage: if (kd.contains(pt))
     * ----------------------------------------------------
     * Returns whether the specified point is contained in the KDTree.
     **/
    bool contains(const Point<N>& pt) const;

    /**void insert(const Point<N>& pt, const ElemType& value);
     * Usage: kd.insert(v, "This value is associated with v.");
     * ----------------------------------------------------
     * Inserts the point pt into the KDTree, associating it with the specified
     * value. If the element already existed in the tree, the new value will
     * overwrite the existing one.
     **/

    void insert(const Point<N>& pt, const ElemType& value);

    /** Build Tree*/
    void buildTree(const vector<vector<double>>& dataPointsVec, const vector<ElemType>& indices);

    /** ElemType& operator[](const Point<N>& pt);
     * Usage: kd[v] = "Some Value";
     * ----------------------------------------------------
     * Returns a reference to the value associated with point pt in the KDTree.
     * If the point does not exist, then it is added to the KDTree using the
     * default value of ElemType as its key.
     **/
    ElemType& operator[](const Point<N>& pt);

    /** ElemType& at(const Point<N>& pt);
     * const ElemType& at(const Point<N>& pt) const;
     * Usage: cout << kd.at(v) << endl;
     * ----------------------------------------------------
     * Returns a reference to the key associated with the point pt. If the point
     * is not in the tree, this function throws an out_of_range exception.
     **/
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;

    /** function for saving the KD tree to file **/
    bool save_tree_to_file(const char* fileName);

    /** vector<ElemType> kNNValue(const Point<N>& key, size_t k) const
     * Usage: k nearest ;
     * ----------------------------------------------------
     * Given a point key and an integer k, finds the k nearest neighbors in the KDTree
     * nearest to key
     **/
   vector<ElemType> kNNValues(const Point<N>& key, size_t k) const;

   /** vector<ElemType> BBFkNNValue(const Point<N>& key, size_t k) const
     * Usage: approximate k nearest neigbhour ;
     * ----------------------------------------------------
     * Given a point key and an integer k, finds the k nearest neighbors in the KDTree
     * nearest to key
     **/
    vector<ElemType> BBFkNNValues(const Point<N>& query_point, size_t k, size_t maxEpoch) const;


private:

  /** implementation details **/

    struct TreeNode {

        Point<N> key;
        ElemType value;
        size_t level;
        //int leaf;   /** 1 if node is a leaf, 0 otherwise */

        TreeNode* left;
        TreeNode* right;
    };

    typedef priority_queue<TreeNode_distance<TreeNode*>, vector<TreeNode_distance<TreeNode*>>, greater<TreeNode_distance<TreeNode*>>> NodeMinPQ;

    TreeNode* root;

    /** The number of elements currently stored **/
    size_t numElements;

    /** Takes in a node and recursively delete the subtree it represents
     going from its children up **/
    void deleteTreeNode(TreeNode* currentTreeNode);

    /** Helper function for copying KDTree **/
    TreeNode* copyTree(TreeNode* rootNode);

    /** Helper function for saving tree to disk **/
    void save_tree_helper(FILE *pf, TreeNode* node);

    /** KNNValue helper function of building BoundedPQueue for the KNNValue  */
    void KNNValueshelper(const Point<N>&key, BoundedPQueue<TreeNode*>& kNearestPQ, TreeNode* currentNode) const;

    /** function for traversing to the leaf＊*/
    TreeNode* bbfExploreToLeaf(const Point<N>& pt, TreeNode* root, NodeMinPQ& minPQ) const;

};

/** KDTree class implementation details */

/** Constructor **/
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
    numElements = 0;
    root = NULL;
}

/** Destructor **/
template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    deleteTreeNode(root);
    numElements = 0;
}

/** Helper function of deleting the current TreeNode**/
template <size_t N, typename ElemType>
void KDTree<N, ElemType>::deleteTreeNode(TreeNode* currentNode){

    if(currentNode == NULL) return;
    /**Recursion**/
    deleteTreeNode(currentNode->left);
    deleteTreeNode(currentNode->right);
    delete currentNode;
}

/**Copy Constructor **/
template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs) {

    root = copyTree(rhs.root);
    numElements = rhs.numElements;
}

/** Assignment operator, clears old tree if not the same tree and copies
 *  the 'other' tree into the new tree
 */
template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs){

    if (this !=&rhs)
    {
        deleteTreeNode(this->root);
        root=copyTree(rhs.root);
        numElements = rhs.numElements;
    }

    return *this;
}

/** CopyTree **/
template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::TreeNode* KDTree<N, ElemType>::copyTree(TreeNode* rootNode){

    if(rootNode==NULL) return NULL;

    TreeNode* rootNodeCopy = new TreeNode;

    rootNodeCopy->key = rootNode->key;
    rootNodeCopy->value = rootNode->value;
    rootNodeCopy->level = rootNode->level;

    /** Recursion**/
    rootNodeCopy->left = copyTree(rootNode->left);
    rootNodeCopy->right = copyTree(rootNode->right);

    return rootNodeCopy;
}


template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {

    return N;
}

/** Return the number of elements currently stored **/
template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const{

    return numElements;
}

/** Returns whether the it's empty**/
template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const{

    if(numElements==0)
    return true;
    else
    return false;
}

/**
 * Returns whether the specified Point is contained in the KDTree
 **/
template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const{

    TreeNode* currentNode = root;
    while(currentNode!=NULL)
    {
        if(currentNode->key==pt)
        {
            return true;
        }

        /** compares the parts of the pt to determine which subtree to look into next, image N=2, that is binary search tree **/
        if(pt[currentNode->level % N] >= currentNode->key[currentNode->level %N])
        {
            currentNode = currentNode->right;
        }
        else
        {
            currentNode = currentNode->left;
        }

    }

    return false;
}


/** Insert the specified point into the KDTree with associated value. If the point already exists in the KDTree, the old value is overwritten **/
template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {


    TreeNode* currentNode = root;
    TreeNode* prevNode = NULL;
    size_t level = 0;

    while(currentNode!=NULL)
    {
        ++level;
        if(pt==currentNode->key)
        {
            currentNode->value = value;
        }

        size_t dim = currentNode->level % N;

        // If pt[dim] >= currentNode->key[dim],insert the node to the right node.

        if(pt[dim] < currentNode->key[dim]){
            prevNode = currentNode;
            currentNode = currentNode->left;
        }
        else {
            prevNode = currentNode;
            currentNode = currentNode->right;
        }
    }


    ++numElements;

    TreeNode* newNode = new TreeNode;
    newNode->key = pt;
    newNode->value = value;
    newNode->left = NULL;
    newNode->right = NULL;

    if(currentNode == root){
        root = newNode;
    }
    else {
        if(pt[prevNode->level % N]<prevNode->key[prevNode->level % N])
        {
            prevNode->left = newNode;
        }
        else
        {
            prevNode->right = newNode;
        }
    }
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::buildTree(const vector<vector<double>>& dataPointsVec, const vector<ElemType>& indices)
{
    int numberOfElements = dataPointsVec.size();
    int dim = dataPointsVec[1].size();

    double dataPoints[numberOfElements][dim];
    for(int i=0; i<numberOfElements; i++)
    {
        for(int j=0; j<dim; j++)
        {
             dataPoints[i][j]=dataPointsVec[i][j];
        }
    }

     for(int i=0; i<dataPointsVec.size(); i++)
    {

       insert(PointFromRange<N>(dataPoints[i], dataPoints[i] + N), indices[i]);

    }
}


/** Returns a reference to the value associated with the point pt. If the point does not exist in the KDTree, it is added with
 * the default value of ElemType as its value, and a reference to this new value is returned. This is the same behavior as the
 * STL map's operator[]
 * Note that this function does not have a const overload because the function may mutate the tree
 */
template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt){

    TreeNode* currentNode = root;
    TreeNode* prevNode = NULL;
    size_t level = 0;
    while (currentNode != NULL){
        ++level;

        if(pt==currentNode->key){
            return currentNode->value;
        }

        size_t dim = currentNode->level % N;

        if(pt[dim] < currentNode->key[dim]) {
            prevNode = currentNode;
            currentNode = currentNode->left;
        }
        else {
            prevNode = currentNode;
            currentNode = currentNode->right;
        }
    }

    ++numElements;

    //Make the new node to insert into the KDTree
    TreeNode* newNode = new TreeNode;
    newNode->key = pt;
    newNode->value = ElemType();
    newNode->level = level;
    newNode->left = NULL;
    newNode->right = NULL;

    if(currentNode == root) {
        root = newNode;
        return newNode->value;
    }
    else {
        if(pt[prevNode->level % N] >= prevNode->key[prevNode->level % N])
        {
            prevNode->right = newNode;
        }
        else
        {
            prevNode->left = newNode;
        }

        return newNode->value;
    }

}

/** Returns a reference to the value associated with the point pt, if it exists. If the point is not in the tree,
 * then this function throws an out_of_range exception
 **/
template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt){

    TreeNode* currentNode = root;
    while (currentNode != NULL) {
        if(currentNode->key == pt) {
            return currentNode->value;
        }

        //Find the dim
        size_t dim = currentNode->level % N;

        //compare the approximate indices
        if(pt[dim] < currentNode->key[dim]) {
            currentNode = currentNode->left;
        }
        else {
            currentNode = currentNode->right;
        }
    }

    throw out_of_range("That point does not exist");
}

/** This function is const-overloaded, since it does not change the tree **/
template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {

     TreeNode* currentNode = root;
     while(currentNode != NULL) {
        if(currentNode->key == pt) {
            return currentNode->value;
        }

        size_t dim = currentNode->level % N;

        if(pt[dim] < currentNode ->key[dim]) {
            currentNode = currentNode->left;
        }
        else {
            currentNode = currentNode->right;
        }
     }

    throw out_of_range("That point does not exist");

}

template <size_t N, typename ElemType>
void  KDTree<N, ElemType>::save_tree_helper(FILE *pf, TreeNode* node)
{
    if (!node) {
        fprintf(pf, "#\n");
        return;
    }
    size_t dim = node->level % N;
    // write current node
    fprintf(pf, "%u\t %lf\n", node->level%N, node->key[dim]);

    save_tree_helper(pf, node->left);
    save_tree_helper(pf, node->right);
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::save_tree_to_file(const char* fileName)
{
    assert(root);
    FILE *pf = fopen(fileName, "w");
    if (!pf) {
        cout<<"can not open file "<<fileName<<endl;
        return false;
    }

    fprintf(pf, "level\t split_value_\n");

    save_tree_helper(pf, root);
    fclose(pf);
    return true;
}

/** Exact Query of K Nearest-Neigbour **/
template <size_t N, typename ElemType>
vector<ElemType> KDTree<N, ElemType>::kNNValues(const Point<N>& key, size_t k) const {

    BoundedPQueue<TreeNode*> kNearestPQ(k);
    KNNValueshelper(key, kNearestPQ, root);

    vector<ElemType> kNNValues;
    while(!kNearestPQ.empty())
    {
        kNNValues.push_back((kNearestPQ.dequeueMin())->value);
    }
    return kNNValues;
}

/**
 * KNNValueshelper(key, bpq, currentNode)
 * key--query point
 * A helper function of building the bounded priority queue of points nearest to the query point in the KDTree
 **/

template<size_t N, typename ElemType>
void KDTree<N, ElemType>::KNNValueshelper(const Point<N>& key, BoundedPQueue<TreeNode*>& kNearestPQ, TreeNode* currentNode) const {

    if (currentNode == NULL) return;

    kNearestPQ.enqueue(currentNode, Distance(currentNode->key, key));

    size_t dim = currentNode->level % N;

    //like Binary Search Tree, if key[dim] < currentNode->key[dim], turn to the left of the tree
    if(key[dim] < currentNode->key[dim])
    {
        KNNValueshelper(key, kNearestPQ, currentNode->left);

        // If the query hypersphere crosses the splitting plane, check the other subtree
        if(kNearestPQ.size() < kNearestPQ.maxSize() || fabs(currentNode->key[dim] - key[dim]) < kNearestPQ.worst() ) {

            KNNValueshelper(key, kNearestPQ, currentNode->right);
        }
    }
    else //like Binary Search Tree, if key[dim] >= currentNode->key[dim], turn to the right of the tree
    {
        KNNValueshelper(key, kNearestPQ, currentNode->right);

        // If the hypersphere crosses the splitting plane, check the other subtree
        if(kNearestPQ.size() < kNearestPQ.maxSize() || fabs(currentNode->key[dim] - key[dim]) < kNearestPQ.worst() ) {

            KNNValueshelper(key, kNearestPQ, currentNode->left);
        }

    }
}



/** Function for traversing the tree**/
template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::TreeNode* KDTree<N, ElemType>::bbfExploreToLeaf(const Point<N>& pt, TreeNode* root, NodeMinPQ& untraversed_minPQ) const
{
    TreeNode* untraversed;   // untraversed storing the untraversed TreeNode* on the KD Tree
    TreeNode* currentNode = root;

    double value;
    size_t dim;

    while(currentNode!=NULL && currentNode->left!=NULL && currentNode->right!=NULL)   //currentNode->left!=NULL && currentNode->right!=NULL signifies that currentNode is not a leaf
    {
        //partition dimension and value;
        dim = currentNode->level % N;
        value = currentNode->value;


        // go to a child and preserve the other
        if(pt[dim] < currentNode->key[dim])
        {
            untraversed = currentNode;
            currentNode = currentNode->left;
        }
        else
        {
            untraversed = currentNode->left;
            currentNode = currentNode->right;
        }

        if(untraversed!=NULL)
        {
            TreeNode_distance<TreeNode*> di(untraversed, fabs(pt[untraversed->dim]-untraversed->key[dim]));
            untraversed_minPQ.push(di);
        }


    }

    return currentNode;
}

/**
    * Search for approximate k nearest neighbours using Best-Bin-First (BBF) approach
    * @param key        Query point data
    * @param k          number of nearest neighbour returned
    * @param maxEpoch   maximum search epoch

template <size_t N, typename ElemType>
vector<ElemType> KDTree<N, ElemType>::BBFkNNValues(const Point<N>& query_point, size_t k, size_t maxEpoch) const{

    BoundedPQueue<TreeNode*> kNearestPQ(k);

    size_t epoch = 0;

    double max_distance = numeric_limits<double>::max();

    TreeNode* cur_node = root;

    NodeMinPQ priority_unexplored_points;

    TreeNode_distance<TreeNode*> di(root, 0);
    priority_unexplored_points.push(di);

    priority_queue<distance_index> priority_points;

    while(!priority_unexplored_points.empty() && epoch < maxEpoch)
    {

        bbfExploreToLeaf(query_point, cur_node, priority_unexplored_points);

        kNearestPQ.enqueue(priority_unexplored_points.top().node_, priority_unexplored_points.top().distance_);

        ++epoch;

    }


    vector<ElemType> bbfkNNValues;
    while(!kNearestPQ.empty())
    {
         kNNValues.push_back((kNearestPQ.dequeueMin())->value);
    }
    return bbfkNNValues;
}*/


#endif // KDTREE_INCLUDED
