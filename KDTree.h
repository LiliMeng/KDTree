/**
 * File: KDTree.h
 * Author Lili Meng (lilimeng1103@gmail.com) based on the starter code by Keith Schwarz (htiek@cs.stanford.edu)
 * ------------------------
 * Perform constructing trees, efficient exact query for k-nearest neighbors based on Bounded priority queue kd-tree,
 * Best-Bin-First query for k-nearest points based on BPQ.
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 **/

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include "Point.h"
#include "BoundedPQueue.h"
#include <stdexcept>
#include <cmath>

using namespace std;

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

    /** ElemType kNNValue(const Point<N>& key, size_t k) const
     * Usage: cout << kd.kNNValue(v, 3) << endl;
     * ----------------------------------------------------
     * Given a point v and an integer k, finds the k points in the KDTree
     * nearest to v and returns the most common value associated with those
     * points. In the event of a tie, one of the most frequent value will be
     * chosen.
     **/
    ElemType kNNValue(const Point<N>& key, size_t k) const;

    /*
    ElemType BBFkNNValue(const Point<N>& key, size_t k) const;
    */
private:
    /** implementation details **/

    struct TreeNode {

        Point<N> key;
        ElemType value;
        size_t level;

        TreeNode* left;
        TreeNode* right;
    };

    TreeNode* root;

    /** The number of elements currently stored **/
    size_t numElements;

    /** Takes in a node and recursively delete the subtree it represents
     going from its children up **/
    void deleteTreeNode(TreeNode* currentTreeNode);

    /** Helper function for copying KDTree **/
    TreeNode* copyTree(TreeNode* rootNode);

    /** Recursive helper function for the KNNValue function */
    void KNNValueRecurse(const Point<N>&key, BoundedPQueue<TreeNode*>& nearestPQ, TreeNode* currentNode) const;

    /** A helper function that returns the most commonly occuring value
     stored in the Nodes of a Node* PQ **/
    ElemType FindMostCommonValueInPQ(BoundedPQueue<TreeNode*> nearestPQ) const;

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

        size_t keyIndex = currentNode->level % N;

        if(pt[keyIndex] < currentNode->key[keyIndex]){
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
        if(pt[prevNode->level % N]>=prevNode->key[prevNode->level % N])
        {
            prevNode->right = newNode;
        }
        else
        {
            prevNode->left = newNode;
        }
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

        size_t keyIndex = currentNode->level % N;

        if(pt[keyIndex] < currentNode->key[keyIndex]) {
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

        //Find the keyIndex
        size_t keyIndex = currentNode->level % N;

        //compare the approximate indices
        if(pt[keyIndex] < currentNode->key[keyIndex]) {
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

        size_t keyIndex = currentNode->level % N;

        if(pt[keyIndex] < currentNode ->key[keyIndex]) {
            currentNode = currentNode->left;
        }
        else {
            currentNode = currentNode->right;
        }
     }

    throw out_of_range("That point does not exist");

}

/** Exact Query of K Nearest-Neigbour **/
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& key, size_t k) const {

    BoundedPQueue<TreeNode*> nearestPQ(k);
    KNNValueRecurse(key, nearestPQ, root);

    return FindMostCommonValueInPQ(nearestPQ);
}

/*
 * KNNValueRecurse(pt, bpq, currentNode)
 * A recursive helper function which builds a bounded priority queue of the points
 * nearest to the entered point in the KDTree
 */
template<size_t N, typename ElemType>
void KDTree<N, ElemType>::KNNValueRecurse(const Point<N>& key, BoundedPQueue<TreeNode*>& nearestPQ, TreeNode* currentNode) const {

    if (currentNode == NULL) return;

    nearestPQ.enqueue(currentNode, Distance(currentNode->key, key));

    //Recursion
    size_t keyIndex = currentNode->level % N;
    if(key[keyIndex] < currentNode->key[keyIndex]) {
        KNNValueRecurse(key, nearestPQ, currentNode->left);

        //If the hypersphere crosses the splitting plane check the other subtree
        if(nearestPQ.size()!=nearestPQ.maxSize() || fabs(currentNode->key[keyIndex] - key[keyIndex]) < nearestPQ.worst() ) {
            KNNValueRecurse(key, nearestPQ, currentNode->right);
        }
    }
    else {

        KNNValueRecurse(key, nearestPQ, currentNode->right);

        //If the hypersphere crosses the splitting plane check the other subtree
        if((nearestPQ.size() != nearestPQ.maxSize()) || fabs(currentNode->key[keyIndex] - key[keyIndex]) < nearestPQ.worst() ){
            KNNValueRecurse(key, nearestPQ, currentNode->left);
        }

    }

}

/*
 * FindMostCommonValueInPQ(bpq)
 * Takes in a bounded priority queue of Node* in the KDTree and
 * returns the most common value stored in the nodes.
 */
template<size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::FindMostCommonValueInPQ(BoundedPQueue<TreeNode*> nearestPQ) const{
    multiset<ElemType> values;
    while(!nearestPQ.empty()) {
        values.insert((nearestPQ.dequeueMin())->value);
    }

    ElemType best;
    size_t bestFrequency = 0;
    for(typename multiset<ElemType>::iterator iter = values.begin(); iter!=values.end(); ++iter){
        if(values.count(*iter) > bestFrequency) {
            best = *iter;
            bestFrequency = values.count(*iter);
        }
    }
   return best;
}

/*
/** Best-Bin-First Query for K Nearest-Neigbour
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::BBFkNNValue(const Point<N>& key, size_t k) const {

    BoundedPQueue<TreeNode*> nearestPQ(k);
    KNNValueRecurse(key, nearestPQ, root);

    return FindMostCommonValueInPQ(nearestPQ);

}**/

#endif // KDTREE_INCLUDED/**
 * File: KDTree.h
 * Author Lili Meng (lilimeng1103@gmail.com) based on the starter code by Keith Schwarz (htiek@cs.stanford.edu)
 * ------------------------
 * Perform constructing trees, efficient exact query for k-nearest neighbors based on Bounded priority queue kd-tree,
 * Best-Bin-First query for k-nearest points based on BPQ.
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 **/

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include "Point.h"
#include "BoundedPQueue.h"
#include <stdexcept>
#include <cmath>

using namespace std;

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

    /** ElemType kNNValue(const Point<N>& key, size_t k) const
     * Usage: cout << kd.kNNValue(v, 3) << endl;
     * ----------------------------------------------------
     * Given a point v and an integer k, finds the k points in the KDTree
     * nearest to v and returns the most common value associated with those
     * points. In the event of a tie, one of the most frequent value will be
     * chosen.
     **/
    ElemType kNNValue(const Point<N>& key, size_t k) const;

    /*
    ElemType BBFkNNValue(const Point<N>& key, size_t k) const;
    */
private:
    /** implementation details **/

    struct TreeNode {

        Point<N> key;
        ElemType value;
        size_t level;

        TreeNode* left;
        TreeNode* right;
    };

    TreeNode* root;

    /** The number of elements currently stored **/
    size_t numElements;

    /** Takes in a node and recursively delete the subtree it represents
     going from its children up **/
    void deleteTreeNode(TreeNode* currentTreeNode);

    /** Helper function for copying KDTree **/
    TreeNode* copyTree(TreeNode* rootNode);

    /** Recursive helper function for the KNNValue function */
    void KNNValueRecurse(const Point<N>&key, BoundedPQueue<TreeNode*>& nearestPQ, TreeNode* currentNode) const;

    /** A helper function that returns the most commonly occuring value
     stored in the Nodes of a Node* PQ **/
    ElemType FindMostCommonValueInPQ(BoundedPQueue<TreeNode*> nearestPQ) const;

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

        size_t keyIndex = currentNode->level % N;

        if(pt[keyIndex] < currentNode->key[keyIndex]){
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
        if(pt[prevNode->level % N]>=prevNode->key[prevNode->level % N])
        {
            prevNode->right = newNode;
        }
        else
        {
            prevNode->left = newNode;
        }
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

        size_t keyIndex = currentNode->level % N;

        if(pt[keyIndex] < currentNode->key[keyIndex]) {
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

        //Find the keyIndex
        size_t keyIndex = currentNode->level % N;

        //compare the approximate indices
        if(pt[keyIndex] < currentNode->key[keyIndex]) {
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

        size_t keyIndex = currentNode->level % N;

        if(pt[keyIndex] < currentNode ->key[keyIndex]) {
            currentNode = currentNode->left;
        }
        else {
            currentNode = currentNode->right;
        }
     }

    throw out_of_range("That point does not exist");

}

/** Exact Query of K Nearest-Neigbour **/
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& key, size_t k) const {

    BoundedPQueue<TreeNode*> nearestPQ(k);
    KNNValueRecurse(key, nearestPQ, root);

    return FindMostCommonValueInPQ(nearestPQ);
}

/*
 * KNNValueRecurse(pt, bpq, currentNode)
 * A recursive helper function which builds a bounded priority queue of the points
 * nearest to the entered point in the KDTree
 */
template<size_t N, typename ElemType>
void KDTree<N, ElemType>::KNNValueRecurse(const Point<N>& key, BoundedPQueue<TreeNode*>& nearestPQ, TreeNode* currentNode) const {

    if (currentNode == NULL) return;

    nearestPQ.enqueue(currentNode, Distance(currentNode->key, key));

    //Recursion
    size_t keyIndex = currentNode->level % N;
    if(key[keyIndex] < currentNode->key[keyIndex]) {
        KNNValueRecurse(key, nearestPQ, currentNode->left);

        //If the hypersphere crosses the splitting plane check the other subtree
        if(nearestPQ.size()!=nearestPQ.maxSize() || fabs(currentNode->key[keyIndex] - key[keyIndex]) < nearestPQ.worst() ) {
            KNNValueRecurse(key, nearestPQ, currentNode->right);
        }
    }
    else {

        KNNValueRecurse(key, nearestPQ, currentNode->right);

        //If the hypersphere crosses the splitting plane check the other subtree
        if((nearestPQ.size() != nearestPQ.maxSize()) || fabs(currentNode->key[keyIndex] - key[keyIndex]) < nearestPQ.worst() ){
            KNNValueRecurse(key, nearestPQ, currentNode->left);
        }

    }

}

/*
 * FindMostCommonValueInPQ(bpq)
 * Takes in a bounded priority queue of Node* in the KDTree and
 * returns the most common value stored in the nodes.
 */
template<size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::FindMostCommonValueInPQ(BoundedPQueue<TreeNode*> nearestPQ) const{
    multiset<ElemType> values;
    while(!nearestPQ.empty()) {
        values.insert((nearestPQ.dequeueMin())->value);
    }

    ElemType best;
    size_t bestFrequency = 0;
    for(typename multiset<ElemType>::iterator iter = values.begin(); iter!=values.end(); ++iter){
        if(values.count(*iter) > bestFrequency) {
            best = *iter;
            bestFrequency = values.count(*iter);
        }
    }
   return best;
}

/*
/** Best-Bin-First Query for K Nearest-Neigbour
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::BBFkNNValue(const Point<N>& key, size_t k) const {

    BoundedPQueue<TreeNode*> nearestPQ(k);
    KNNValueRecurse(key, nearestPQ, root);

    return FindMostCommonValueInPQ(nearestPQ);

}**/

#endif // KDTREE_INCLUDED
