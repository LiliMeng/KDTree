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
#include "KDTree_math.h"
#include <stdexcept>
#include <cmath>
#include <queue>

using namespace std;

/// key value pair. Not use the std::pair because there is default
/// "<" operator overloading for it and it is too heavy for operation.
template <class T, class V = double>
struct KeyValue
{
    T key;
    V value;
    KeyValue(const T k, const V v) : key(k), value(v){}
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


    ElemType BBFkNNValue(const Point<N>& key, size_t k, size_t max_epoch) const;

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

     // typedef to avoid ugly long declaration
    //typedef stack<TreeNode*> NodeStack;
    typedef KeyValue<TreeNode*> NodeBind;
    typedef priority_queue<NodeBind, vector<NodeBind>, greater<NodeBind> > NodeMinPQ;
    typedef KeyValue<Point<N>> PointBind;
    typedef priority_queue<PointBind, vector<PointBind> > PointMaxPQ;

    TreeNode* root;

    /** The number of elements currently stored **/
    size_t numElements;

    /** Takes in a node and recursively delete the subtree it represents
     going from its children up **/
    void deleteTreeNode(TreeNode* currentTreeNode);

    /** Helper function for copying KDTree **/
    TreeNode* copyTree(TreeNode* rootNode);

    /** function for traversing to the leaf＊*/
    TreeNode* traverse_to_leaf(Point<N> pt, TreeNode* root, NodeMinPQ& container);

    /** KNNValue helper function of building BoundedPQueue for the KNNValue  */
    void KNNValuehelper(const Point<N>&key, BoundedPQueue<TreeNode*>& nearestPQ, TreeNode* currentNode) const;

   /**
    A helper function that returns the best value
     stored in the Nodes of a TreeNode* BPQ **/
    ElemType FindBestValueInBPQ(BoundedPQueue<TreeNode*> nearestPQ) const;

    /**
     BBFKNNValue helper function of building BoundedPQueue for BBFKNNValue
    void BBFKNNValuehelper(const Point<N>& key, BoundedPQueue<TreeNode*>& nearestPQ, TreeNode* currentNode) const;

     A helper function that returns the best value stored in BBF BPQ
    ElemType BBFFindBestValueINBPQ(BoundedPQueue<TreeNode*> nearestPQ);
    **/

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
    KNNValuehelper(key, nearestPQ, root);

    return FindBestValueInBPQ(nearestPQ);
}


/**
 * KNNValuehelper(key, bpq, currentNode)
 * key--query point
 * A helper function of building the bounded priority queue of points nearest to the query point in the KDTree
 **/

template<size_t N, typename ElemType>
void KDTree<N, ElemType>::KNNValuehelper(const Point<N>& key, BoundedPQueue<TreeNode*>& nearestPQ, TreeNode* currentNode) const {

    if (currentNode == NULL) return;

    nearestPQ.enqueue(currentNode, Distance(currentNode->key, key));

    size_t dim = currentNode->level % N;

    //like Binary Search Tree, if key[dim] < currentNode->key[dim], turn to the left of the tree
    if(key[dim] < currentNode->key[dim])
    {
        KNNValuehelper(key, nearestPQ, currentNode->left);

        // If the query hypersphere crosses the splitting plane, check the other subtree
        if(nearestPQ.size() < nearestPQ.maxSize() || fabs(currentNode->key[dim] - key[dim]) < nearestPQ.worst() ) {

            KNNValuehelper(key, nearestPQ, currentNode->right);
        }
    }
    else //like Binary Search Tree, if key[dim] >= currentNode->key[dim], turn to the right of the tree
    {
        KNNValuehelper(key, nearestPQ, currentNode->right);

        // If the hypersphere crosses the splitting plane, check the other subtree
        if(nearestPQ.size() < nearestPQ.maxSize() || fabs(currentNode->key[dim] - key[dim]) < nearestPQ.worst() ) {

            KNNValuehelper(key, nearestPQ, currentNode->left);
        }

    }
}

/**
 * FindBestValueInBPQ
 * Takes in a bounded Priority Queue of TreeNode* in the KDTree and
 * returns the best value
 */
template<size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::FindBestValueInBPQ(BoundedPQueue<TreeNode*> nearestPQ) const {

    multiset<ElemType> values;
    while(!nearestPQ.empty()) {
        values.insert((nearestPQ.dequeueMin())->value);
    }

    ElemType bestMatch;
    size_t bestOccurence = 0;
    for(auto iter = values.begin(); iter!=values.end(); ++iter) {
        if(values.count(*iter) > bestOccurence)
        {
            bestMatch = * iter;
            bestOccurence = values.count(*iter);
        }
    }

    return bestMatch;
}

/** Function for traversing the tree **/
template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::TreeNode* KDTree<N, ElemType>::traverse_to_leaf(Point<N> pt, TreeNode* root, NodeMinPQ& minPQ)
{
        TreeNode* other;
        TreeNode* currentNode = root;

        double value=currentNode->value;
        size_t dim = currentNode->level % N;

        while(currentNode!=NULL && currentNode->left==NULL && currentNode->right==NULL)
        {
            // partition dimension and value
            dim = currentNode->level % N;
            value = currentNode->value;

            // go to a child and preserve the other
            if(pt[dim] <= value)
            {
                currentNode = currentNode->left;
                other = currentNode->right;
            }
            else
            {
                currentNode = currentNode->right;
                other = currentNode->left;
            }

            if(other!=NULL)
            {
                minPQ.push(NodeBind(other, abs(other->value - pt[other->level % N])));
            }
        }

        return currentNode;
    }

/**
     * Search for approximate k nearest neighbours using the
     * Best Bin First approach.
     *
     * @param Point<N>   query point data in array form
     * @param k          number of nearest neighbour returned
     * @param max_epoch  maximum of epoch of search
     *
     * @return
     */
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::BBFkNNValue(const Point<N>& pt, size_t k, size_t max_epoch) const {

    ElemType result;

    size_t epoch = 0;

    TreeNode* node;

    // checklist for backtracking;
    NodeMinPQ minPQ;

    // min priority queue to keep top k largest
    PointMaxPQ max_pq;

    double cur_best = numeric_limits<double>::max();

    double dist = 0;

    minPQ.push(NodeBind(this->root, 0));

    while(!minPQ.empty() && epoch < max_epoch)
    {
        node = minPQ.top().key;
        minPQ.pop();

        // find leaf and push unprocessed to minPQ
        node = this->traverse_to_leaf(pt, node, minPQ);

        for(size_t i=0; i< node->n; ++i)
        {
            dist = euclidean(node->points[i].data, pt,
                                       N,false);

            if(dist < cur_best)
            {
                // maintain the bounded min priority queue
                if(max_pq.size()==k)
                {
                    //update current best
                    max_pq.pop();
                    max_pq.push(KeyValue<Point<N>>(node->pt[i], dist));
                    cur_best = max_pq.top().value;
                }

               // the special point here is that we need to set best
               // distance to the distance value of largest smallest
               // feature
                else if(max_pq.size() == k-1)
                {
                    max_pq.push(KeyValue<Point<N>>(node->features[i], dist));
                    cur_best = max_pq.top().value;
                }
                else
                {
                    max_pq.push(KeyValue<Point<N>>(node->features[i], dist));
                }
            }
        }
            ++epoch;
    }

        // finally pass results to returned result
        const size_t detected = max_pq.size();
        for(size_t i = 0; i < detected ; ++i)
        {
            result.push_back(max_pq.top().key);
            max_pq.pop();
        }

    return result;
}

#endif // KDTREE_INCLUDED
