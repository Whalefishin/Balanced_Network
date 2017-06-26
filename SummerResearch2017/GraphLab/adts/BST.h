#pragma once

#include <utility>
#include <vector>
#include "dictionary.h"

using std::pair;
using std::vector;

/**
 * The interface for a binary search tree.  Note that it is a template
 * definition based on a type for key (K) and type for value (V).  A BST is a
 * dictionary structure; that is, elements are indexed by key and all keys must
 * be unique.  This interface therefore inherits from the Dictionary interface.
 * @tparam K The type of keys in the BST.
 * @tparam V The type of values in the BST.
 */
template <typename K, typename V>
class BST : public Dictionary<K,V> {
public:
    virtual ~BST() {};

    /**
     * Returns the size of the BST.
     * @return The number of key-value pairs in the data structure.
     */
    virtual int getSize() = 0;

    /**
     * Returns true if the tree is empty.
     * @return true if there are no elements in the BST.
     */
    virtual bool isEmpty() = 0;

    /**
     * Inserts the key-value pair into the tree
     * @param key The key for the new mapping.
     * @param value The value to associate with that key.
     * @throws runtime_error If they key already exists.
     */
    virtual void insert(K key, V value) = 0;

    /**
     * Finds the element indexed by the given key and updates its value to the
     * provided value parameter.
     * @param key The key of the mapping to update.
     * @param value The new value to associate with that key.
     * @throws runtime_error if the key is not found in the BST.
     */
    virtual void update(K key, V value) = 0;

    /**
     * Returns the value associated with the given key
     * @param key The key of the mapping to find.
     * @return The value associated with that key.
     * @throws runtime_error If the key is not found in the BST.
     */
    virtual V get(K key) = 0;

    /**
     * Determines if a given key exists in a mapping in this BST.
     * @param key The key to look for.
     * @return true if a mapping in this BST has that key.
     */
    virtual bool contains(K key) = 0;

    /**
     * Deletes the element with given key from the tree.
     * @param key The key of the mapping to remove.
     * @throws runtime_error If they key was not already in this BST.
     */
    virtual void remove(K key) = 0;

    /**
     * Obtains a vector containing all keys in this BST.
     * @return An STL vector object containing the keys of every mapping in this
     *         BST (in no particular order).
     */
    virtual vector<K> getKeys() = 0;

    /**
     * Obtains a vector containing all key-value pairs in this BST.
     * @return An STL vector object containing the key-value pairs in every
     *         mapping in this BST (in no particular order).
     */
    virtual vector<pair<K,V>> getItems() = 0;

    ////////////////////////////////////////////////////////////////////////////
    // The following methods are unique to the BST interface

    /**
     * Returns a height for the tree (i.e., largest depth for any leaf node).
     * @return The height of this tree (or -1 if the tree is empty).
     */
    virtual int getHeight() = 0;

    /**
    * Returns the smallest key in this tree.
    * @return The minimum key in the BST.
    * @throws runtime_error If this BST is empty.
    */
    virtual K getMinKey() = 0;

    /**
     * Returns the largest key in this tree.
     * @return The maximum key in the BST.
     * @throws runtime_error If this BST is empty.
     */
    virtual K getMaxKey() = 0;

    /**
     * Obtains a pre-order traversal of the key-value pairs in this BST.
     * @return An STL vector containing all key-value pairs in this BST.  This
     *         vector is guaranteed to return the elements in a pre-order
     *         traversal.
     */
    virtual vector<pair<K,V>> traversePreOrder() = 0;

    /**
     * Obtains a in-order traversal of the key-value pairs in this BST.
     * @return An STL vector containing all key-value pairs in this BST.  This
     *         vector is guaranteed to return the elements in a in-order
     *         traversal.
     */
    virtual vector<pair<K,V>> traverseInOrder() = 0;

    /**
     * Obtains a post-order traversal of the key-value pairs in this BST.
     * @return An STL vector containing all key-value pairs in this BST.  This
     *         vector is guaranteed to return the elements in a post-order
     *         traversal.
     */
    virtual vector<pair<K,V>> traversePostOrder() = 0;

    /**
     * Obtains a level-order traversal of the key-value pairs in this BST.
     * @return An STL vector containing all key-value pairs in this BST.  This
     *         vector is guaranteed to return the elements in a level-order
     *         traversal.
     */
    virtual vector<pair<K,V>> traverseLevelOrder() = 0;
};

