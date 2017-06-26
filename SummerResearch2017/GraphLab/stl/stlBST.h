#pragma once

#include <map>
#include <stdexcept>
#include <vector>

#include "../adts/dictionary.h"

using std::map;
using std::runtime_error;
using std::vector;

template <typename K, typename V>
class STLBST : public Dictionary<K,V> {
public:
    int getSize();
    bool isEmpty();
    void insert(K key, V value);
    void update(K key, V value);
    V get(K key);
    bool contains(K key);
    void remove(K key);
    std::vector<K> getKeys();
    std::vector<pair<K,V> > getItems();
private:
    map<K,V> actualDictionary;
};

template <typename K, typename V>
int STLBST<K,V>::getSize() {
    return this->actualDictionary.size();
}

template <typename K, typename V>
bool STLBST<K,V>::isEmpty() {
    return this->actualDictionary.empty();
}

template <typename K, typename V>
void STLBST<K,V>::insert(K key, V value) {
    if (this->contains(key)) {
        throw runtime_error("STLBST::insert: key already present");
    }
    this->actualDictionary[key] = value;
}

template <typename K, typename V>
void STLBST<K,V>::update(K key, V value) {
    if (!this->contains(key)) {
        throw runtime_error("STLBST::update: key not present");
    }
    this->actualDictionary[key] = value;
}

template <typename K, typename V>
V STLBST<K,V>::get(K key) {
    if (!this->contains(key)) {
        throw runtime_error("STLBST::get: key not present");
    }
    return this->actualDictionary[key];
}

template <typename K, typename V>
bool STLBST<K,V>::contains(K key) {
    try {
        this->actualDictionary.at(key); // check to see if the key exists
        return true;
    } catch (std::out_of_range& e) {
        return false;
    }
}

template <typename K, typename V>
void STLBST<K,V>::remove(K key) {
    if (!this->contains(key)) {
        throw runtime_error("STLBST::remove: key not present");
    }
    this->actualDictionary.erase(key);
}

template <typename K, typename V>
std::vector<K> STLBST<K,V>::getKeys() {
    vector<K> vec;
    for (auto it = actualDictionary.begin();
            it != actualDictionary.end();
            it++) {
        vec.push_back(it->first);
    }
    return vec;
}

template <typename K, typename V>
std::vector<pair<K,V> > STLBST<K,V>::getItems() {
    vector<pair<K,V> > vec;
    for (auto it = actualDictionary.begin();
            it != actualDictionary.end();
            it++) {
        vec.push_back(pair<K,V>(it->first, it->second));
    }
    return vec;
}
