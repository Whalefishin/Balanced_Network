#pragma once

#include <stdexcept>
#include <queue>
#include <vector>
#include <utility>

using std::pair;
using std::runtime_error;
using std::vector;

#include "../adts/priorityQueue.h"

template <typename T, typename U>
class FirstLess {
public:
    bool operator()(pair<T,U> a, pair<T,U> b);
};

template <typename T, typename U>
bool FirstLess<T,U>::operator()(pair<T,U> a, pair<T,U> b) {
    return a.first < b.first;
}

template <typename P, typename V>
class STLPriorityQueue : public PriorityQueue<P,V> {
public:
    STLPriorityQueue();
    STLPriorityQueue(vector<pair<P,V> > values);
    ~STLPriorityQueue();
    void insert(P priority, V value);
    V removeMax();
    V getMax();
    P getMaxPriority();
    int  getSize();
    bool isEmpty();
private:
    std::priority_queue<pair<P,V>, vector<pair<P,V> >, FirstLess<P,V> >*
        actualPriorityQueue;
};


template <typename P, typename V>
STLPriorityQueue<P,V>::STLPriorityQueue() {
    actualPriorityQueue =
        new std::priority_queue<
            pair<P,V>,
            vector<pair<P,V> >, FirstLess<P,V> >();
}

template <typename P, typename V>
STLPriorityQueue<P,V>::STLPriorityQueue(vector<pair<P,V> > values) {
    actualPriorityQueue =
        new std::priority_queue<
            pair<P,V>,
            vector<pair<P,V> >, FirstLess<P,V> >(
                FirstLess<P,V>(), values);
}

template <typename P, typename V>
STLPriorityQueue<P,V>::~STLPriorityQueue() {
    delete actualPriorityQueue;
}

template <typename P, typename V>
void STLPriorityQueue<P,V>::insert(P priority, V value) {
    actualPriorityQueue->push(pair<P,V>(priority, value));
}

template <typename P, typename V>
V STLPriorityQueue<P,V>::removeMax() {
    if (actualPriorityQueue->empty()) {
        throw runtime_error("STLPriorityQueue::removeMax(): empty prio queue");
    }
    V v = actualPriorityQueue->top().second;
    actualPriorityQueue->pop();
    return v;
}

template <typename P, typename V>
V STLPriorityQueue<P,V>::getMax() {
    if (actualPriorityQueue->empty()) {
        throw runtime_error("STLPriorityQueue::getMax(): empty prio queue");
    }
    return actualPriorityQueue->top().second;
}

template <typename P, typename V>
P STLPriorityQueue<P,V>::getMaxPriority() {
    if (actualPriorityQueue->empty()) {
        throw runtime_error("STLPriorityQueue::getMaxPriority(): empty prio queue");
    }
    return actualPriorityQueue->top().first;
}

template <typename P, typename V>
int  STLPriorityQueue<P,V>::getSize() {
    return actualPriorityQueue->size();
}

template <typename P, typename V>
bool STLPriorityQueue<P,V>::isEmpty() {
    return actualPriorityQueue->empty();
}
