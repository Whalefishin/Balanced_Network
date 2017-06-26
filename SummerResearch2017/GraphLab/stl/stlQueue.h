#pragma once

#include <queue>
#include <stdexcept>

#include "../adts/queue.h"

using std::runtime_error;
using std::queue;

template <typename T>
class STLQueue : public Queue<T> {
public:
    void enqueue(T item);
    T dequeue();
    int getSize();
    bool isEmpty();
    T getFront();
private:
    queue<T> actualQueue;
};

template <typename T>
void STLQueue<T>::enqueue(T item) {
    actualQueue.push(item);
}

template <typename T>
T STLQueue<T>::dequeue() {
    if (actualQueue.empty()) {
        throw runtime_error("dequeue: empty queue");
    }
    T value = actualQueue.front();
    actualQueue.pop();
    return value;
}

template <typename T>
int STLQueue<T>::getSize() {
    return actualQueue.size();
}

template <typename T>
bool STLQueue<T>::isEmpty() {
    return actualQueue.empty();
}

template <typename T>
T STLQueue<T>::getFront() {
    if (actualQueue.empty()) {
        throw runtime_error("getFront: empty queue");
    }
    return actualQueue.front();
}

