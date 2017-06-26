#pragma once

#include <iterator>
#include <list>
#include <stdexcept>

#include "../list.h"

using std::list;
using std::runtime_error;

template <typename T>
class STLList : public List<T> {
public:
    int getSize();
    bool isEmpty();
    T peekHead();
    T peekTail();
    T get(int i);
    void insertAtHead(T value);
    void insertAtTail(T value);
    T removeHead();
    T removeTail();
private:
    list<T> actualList;
};

template <typename T>
int STLList<T>::getSize() {
    return actualList.size();
}

template <typename T>
bool STLList<T>::isEmpty() {
    return actualList.empty();
}

template <typename T>
T STLList<T>::peekHead() {
    if (actualList.empty()) {
        throw runtime_error("peekHead: empty list");
    } else {
        return actualList.front();
    }
}

template <typename T>
T STLList<T>::peekTail() {
    if (actualList.empty()) {
        throw runtime_error("peekTail: empty list");
    } else {
        return actualList.back();
    }
}

template <typename T>
T STLList<T>::get(int i) {
    if (i<0 || i>=actualList.size()) {
        throw runtime_error("get: invalid index");
    }
    return *(std::next(actualList.begin(), i));
}

template <typename T>
void STLList<T>::insertAtHead(T value) {
    actualList.push_front(value);
}

template <typename T>
void STLList<T>::insertAtTail(T value) {
    actualList.push_back(value);
}

template <typename T>
T STLList<T>::removeHead() {
    if (actualList.empty()) {
        throw runtime_error("removeHead: empty list");
    } else {
        T value = peekHead();
        actualList.pop_front();
        return value;
    }
}

template <typename T>
T STLList<T>::removeTail() {
    if (actualList.empty()) {
        throw runtime_error("removeTail: empty list");
    } else {
        T value = peekTail();
        actualList.pop_back();
        return value;
    }
}


