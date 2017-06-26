#pragma once

#include <stack>
#include <stdexcept>

#include "../adts/stack.h"

using std::runtime_error;
using std::stack;

template <typename T>
class STLStack : public Stack<T> {
public:
    void push(T item);
    T pop();
    int getSize();
    bool isEmpty();
    T getTop();
private:
    stack<T> actualStack;
};

template <typename T>
void STLStack<T>::push(T item) {
    actualStack.push(item);
}

template <typename T>
T STLStack<T>::pop() {
    if (actualStack.empty()) {
        throw runtime_error("pop: empty stack");
    }
    T value = actualStack.top();
    actualStack.pop();
    return value;
}

template <typename T>
int STLStack<T>::getSize() {
    return actualStack.size();
}

template <typename T>
bool STLStack<T>::isEmpty() {
    return actualStack.empty();
}

template <typename T>
T STLStack<T>::getTop() {
    if (actualStack.empty()) {
        throw runtime_error("getTop: empty stack");
    }
    return actualStack.top();
}

