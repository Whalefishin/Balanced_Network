#pragma once

/**
 * The Stack is a pure virtual (abstract) interface of a stack.  Stack
 * implementations should inherit from and implement this interface.
 * @tparam The type of elements stored in this stack.
 */
template <typename T>
class Stack {
public:
    virtual ~Stack() { };

    /**
     * Adds an element to this stack.
     * @param element The element to add to the stack.
     */
    virtual void push(T item) = 0;

    /**
     * Removes the top element from the stack.
     * @return The element that was removed.
     * @throws runtime_error If the stack is empty.
     */
    virtual T pop() = 0;

    /**
     * Determines the size of the stack.
     * @return The number of elements in the stack.
     */
    virtual int getSize() = 0;

    /**
     * Determines whether the stack is empty.
     * @return true if the stack is empty; false otherwise.
     */
    virtual bool isEmpty() = 0;

    /**
     * Retrieves the top element of this stack (without removing it).
     * @return The top element of this stack.
     * @throws runtime_error If the stack is empty.
     */
    virtual T getTop() = 0;

    // You can safely ignore the following code.  This eliminates some default
    // class routines, preventing you from using them accidentally.
    // Specifically, we are disabling the use of the copy constructor and the
    // copy assignment operator.  You can read more here:
    //   http://www.cplusplus.com/articles/y8hv0pDG/
public:
    Stack() { }
private:
    Stack(const Stack& other) = delete;
    Stack& operator=(const Stack& other) = delete;
};

