#pragma once

/**
 * The Queue is a pure virtual (abstract) interface of a queue.  Queue
 * implementations should inherit from and implement this interface.
 * @tparam T The type of elements stored in the queue.
 */
template <typename T>
class Queue {
public:
    virtual ~Queue() { };

    /**
     * Adds an element to the queue.
     * @param item The element to enqueue.
     */
    virtual void enqueue(T element) = 0;

    /**
     * Removes an element from the queue.
     * @return The element that was dequeued.
     * @throws runtime_error If the queue was empty.
     */
    virtual T dequeue() = 0;

    /**
     * Determines the size of this queue.
     * @return The number of elements in the queue.
     */
    virtual int getSize() = 0;

    /**
     * Determines whether this queue is empty.
     * @return true if the queue is empty; false if it is not.
     */
    virtual bool isEmpty() = 0;

    /**
     * Retrieves the front element from the queue (without removing it).
     * @return The next element in the queue.
     * @throws runtime_error If the queue is empty.
     */
    virtual T getFront() = 0;

    // You can safely ignore the following code.  This eliminates some default
    // class routines, preventing you from using them accidentally.
    // Specifically, we are disabling the use of the copy constructor and the
    // copy assignment operator.  You can read more here:
    //   http://www.cplusplus.com/articles/y8hv0pDG/
public:
    Queue() { }
private:
    Queue(const Queue& other) = delete;
    Queue& operator=(const Queue& other) = delete;
};

