#pragma once

/**
 * The PriorityQueue is a pure virtual (abstract) interface of a maximum
 * priority queue.
 * @tparam P The type of priorities of queue elements.
 * @tparam V The type of the queue elements themselves.
 */
template <typename P, typename V>
class PriorityQueue {
public:
    virtual ~PriorityQueue() { };

    /**
     * Adds an item with given priority and value to this priority queue.
     * @param priority The priority of the item being added.
     * @param value The value to store in the queue at that priority.
     */
    virtual void insert(P priority, V value) = 0;

    /**
     * Removes the item with maximum priority and returns its value.
     * @return The value of the removed item.
     * @throws runtime_error if there are no items in the priority queue.
     */
    virtual V removeMax() = 0;

    /**
     * Retrieves the item with maximum priority in this queue (without removing
     * it).
     * @return The value of the item with maximum priority.
     * @throws runtime_error If there are no items in the priority queue.
     */
    virtual V getMax() = 0;

    /**
     * Retrieves the maximum priority of items in the queue.
     * @return The maximum priority of any item in the queue.
     * @throws runtime_error if there are no items in the priority queue.
     */
    virtual P getMaxPriority() = 0;

    /**
     * Determines the number of elements in the priority queue.
     * @return The number of elements in the priority queue.
     */
    virtual int getSize() = 0;

    /**
     * Determines whether this priority queue is empty.
     * @return true if the priority
     */
    virtual bool isEmpty() = 0;

    // You can safely ignore the following code.  This eliminates some default
    // class routines, preventing you from using them accidentally.
    // Specifically, we are disabling the use of the copy constructor and the
    // copy assignment operator.  You can read more here:
    //   http://www.cplusplus.com/articles/y8hv0pDG/
public:
    PriorityQueue() { }
private:
    PriorityQueue(const PriorityQueue& other) = delete;
    PriorityQueue& operator=(const PriorityQueue& other) = delete;
};

