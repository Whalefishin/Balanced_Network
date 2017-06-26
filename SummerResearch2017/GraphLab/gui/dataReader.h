#pragma once

#include <functional>
#include <istream>
#include <string>
#include <utility>
#include <vector>

class DataReader {
public:
    DataReader(std::istream& input, bool debug=false);

    /**
     * Reads an integer.  Integers in the stream consist of an 'i' type code,
     * an optional negative sign, a non-empty sequence of digits, and a dollar
     * sign.
     * @return The integer which was read.
     * @throws FormatError If the format of the stream is invalid.
     */
    int readInt();

    /**
     * Reads a string.  Strings in the stream consist of an 's' type code, an
     * integer length, a number of characters equal to that length, and a '$'
     * terminator.
     * @return The string which was read.
     * @throws FormatError If the format of the stream is invalid.
     */
    std::string readString();

    /**
     * Reads a pair.  Pairs in the stream consist of a 'p' type code, the two
     * elements in the pair, and a '$' terminator.
     */
    template <typename T, typename U>
    std::pair<T,U> readPair(std::function<T ()> firstReader, std::function<U ()> secondReader);

    /**
     * Reads a list.  Lists in the stream consist of an 'l' type code, the list
     * start symbol '[', a series of list elements, and a list end symbol ']'.
     * A list element is the symbol ',' followed by the contents of that
     * element.
     */
    template <typename T>
    std::vector<T> readList(std::function<T ()> elementReader);

protected:
    /** Reads a single character from the underlying stream. */
    char readCharFromStream();
    /**
     * Requires that the next character in the stream be the one specified.
     * If this is not the case, an exception is thrown with the given message.
     * @param c The character to expect.
     * @param message The message to produce.
     * @throws FormatError If the format of the stream is invalid.
     */
    void requireCharFromStream(char c, std::string message);

private:
    std::istream& input;
    bool debug;
};

#include "dataReader-inl.h"

