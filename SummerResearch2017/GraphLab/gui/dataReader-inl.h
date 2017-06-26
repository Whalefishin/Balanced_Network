#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include "dataReader.h"
#include "formatError.h"

template <typename T, typename U>
std::pair<T,U> DataReader::readPair(std::function<T ()> firstReader, std::function<U ()> secondReader) {
    requireCharFromStream('p', "Expected pair");
    T t = firstReader();
    U u = secondReader();
    requireCharFromStream('$', "Expected $ terminator for pair");
    return std::pair<T,U>(t,u);
}

template <typename T>
std::vector<T> DataReader::readList(std::function<T ()> elementReader) {
    requireCharFromStream('l', "Expected list");
    requireCharFromStream('[', "Expected list start symbol");
    std::vector<T> result;
    char c = readCharFromStream();
    while (c == ',') {
        T t = elementReader();
        result.push_back(t);
        c = readCharFromStream();
    }
    if (c != ']') {
        throw FormatError("Expected list stop symbol");
    }
    return result;
}
