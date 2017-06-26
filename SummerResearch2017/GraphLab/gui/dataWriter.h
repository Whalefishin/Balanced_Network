#pragma once

#include <functional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "../adts/dictionary.h"
#include "../adts/edge.h"
#include "../adts/graph.h"
#include "../adts/list.h"

class DataWriter {
public:
    DataWriter(std::ostream& output, bool debug=false);

    void writeInt(int i);
    void writeString(std::string s);

    template<typename T, typename U>
    void writePair(
        std::pair<T,U> pair,
        std::function<void (T)> firstWriter,
        std::function<void (U)> secondWriter);

    template<typename T>
    void writeList(List<T>* elements, std::function<void (T)> elementWriter);

    template<typename T>
    void writeVector(
        std::vector<T> elements, std::function<void (T)> elementWriter);

    template<typename K, typename V>
    void writeDictionary(
        Dictionary<K,V>* dictionary,
        std::function<void (K)> keyWriter,
        std::function<void (V)> valueWriter);

    template<typename V, typename E, typename W>
    void writeEdge(
        Edge<V,E,W> graph,
        std::function<void (V)> vertexWriter,
        std::function<void (E)> edgeDataWriter,
        std::function<void (W)> weightWriter);

    template<typename V, typename E, typename W>
    void writeGraph(
        Graph<V,E,W>* graph,
        std::function<void (V)> vertexWriter,
        std::function<void (E)> edgeDataWriter,
        std::function<void (W)> weightWriter);
protected:
    void writeCharToStream(char c);
    void writeStringToStream(std::string s);
    void flushStream();

private:
    std::ostream& output;
    bool debug;
};

#include "dataWriter-inl.h"

