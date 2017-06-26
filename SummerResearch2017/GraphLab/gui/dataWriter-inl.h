#include <functional>
#include <ostream>
#include <string>
#include <vector>
#include <utility>

#include "../adts/dictionary.h"
#include "../adts/edge.h"
#include "../adts/graph.h"
#include "../adts/list.h"
#include "dataWriter.h"

template<typename T, typename U>
void DataWriter::writePair(std::pair<T,U> pair,
        std::function<void (T)> firstWriter,
        std::function<void (U)> secondWriter) {
    writeCharToStream('p');
    firstWriter(pair.first);
    secondWriter(pair.second);
    writeCharToStream('$');
}

template<typename T>
void DataWriter::writeList(List<T>* elements,
        std::function<void (T)> elementWriter) {
    writeCharToStream('l');
    writeCharToStream('[');
    for (int i=0;i<elements->getSize();i++) {
        writeCharToStream(',');
        elementWriter(elements->get(i));
    }
    writeCharToStream(']');
}

template<typename T>
void DataWriter::writeVector(vector<T> elements,
        std::function<void (T)> elementWriter) {
    writeCharToStream('l');
    writeCharToStream('[');
    for (int i=0;i<elements.size();i++) {
        writeCharToStream(',');
        elementWriter(elements[i]);
    }
    writeCharToStream(']');
}

template<typename K, typename V>
void DataWriter::writeDictionary(Dictionary<K,V>* dictionary,
        std::function<void (K)> keyWriter,
        std::function<void (V)> valueWriter) {
    writeCharToStream('d');
    writeCharToStream('{');
    vector<pair<K,V>> items = dictionary->getItems();
    writeVector<pair<K,V>>(items,
        [&](pair<K,V> pair) {
            return writePair(pair,keyWriter,valueWriter);
        } );
    writeCharToStream('}');
}

template<typename V, typename E, typename W>
void DataWriter::writeEdge(Edge<V,E,W> edge,
        std::function<void (V)> vertexWriter,
        std::function<void (E)> edgeDataWriter,
        std::function<void (W)> weightWriter) {
    writeCharToStream('e');
    vertexWriter(edge.source);
    edgeDataWriter(edge.label);
    weightWriter(edge.weight);
    vertexWriter(edge.target);
    writeCharToStream('$');
}

template<typename V, typename E, typename W>
void DataWriter::writeGraph(Graph<V,E,W>* graph,
        std::function<void (V)> vertexWriter,
        std::function<void (E)> edgeDataWriter,
        std::function<void (W)> weightWriter) {
    writeCharToStream('g');
    writeCharToStream('{');

    vector<V> vertices = graph->getVertices();
    writeVector<V>(vertices, vertexWriter);

    vector<Edge<V,E,W>> edges = graph->getEdges();
    writeVector<Edge<V,E,W>>(edges,
        [&](Edge<V,E,W> edge) {
            return writeEdge(edge, vertexWriter, edgeDataWriter, weightWriter);
        });

    writeCharToStream('}');
}
