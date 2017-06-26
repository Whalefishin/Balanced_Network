/*
 * Defines AdjacencyUndirectedGraph class and its methods
 *
 * Bill Huang & Genji Kawakita
 *
 */


#pragma once

#include "adjacencyListGraph.h"

/**
 * AdjacencyListUndirectedGraph is a class representing
 *   undirected, weighted graphs.
 *
 * It is a subclass of AdjacencyListGraph for the purpose of code reuse.
 * @tparam V a type to represent the vertex labels
 * @tparam E a type to represent the edge labels
 * @tparam W a type to represent the weight (usually an int, float, or double)
 */
template <typename V, typename E, typename W>
class AdjacencyListUndirectedGraph : public AdjacencyListGraph<V,E,W> {
public:
    void insertEdge(V src, V dest, E label, W weight);
    void removeEdge(V src, V dest);
    vector<Edge<V,E,W>> getEdges();
};



/**
 * Adds an undirected edge to this graph.
 * @param source The source vertex for the edge.
 * @param destination The destination vertex for the edge.
 * @param label The label of the edge to add.
 * @param weight The weight of the edge to add.
 * @throws runtime_error If an edge already exists between the given source
 *                       and target.
 */
template <typename V, typename E, typename W>
void AdjacencyListUndirectedGraph<V,E,W>::insertEdge(V src, V dest, E label, W weight) {
  // insert the edge as two directional edges.
  // This lets getNeighbors to work without changes
  AdjacencyListGraph<V,E,W>::insertEdge(src, dest, label, weight);
  AdjacencyListGraph<V,E,W>::insertEdge(dest, src, label, weight);
}

/**
 * Removes an edge from this graph.
 * @param source The source vertex for the edge.
 * @param destination The destination vertex for the edge.
 * @throws runtime_error If the edge does not exist.
 */
template <typename V, typename E, typename W>
void AdjacencyListUndirectedGraph<V,E,W>::removeEdge(V src, V dest) {
  AdjacencyListGraph<V,E,W>::removeEdge(src,dest);
  AdjacencyListGraph<V,E,W>::removeEdge(dest,src);
}

/**
 * Retrieves all edges from this graph.
 * @return A vector containing every edge in this graph.
 */
template <typename V, typename E, typename W>
vector<Edge<V,E,W>> AdjacencyListUndirectedGraph<V,E,W>::getEdges() {
  vector<Edge<V,E,W>> results = AdjacencyListGraph<V,E,W>::getEdges();

  for (int i=0;i<results.size();i++){
    for (int j=i+1;j<results.size();j++){
      if (results[i].source == results[j].target && results[i].target == results[j].source){
        results.erase(results.begin()+j,results.begin()+j+1);
      }
    }
  }
  return results;
}
