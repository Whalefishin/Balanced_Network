/*
 * Defines graph searching algorithms
 * 
 * Bill Huang & Genji Kawakita
 */

#include <stdexcept>

#include "stl/stlStack.h"
#include "adjacencyListUndirectedGraph.h"
#include "adjacencyListGraph.h"
#include "stl/stlPriorityQueue.h"
#include "stl/stlQueue.h"



template <typename V, typename E, typename W>
bool reachableDFS(V src, V dest, Graph<V,E,W>* g) {
  if (!g->containsVertex(src)){
    throw runtime_error("Source vertex is not in the graph!");
  }
  if (!g->containsVertex(dest)){
    throw runtime_error("Destination vertex is not in the graph!");
  }

  STLStack<V> stack;
  STLHashTable<V,bool> reachableVertex;

  stack.push(src);

  while(!stack.isEmpty()){
    V current = stack.pop();
    if (!reachableVertex.contains(current)){
    reachableVertex.insert(current,true);
}
    if (current == dest){
      return true;
    }
    else{
      vector<V> neighbors = g->getNeighbors(current);
      for (int j=0;j<neighbors.size();j++){
        if (!reachableVertex.contains(neighbors[j])){
          stack.push(neighbors[j]);
        }
      }
    }
  }
  return false;
}

template <typename V, typename E, typename W>
vector<V> shortestLengthPathBFS(V src, V dest, Graph<V,E,W>* g) {
  if (!g->containsVertex(src)){
    throw runtime_error("Source vertex is not in the graph!");
  }
  if (!g->containsVertex(dest)){
    throw runtime_error("Destination vertex is not in the graph!");
  }

  STLHashTable<V,V> pathTable;

  STLHashTable<V,int> distances;
  vector<V> vertices = g->getVertices();

  for (int i=0;i<vertices.size();i++){
    distances.insert(vertices[i],-1);
  }
  distances.update(src,0);

  STLQueue<V> queue;
  queue.enqueue(src);

  while(!queue.isEmpty()){
    V current = queue.dequeue();
    int current_distance = distances.get(current);
    vector<V> neighbors = g->getNeighbors(current);

    for (int j=0;j<neighbors.size();j++){
      if (distances.get(neighbors[j])==-1){
        distances.update(neighbors[j],current_distance+1);
        queue.enqueue(neighbors[j]);
        pathTable.insert(neighbors[j],current);
      }
    }
  }

  if (distances.get(dest) == -1){
    throw runtime_error("Cannot get from source vertex to destination vertex");
  }


  vector<V> shortestPath_reversed;

  shortestPath_reversed.push_back(dest);
  V currentNode = pathTable.get(dest);

  while (currentNode!=src){
    shortestPath_reversed.push_back(currentNode);
    currentNode = pathTable.get(currentNode);
  }

  shortestPath_reversed.push_back(src);


  vector<V> shortestPath;
  int size = shortestPath_reversed.size();
  for (int k=0;k<size;k++){
    shortestPath.push_back(shortestPath_reversed[size-k-1]);
  }

  return shortestPath;
}

template <typename V, typename E, typename W>
Dictionary<V,W>* singleSourceShortestPath(V src, Graph<V,E,W>* g) {
  if (!g->containsVertex(src)){
    throw runtime_error("The source vertex is not in the graph!");
  }

  vector<V> vertices = g->getVertices();
  Dictionary<V,int>* distances = new STLHashTable<V,int>;

  STLPriorityQueue<int,V> priority_queue;

  //insert the value 1 as infinity. This works because we will
  //insert other weights as negatives, so a default value of 1
  //is greater than any updated value.


  distances->insert(src,0);
  priority_queue.insert(0,src);

  while(!priority_queue.isEmpty()){
    V current = priority_queue.removeMax();
    vector<V> neighbors = g->getNeighbors(current);

    for (int k=0;k<neighbors.size();k++){
      int newCost = -1*(distances->get(current) + g->getEdge(current,neighbors[k]).weight);
      if (distances->contains(neighbors[k])){
        if ((-1)*newCost < distances->get(neighbors[k])){
          distances->update(neighbors[k],(-1)*newCost);
          priority_queue.insert(newCost,neighbors[k]);
        }
      }
      else{
        distances->insert(neighbors[k],(-1)*newCost);
        priority_queue.insert(newCost,neighbors[k]);
      }
    }

  }
  return distances;
}
