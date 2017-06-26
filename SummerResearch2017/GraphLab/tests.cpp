#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include <UnitTest++/UnitTest++.h>

#include "adts/dictionary.h"
#include "adts/graph.h"
#include "adjacencyListGraph.h"
#include "adjacencyListUndirectedGraph.h"
#include "graphAlgorithms.h"
#include "ioUtils.h"


using namespace std;

template <typename V, typename E, typename W, template <class,class,class> class G = AdjacencyListUndirectedGraph>
Graph<V,E,W>* makeGraph(vector<V> vertices, vector<Edge<V,E,W>> edges) {
    Graph<V,E,W>* graph = new G<V,E,W>();
    for (int i=0;i<vertices.size();i++) {
        graph->insertVertex(vertices.at(i));
    }
    for (int i=0;i<edges.size();i++) {
        Edge<V,E,W> edge = edges.at(i);
        graph->insertEdge(edge.source, edge.target, edge.label, edge.weight);
    }
    return graph;
}

TEST(dfsTwoVertexGraph) {
    // This creates the *undirected* graph:
    //    "1"  --- 1,true ---  "2"
    Graph<string,bool,int>* graph = makeGraph<string,bool,int>(
        { "1", "2" }, { { "1", true, 1, "2" } });
    CHECK(reachableDFS(string("1"),string("2"),graph));
    CHECK(reachableDFS(string("2"),string("1"),graph));
    delete graph;
}

TEST(dfsTwoDisconnectedVertices) {
    // This creates the *undirected* graph:
    //    "1"                  "2"
    Graph<string,bool,int>* graph = makeGraph<string,bool,int>(
        { "1", "2" }, { });
    CHECK(!reachableDFS(string("1"),string("2"),graph));
    CHECK(!reachableDFS(string("2"),string("1"),graph));
    delete graph;
}

TEST(dfsThreeVertexDirectedGraph) {
    // This call makes a *directed* graph!
    Graph<string,bool,int>* graph = makeGraph<string,bool,int,AdjacencyListGraph>(
        { "1", "2", "3" },
        { { "1", true, 1, "2" },
          { "1", true, 1, "3" },
        });
    CHECK(reachableDFS(string("1"),string("1"),graph));
    CHECK(reachableDFS(string("1"),string("2"),graph));
    CHECK(reachableDFS(string("1"),string("3"),graph));
    CHECK(!reachableDFS(string("3"),string("1"),graph));
    CHECK(!reachableDFS(string("2"),string("1"),graph));
    CHECK(!reachableDFS(string("2"),string("3"),graph));
    delete graph;
}

TEST(dfsSocialNetwork25) {
    Graph<string,int,int>* graph =
        readRailwayGraph("test_data/socialNetwork_25.graph");
    vector<string> vertices = graph->getVertices();
    CHECK(reachableDFS(vertices[3],vertices[17],graph));
    CHECK(reachableDFS(vertices[18],vertices[9],graph));
    CHECK(reachableDFS(vertices[22],vertices[12],graph));
    CHECK(reachableDFS(vertices[1],vertices[18],graph));
    delete graph;
}

TEST(dfsSocialNetwork1000) {
    Graph<string,int,int>* graph =
        readRailwayGraph("test_data/socialNetwork_1000.graph");
    vector<string> vertices = graph->getVertices();
    CHECK(reachableDFS(vertices[34],vertices[71],graph));
    CHECK(reachableDFS(vertices[18],vertices[92],graph));
    CHECK(reachableDFS(vertices[46],vertices[15],graph));
    CHECK(reachableDFS(vertices[8],vertices[66],graph));
    delete graph;
}

TEST(dfsSocialNetworkDisconnected20) {
    Graph<string,int,int>* graph =
        readRailwayGraph("test_data/socialNetworkDisconnected_20.graph");
    vector<string> vertices = graph->getVertices();
    CHECK(reachableDFS(string("Kevin Bacon"), string("Pigby Strangepork"),graph));
    CHECK(reachableDFS(string("Pigby Strangepork"), string("Nanny Summer"),graph));
    CHECK(reachableDFS(string("Hogzilla Bland"), string("Piglet Corpus"),graph));
    CHECK(reachableDFS(string("Hamlet Pigasus"), string("Toot Bacon"),graph));
    CHECK(!reachableDFS(string("Kevin Bacon"), string("Toot Bacon"),graph));
    CHECK(!reachableDFS(string("Pigby Strangepork"), string("Hogzilla Bland"),graph));
    CHECK(!reachableDFS(string("Nanny Summer"), string("Piglet Corpus"),graph));
    CHECK(!reachableDFS(string("Hamlet Pigasus"), string("Pigby Strangepork"),graph));
    delete graph;
}

#define CHECK_EQUAL_PATH(expectedExpr,actualExpr) \
    { \
        vector<string> expected = expectedExpr; \
        vector<string> actual = actualExpr; \
        string expectedPath = ""; \
        for (int i=0;i<expected.size();i++) { \
            if (expectedPath.size() > 0) { \
                expectedPath += " -- "; \
            } \
            expectedPath += expected.at(i); \
        } \
        string actualPath = ""; \
        for (int i=0;i<actual.size();i++) { \
            if (actualPath.size() > 0) { \
                actualPath += " -- "; \
            } \
            actualPath += actual.at(i); \
        } \
        CHECK_EQUAL( \
            "path = " + expectedPath, \
            "path = " + actualPath); \
    }

TEST(bfsThreeVertexClique) {
    // This call makes an undirected graph.
    Graph<string,bool,int>* graph = makeGraph<string,bool,int>(
        { "1", "2", "3" },
        { { "1", true, 1, "2" },
          { "1", true, 1, "3" },
          { "2", true, 1, "3" },
        });
    vector<string> path = {"1","3"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(string("1"),string("3"),graph));
    path = {"1","2"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(string("1"),string("2"),graph));
    delete graph;
}

TEST(bfsFiveVertex) {
    // This call makes an undirected graph.
    Graph<string,bool,int>* graph = makeGraph<string,bool,int>(
        { "start", "1", "2", "finish", "alt" },
        { { "start", true, 1, "1" },
          { "1", true, 1, "2" },
          { "2", true, 1, "finish" },
          { "start", true, 1, "alt" },
          { "alt", true, 1, "finish" }
        });
    vector<string> path = {"start", "alt", "finish"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(string("start"),string("finish"),graph));
    delete graph;
}

TEST(bfsFiveVertexDirected) {
    // This call makes a directed graph.
    Graph<string,bool,int>* graph = makeGraph<string,bool,int,AdjacencyListGraph>(
        { "start", "1", "2", "finish", "alt" },
        { { "start", true, 1, "1" },
          { "1", true, 1, "2" },
          { "2", true, 1, "finish" },
          { "finish", true, 1, "alt" },
          { "alt", true, 1, "start" }
        });
    vector<string> path = {"start", "1", "2", "finish"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(string("start"),string("finish"),graph));
    delete graph;
}

TEST(bfsSocialNetwork25) {
    Graph<string,int,int>* graph = readRailwayGraph("test_data/socialNetwork_25.graph");
    vector<string> path = {"Jaimie Dengar", "Kevin Bacon"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(
        string("Jaimie Dengar"), string("Kevin Bacon"), graph));
    path = {"Sam Curlytail", "Ace Watson", "Kevin Bacon"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(
        string("Sam Curlytail"), string("Kevin Bacon"), graph));
    delete graph;
}

TEST(bfsSocialNetwork1000) {
    Graph<string,int,int>* graph = readRailwayGraph("test_data/socialNetwork_1000.graph");
    vector<string> path = {"Jaimie Dengar", "Kevin Bacon"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(
        string("Jaimie Dengar"), string("Kevin Bacon"), graph));
    path = {"Sam Curlytail", "Waddles Hog", "Kevin Bacon"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(
        string("Sam Curlytail"), string("Kevin Bacon"), graph));
    path = {"Napoleon Boar", "Huxley Poppleton", "Kevin Bacon"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(
        string("Napoleon Boar"), string("Kevin Bacon"), graph));
    path = {"Lester Strangepork", "Kevin Bacon"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(
        string("Lester Strangepork"), string("Kevin Bacon"), graph));
    path = {"Annie Sue Saddlebottom", "Truffles Floyd", "Kevin Bacon"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(
        string("Annie Sue Saddlebottom"), string("Kevin Bacon"), graph));
    delete graph;
}

TEST(bfsSocialNetworkDisconnected20) {
    Graph<string,int,int>* graph = readRailwayGraph("test_data/socialNetworkDisconnected_20.graph");
    vector<string> path = {"Piglet Corpus", "Rooter Ziffel", "Hogzilla Bland"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(
        string("Piglet Corpus"), string("Hogzilla Bland"), graph));
    path = {"Freddy Wen", "Kevin Bacon"};
    CHECK_EQUAL_PATH(path, shortestLengthPathBFS(
        string("Freddy Wen"), string("Kevin Bacon"), graph));
    CHECK_THROW(shortestLengthPathBFS(
        string("Edessa Pork"), string("Kevin Bacon"), graph),
        std::exception);
    delete graph;
}

TEST(ssspFiveVertex) {
    // This call makes an undirected graph.
    Graph<string,bool,int>* graph = makeGraph<string,bool,int>(
        { "start", "1", "2", "finish", "alt" },
        { { "start", false, 3, "1" },
          { "1", false, 4, "2" },
          { "2", false, 2, "finish" },
          { "start", false, 5, "alt" },
          { "alt", false, 6, "finish" }
        });
    Dictionary<string,int>* results =
        singleSourceShortestPath(string("start"), graph);
    CHECK_EQUAL(9, results->get("finish"));
    CHECK_EQUAL(7, results->get("2"));
    CHECK_EQUAL(5, results->get("alt"));
    delete results;
    delete graph;
}

TEST(ssspFiveVertexDirected) {
    // This call makes a directed graph.
    Graph<string,bool,int>* graph = makeGraph<string,bool,int,AdjacencyListGraph>(
        { "start", "1", "2", "finish", "alt" },
        { { "finish", false, 3, "2" },
          { "2", false, 4, "1" },
          { "1", false, 2, "start" },
          { "start", false, 5, "alt" },
          { "alt", false, 6, "finish" }
        });
    Dictionary<string,int>* results =
        singleSourceShortestPath(string("start"), graph);
    CHECK_EQUAL(11, results->get("finish"));
    CHECK_EQUAL(14, results->get("2"));
    CHECK_EQUAL(5, results->get("alt"));
    delete results;
    delete graph;
}

TEST(ssspTown8) {
    // Creates an undirected graph
    Graph<string,int,int>* graph = readRailwayGraph("test_data/town_8.graph");
    Dictionary<string,int>* results =
        singleSourceShortestPath(string("Pigma's Jewelry Museum"), graph);
    CHECK_EQUAL(15, results->get("Arnold's Convenience Stand"));
    CHECK_EQUAL(1, results->get("Petunia's Bowling Town"));
    CHECK_EQUAL(16, results->get("Sam's Photo Retailer"));
    delete results;
    delete graph;
}

TEST(ssspTown100) {
    // Creates an undirected graph
    Graph<string,int,int>* graph = readRailwayGraph("test_data/town_100.graph");
    Dictionary<string,int>* results =
        singleSourceShortestPath(string("Cochon's Movie Vault"), graph);
    CHECK_EQUAL(24, results->get("Hamhock's Trinket Ltd"));
    CHECK_EQUAL(22, results->get("Toby's Furniture Warehouse"));
    CHECK_EQUAL(22, results->get("Saddlebottom's Photo Trough"));
    delete results;
    results = singleSourceShortestPath(string("Babe's Produce Store"), graph);
    CHECK_EQUAL(23, results->get("Arnold's Convenience Stand"));
    CHECK_EQUAL(35, results->get("Takako's Antique Vault"));
    CHECK_EQUAL(25, results->get("Major's Spice Ltd"));
    delete results;
    delete graph;
}

int main() {
    return UnitTest::RunAllTests();
}
