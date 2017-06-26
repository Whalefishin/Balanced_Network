/*
 * Defines RailwayGame class and its methods
 *
 * Bill Huang & Genji Kawakita
 *
 */


#include <string>
#include "stdlib.h"
#include "railwayGame.h"
#include "goal.h"
#include "stl/stlHashTable.h"
#include "adts/edge.h"


using namespace std;

RailwayGame::RailwayGame(){
  turn_parity = 0;
  p1Score = 0;
  p2Score = 0;
  p1Tracks = 35;
  p2Tracks = 50;
  winner_number = 0;
}

RailwayGame::~RailwayGame(){
  delete p1Map;
  delete p2Map;
  delete p1ClaimableEdges;
  delete p2ClaimableEdges;
  for (int i=0;i<p1Goals.size();i++){
    delete p1Goals[i];
    delete p2Goals[i];
  }
}

void RailwayGame::buildGoals(){

  vector<Goal*> goalVector;
  vector<string> goalVertices;
  STLHashTable<int,int>  indexTable;
  //get the vertices in the graph
  vector<string> verticesVector = railwayMap->getVertices();
  int vSize = verticesVector.size();
  //randomlly picks 12 vertices
  while(indexTable.getSize()<12){
    int randomIndex = rand()%vSize;
    if(!indexTable.contains(randomIndex)){
      indexTable.insert(randomIndex, randomIndex);
      goalVertices.push_back(verticesVector[randomIndex]);
    }
  }
  for(int i=0; i<12; i=i+2){
    vector<string> distanceVector = shortestLengthPathBFS(goalVertices[i], goalVertices[i+1], railwayMap);
    // generates random number till we get the vertex that is more than 3 routes away from goalVertices[i]

    while(distanceVector.size() < 4){
      int randomIndex = rand()%vSize;
      if(!indexTable.contains(randomIndex)){
        indexTable.insert(randomIndex, randomIndex);
        goalVertices[i+1] = verticesVector[randomIndex];
        distanceVector = shortestLengthPathBFS(goalVertices[i], goalVertices[i+1], railwayMap);
      }
    }
    Dictionary<string,int>* pointsDictionary = singleSourceShortestPath(goalVertices[i], railwayMap);
    int tempPoints = pointsDictionary->get(goalVertices[i+1]);
    int points = tempPoints/4;
    string message = to_string(points) + " points";
    Goal* g = new Goal(goalVertices[i],goalVertices[i+1],message,points);
    goalVector.push_back(g);
    delete pointsDictionary;
  }

  // each player gets three goals from goalVector
  for(int i=0;i<3;i++){
    p1Goals.push_back(goalVector[i]);
  }
  for(int i=0;i<3;i++){
    p2Goals.push_back(goalVector[i+3]);
  }
}


void RailwayGame::updateGoals() {
  for(int i=0; i<3; i++){
    if(!reachableDFS(p1Goals[i]->location1,p1Goals[i]->location2,p1Map)){
      p1Goals[i]->message = "impossible";
    }
    /**if the path is completed, we need to update p1's score
    if p1Cities contains both location1 and 2 of each goal, the two location must be somehow connected
    */
    if(p1Cities.contains(p1Goals[i]->location1) && p1Cities.contains(p1Goals[i]->location2)){
      p1Score = p1Score + p1Goals[i]->points;
      p1Goals[i]->message = "completed";
    }
  }
  for(int i=0; i<3; i++){
    if(!reachableDFS(p2Goals[i]->location1,p2Goals[i]->location2,p2Map)){
      p2Goals[i]->message = "impossible";
    }
    //if the path is comleted, we need to update p2's score
    if(p2Cities.contains(p2Goals[i]->location1) && p2Cities.contains(p2Goals[i]->location2)){
      p2Score = p2Score + p2Goals[i]->points;
      p2Goals[i]->message = "completed";
    }
  }
}

void RailwayGame::claimEdge(string src, string dest){

  Edge<string,int,int> edge = railwayMap->getEdge(src,dest);

  if (edge.label !=0){
    bottomMessage = "Player: You can't claim that route because it's already taken";
    return;
  }
  //check whose turn it is, and claim the edge for that player
  //if it's possible.


  //if it's player1's turn
  if (trackTurn()==1){
    //if p1 does not own any edges yet, p1 can claim any edge he wants
    //as long as he can afford it.
    if (p1Cities.getSize()==0){
      p1Cities.insert(src,true);
      p1Cities.insert(dest,true);

      p1Tracks = p1Tracks - edge.weight;

      // +2 points for the score because p1 has two more cities
      p1Score = p1Score + 2;

      p2Map->removeEdge(src,dest);
      if (p2ClaimableEdges->containsEdge(src,dest)){
      p2ClaimableEdges->removeEdge(src,dest);
    }

      railwayMap->removeEdge(src,dest);
      railwayMap->insertEdge(src,dest,1,edge.weight);


      vector<Edge<string,int,int>> edgeList1 = p1Map->getOutgoingEdges(src);
      vector<Edge<string,int,int>> edgeList2 = p1Map->getOutgoingEdges(dest);

      //cout<< edgeList1.size() << endl;

      for (int i=0;i<edgeList1.size();i++){
        if (!p1ClaimableEdges->containsEdge(edgeList1[i].source,edgeList1[i].target) && edgeList1[i].label == 0){
          p1ClaimableEdges->insertEdge(edgeList1[i].source,edgeList1[i].target,edgeList1[i].label,edgeList1[i].weight);
        }
      }
      for (int i=0;i<edgeList2.size();i++){
        if (!p1ClaimableEdges->containsEdge(edgeList2[i].source,edgeList2[i].target) && edgeList2[i].label == 0){
          p1ClaimableEdges->insertEdge(edgeList2[i].source,edgeList2[i].target,edgeList2[i].label,edgeList2[i].weight);
        }
      }

      p1ClaimableEdges->removeEdge(src,dest);


      turn_parity++;
      bottomMessage = "Player 2: it is your turn";
      //At this point, p1 has claimed the cities associated with the edge
      //And we also know all the connected edges for p1
    }
    //if p1 has some cities already, then we need to first check if
    //the edge he is trying to claim is in the set of edges that he can
    //claim

    else {
      if (p1ClaimableEdges->containsEdge(edge.source,edge.target)){
        //also check if p1 has enough tracks to pay for the edge
        if (p1Tracks >= edge.weight){
          //after satisfying these conditions, p1 can finally claim
          //the edge as his.
          p1Tracks = p1Tracks - edge.weight;
          //if p1 has claimed an edge, that edge can be seen as
          //being removed in p2's map.
          p2Map->removeEdge(src,dest);
          if (p2ClaimableEdges->containsEdge(src,dest)){
          p2ClaimableEdges->removeEdge(src,dest);
        }

          Edge<string,int,int> newEdge(src,1,edge.weight,dest);

          railwayMap->removeEdge(src,dest);
          railwayMap->insertEdge(src,dest,1,edge.weight);

          //we need to check whether the source or the destination
          //is the new city claimed by p1
          if (!p1Cities.contains(src)){
            p1Score = p1Score+1;
            p1Cities.insert(src,true);
            vector<Edge<string,int,int>> edgeList = p1Map->getOutgoingEdges(src);
            //then we need to update p1's claimable edges.
            for (int i=0;i<edgeList.size();i++){
              if (!p1ClaimableEdges->containsEdge(edgeList[i].source,edgeList[i].target) && edgeList[i].label == 0){
                p1ClaimableEdges->insertEdge(edgeList[i].source,edgeList[i].target,edgeList[i].label,edgeList[i].weight);
              }
            }
          }
          else if (!p1Cities.contains(dest)){
            p1Score = p1Score+1;
            p1Cities.insert(dest,true);
            vector<Edge<string,int,int>> edgeList = p1Map->getOutgoingEdges(dest);
            //then we need to update p1's claimable edges.
            for (int i=0;i<edgeList.size();i++){
              if (!p1ClaimableEdges->containsEdge(edgeList[i].source,edgeList[i].target) && edgeList[i].label == 0){
                p1ClaimableEdges->insertEdge(edgeList[i].source,edgeList[i].target,edgeList[i].label,edgeList[i].weight);
              }
            }
          }
          turn_parity++;
          p1ClaimableEdges->removeEdge(src,dest);

          bottomMessage = "Player 2: it is your turn";
        }
        else {
          bottomMessage = "Player 1: You don't have enough tracks";
        }
      }
      else {
        bottomMessage = "Player 1: You can't claim that route";
      }

    }

  }



  //if it's player2's turn, do exactly the same as we did for p1.
  else if (trackTurn()==2){

    if (p2Cities.getSize()==0){
      p2Cities.insert(src,true);
      p2Cities.insert(dest,true);

      p2Tracks = p2Tracks - edge.weight;

      p2Score = p2Score + 2;

      p1Map->removeEdge(src,dest);
      if (p1ClaimableEdges->containsEdge(src,dest)){
      p1ClaimableEdges->removeEdge(src,dest);
    }

      railwayMap->removeEdge(src,dest);
      railwayMap->insertEdge(src,dest,2,edge.weight);


      vector<Edge<string,int,int>> edgeList1 = p2Map->getOutgoingEdges(src);
      vector<Edge<string,int,int>> edgeList2 = p2Map->getOutgoingEdges(dest);


      for (int i=0;i<edgeList1.size();i++){
        if (!p2ClaimableEdges->containsEdge(edgeList1[i].source,edgeList1[i].target) && edgeList1[i].label == 0){
          p2ClaimableEdges->insertEdge(edgeList1[i].source,edgeList1[i].target,edgeList1[i].label,edgeList1[i].weight);
        }
      }
      for (int i=0;i<edgeList2.size();i++){
        if (!p2ClaimableEdges->containsEdge(edgeList2[i].source,edgeList2[i].target) && edgeList2[i].label == 0){
          p2ClaimableEdges->insertEdge(edgeList2[i].source,edgeList2[i].target,edgeList2[i].label,edgeList2[i].weight);
        }
      }

      p2ClaimableEdges->removeEdge(src,dest);

      turn_parity++;
      bottomMessage = "Player 1: it is your turn";

    }


    else {
      if (p2ClaimableEdges->containsEdge(edge.source,edge.target)){
        if (p2Tracks >= edge.weight){

          p2Tracks = p2Tracks - edge.weight;

          p1Map->removeEdge(src,dest);
          if (p1ClaimableEdges->containsEdge(src,dest)){
          p1ClaimableEdges->removeEdge(src,dest);
        }
          Edge<string,int,int> newEdge(src,2,edge.weight,dest);

          railwayMap->removeEdge(src,dest);
          railwayMap->insertEdge(src,dest,2,edge.weight);

          if (!p2Cities.contains(src)){
            p2Score = p2Score+1;
            p2Cities.insert(src,true);
            vector<Edge<string,int,int>> edgeList = p2Map->getOutgoingEdges(src);

            for (int i=0;i<edgeList.size();i++){
              if (!p2ClaimableEdges->containsEdge(edgeList[i].source,edgeList[i].target) && edgeList[i].label == 0){
                p2ClaimableEdges->insertEdge(edgeList[i].source,edgeList[i].target,edgeList[i].label,edgeList[i].weight);
              }
            }
          }
          else if (!p2Cities.contains(dest)){
            p2Score = p2Score+1;
            p2Cities.insert(dest,true);
            vector<Edge<string,int,int>> edgeList = p2Map->getOutgoingEdges(dest);

            for (int i=0;i<edgeList.size();i++){
              if (!p2ClaimableEdges->containsEdge(edgeList[i].source,edgeList[i].target) && edgeList[i].label == 0){
                p2ClaimableEdges->insertEdge(edgeList[i].source,edgeList[i].target,edgeList[i].label,edgeList[i].weight);
              }
            }
          }
          turn_parity++;

          p2ClaimableEdges->removeEdge(src,dest);



          bottomMessage = "Player 1: it is your turn";
        }
        else {
          bottomMessage = "Player 2: You don't have enough tracks";
        }
      }
      else {
        bottomMessage = "Player 2: You can't claim that route";
      }

    }
  }

}



int RailwayGame::trackTurn(){
  return (turn_parity%2)+1;
}

void RailwayGame::loadMap(Graph<string,int,int>* g){
  railwayMap = g;

  p1Map = new AdjacencyListUndirectedGraph<string,int,int>();
  p2Map = new AdjacencyListUndirectedGraph<string,int,int>();

  p1ClaimableEdges = new AdjacencyListUndirectedGraph<string,int,int>;
  p2ClaimableEdges = new AdjacencyListUndirectedGraph<string,int,int>;

  vector<string> vertices = g->getVertices();
  vector<Edge<string,int,int>> edges = g->getEdges();
  //copy all the vertices into p1&p2 maps as well as the graphs for claimable
  //edges
  for (int i=0;i<vertices.size();i++){
    p1Map->insertVertex(vertices[i]);
    p2Map->insertVertex(vertices[i]);
    p1ClaimableEdges->insertVertex(vertices[i]);
    p2ClaimableEdges->insertVertex(vertices[i]);
  }
  //copy all the edges into p1&p2 maps
  for (int i=0;i<edges.size();i++){
    Edge<string,int,int> edge = edges[i];
    p1Map->insertEdge(edge.source,edge.target,edge.label,edge.weight);
    p2Map->insertEdge(edge.source,edge.target,edge.label,edge.weight);
  }
}

int RailwayGame::getScore(int playerNumber){
  if (playerNumber == 1){
    return p1Score;
  }
  else if (playerNumber ==2 ){
    return p2Score;
  }
  else {
    throw runtime_error("Invalid player number");
  }
}

int RailwayGame::getTracks(int playerNumber){
  if (playerNumber == 1){
    return p1Tracks;
  }
  else if (playerNumber ==2 ){
    return p2Tracks;
  }
  else {
    throw runtime_error("Invalid player number");
  }
}

vector<Goal*> RailwayGame::getGoals(int playerNumber){
  if (playerNumber == 1){
    return p1Goals;
  }
  else if (playerNumber ==2 ){
    return p2Goals;
  }
  else {
    throw runtime_error("Invalid player number");
  }
}

Graph<string,int,int>* RailwayGame::getMap(){
  return railwayMap;
}

string RailwayGame::getMessage(){
  return bottomMessage;
}

void RailwayGame::pass(int playerNumber){
  if (playerNumber == 1){
    p1Tracks = p1Tracks + p1Cities.getSize();
    turn_parity++;
    bottomMessage = "Player 2: it is your turn";
  }
  else if (playerNumber == 2){
    p2Tracks = p2Tracks + p2Cities.getSize();
    turn_parity++;
    bottomMessage = "Player 1: it is your turn";
  }
  else {
    throw runtime_error("Invalid player number");
  }
}

int RailwayGame::getWinner(){
  return winner_number;
}


bool RailwayGame::gameOver(){
  vector<Edge<string,int,int>> edges = railwayMap->getEdges();

  //the game is over if and only if all the edges are claimed.

  for (int i=0;i<edges.size();i++){
    //as long as one edge is unclaimed, the game is not over
    if (edges[i].label ==0){
      return false;
    }
  }

  //if we're here, it means that all the edges are claimed, and we can
  //check to see who is the winner.

    //compare the scores of the two players and see who wins
    if (p1Score > p2Score){
      winner_number = 1;
    }
    else if (p1Score < p2Score){
      winner_number = 2;
    }

    return true;
}
