#include <iostream>
#include <stdlib.h>
#include <string>
#include <time.h>

#include "ioUtils.h"
#include "railwayGUI.h"
#include "railwayGame.h"


int main(int argc, char** argv) {
  if (argc != 4) {
      cerr << "Expected two arguments:" << endl;
      cerr << "  1. Filename for graph data" << endl;
      cerr << "  2. Filename for vertex positions" << endl;
      cerr << "  3. Filename for background image" << endl;
      return 1;
  }

  // Initialize randomizer.  This should happen before any random numbers are
  // selected.
  srand(time(NULL));

  // Get command-line arguments.
  string graphDataFilename = argv[1];
  string vertexPositionFilename = argv[2];
  string backgroundImageFilename = argv[3];

  // Retrieve vertex positions (so we know where each vertex goes on the map).
  Dictionary<std::string, std::pair<int,int>>* vertexPositions;
  try {
      vertexPositions = readVertexPositions(vertexPositionFilename);
  } catch (std::exception& e) {
      cout << "Could not read vertex positions file "
           << vertexPositionFilename << ": " << e.what() << endl;
      return 1;
  }

  // Load the initial graph.
  Graph<string,int,int>* graph;
  try {
      graph = readRailwayGraph(graphDataFilename);
  } catch (std::exception& e) {
      cout << "Could not read graph data file " << graphDataFilename
           << ": " << e.what() << endl;
      delete vertexPositions;
      return 1;
  }

  RailwayGUI gui(backgroundImageFilename, vertexPositions);

  // Tell the GUI about the graph we have.
  gui.updateRailwayMap(graph);

  RailwayGame myGame;
  myGame.loadMap(graph);
  myGame.buildGoals();

  //myGame.updateGoals();
  int p1track = myGame.getTracks(1);
  cout<<p1track<<endl;

  vector<Goal*> p1goalVector = myGame.getGoals(1);

  for(int i=0; i<3; i++){
    Goal* tempGoal = p1goalVector[i];
    cout<<tempGoal->location1<<" "<<tempGoal->location2<<endl;
  }

  int test =rand()%10;
  cout<<test<<endl;
  return 0;
}
