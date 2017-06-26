/*
 * Main function that deals with user inputs
 *
 * Bill HUnag and Genji Kawakita
 *
 */

#include <iostream>
#include <stdlib.h>
#include <string>
#include <time.h>

#include "ioUtils.h"
#include "railwayGUI.h"
#include "railwayGame.h"

using namespace std;

int main(int argc, char** argv) {
    // Check command line arguments and give up if necessary.
    if (argc != 4) {
        cerr << "Expected three arguments:" << endl;
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

    // Create the GUI object here.  This is statically allocated, so the GUI
    // will disappear the moment your program laves the main function.
    RailwayGUI gui(backgroundImageFilename, vertexPositions);

    // Tell the GUI about the graph we have.
    gui.updateRailwayMap(graph);


    RailwayGame myGame;
    myGame.loadMap(graph);

    myGame.buildGoals();


    gui.updatePlayerStatus(1,myGame.getScore(1),myGame.getTracks(1),myGame.getGoals(1));
    gui.updatePlayerStatus(2,myGame.getScore(2),myGame.getTracks(2),myGame.getGoals(2));


    gui.updateMessage("Player 1: it is your turn");


    vector<string> decision = gui.getNextMove();
    while(decision[0] != "close"){

      if(decision[0] == "edge"){
        int playerNumber = myGame.trackTurn();
        string src = decision[1];
        string dest = decision[2];

        myGame.claimEdge(src,dest);
        myGame.updateGoals();

        gui.updateRailwayMap(myGame.getMap());
        gui.updatePlayerStatus(1,myGame.getScore(1),myGame.getTracks(1),myGame.getGoals(1));
        gui.updatePlayerStatus(2,myGame.getScore(2),myGame.getTracks(2),myGame.getGoals(2));
        gui.updateMessage(myGame.getMessage());

        //if the game is over, update the message to reflect the player who won.
        if (myGame.gameOver()){
          if (myGame.getWinner() == 0){
            gui.updateMessage("Wow, you guys have a draw, incredible!");
          }
          else{
          gui.updateMessage("Game is over. The winner is Player" + to_string(myGame.getWinner()));
        }
      }

        decision = gui.getNextMove();
      }
      if(decision[0] == "pass"){
        int playerNumber = myGame.trackTurn();
        myGame.pass(playerNumber);
        gui.updatePlayerStatus(1,myGame.getScore(1),myGame.getTracks(1),myGame.getGoals(1));
        gui.updatePlayerStatus(2,myGame.getScore(2),myGame.getTracks(2),myGame.getGoals(2));
        gui.updateMessage(myGame.getMessage());

        decision = gui.getNextMove();
      }
    }




    // Finally, clean up and exit successfully.
    delete graph;
    delete vertexPositions;
    return 0;
  }
