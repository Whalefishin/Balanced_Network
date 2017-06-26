/*
 * Declares RailwayGame class
 *
 * Bill Huang & Genji Kawakita
 *
 */


#pragma once

#include <vector>
#include "adjacencyListGraph.h"
#include "goal.h"
#include "stl/stlHashTable.h"
#include "graphAlgorithms-inl.h"
#include "stl/stlBST.h"

using namespace std;

/**
 * An instance of this class represents a single game of Railway.  An object of
 * type RailwayGame does not represent the user interface, but does represent
 * the idea of a game of Railway, including the state of the board and all of
 * the information pertaining to both players.  It also includes the logic for
 * making moves in the game and enforces the game's rules.
 */
class RailwayGame {
public:


    RailwayGame();

    ~RailwayGame();

    /**
     * Randomly choose cities to be part of P1/P2's goals. Attach to each goal
     * an appropriate amount of score for its completion.
     * The randomly chosen cities must be at least three routes apart.
     */
    void buildGoals();


    /**
     * Claims the edge given by src and dest for either Player1 or Player2,
     * depending on whose turn it is.
     * After claiming the edge, update the list of edges P1/P2 can claim.
     * Also update the railwayMap to be reflected on the GUI.
     * Also update the scores and the message that will be reflected in the GUI.
     * This is a really long method, in retrospect we should have written some
     * helper functions..
     */
    void claimEdge(string src, string dest);

    /**
     * Update the data members of Goal objects and check reachability of
     * each Goal.
     * Update the message to reflect if a player has completed a goal, or if
     * that goal is no longer possible to achieve.
     */
    void updateGoals();

    /**
    * Returns 1 if it is Player1's turn, 2 if it is Player2's turn.
    * @return a int that indicates whose turn it is
    */
    int trackTurn();

    /**
     * Gets the score of a given player
     * @param playerNumber(1 or 2)
     * @return the value of given player's score
     */
    int getScore(int playerNumber);

    /**
     * Gets the tracks of a given player
     * @param playerNumber(1 or 2)
     * @return the number of given player's tracks
     */
    int getTracks(int playerNumber);

    /**
     * Gets the goals of a given player
     * @param playerNumber(1 or 2)
     * @return vector of pointers to a given player's goals
     */
    vector<Goal*> getGoals(int playerNumber);

    /**
     * Gets message shown at the bottom of the window
     * @return message to be shown on the window
     */
    string getMessage();

    /**
     * Updates player's turn and tracks when that player passes his turn.
     * @param playerNumber (1 or 2)
     */
    void pass(int playerNumber);

    /**
     * Initializes all the five graphs we need, including the railwayMap in the GUI
     * the map in P1/P2's perspective, and the map that contains all the claimable
     * edges for P1/P2
     * @param the map to load
     */
    void loadMap(Graph<string,int,int>* g);

    /**
     * Returns the railwayMap to be reflected in the GUI
     * @return message to be shown on the window
     */
    Graph<string,int,int>* getMap();

    /**
     * Returns the winner's player number.
     * @return 1 if P1 has won, 2 if P2 has won.
     */
    int getWinner();

    /**
     * Checks if the game is over.
     * @return true if the game is over. False otherwise.
     */
    bool gameOver();



private:
    int p1Score;
    int p2Score;
    int p1Tracks;
    int p2Tracks;
    vector<Goal*> p1Goals;
    vector<Goal*> p2Goals;
    Graph<string,int,int>* railwayMap;
    Graph<string,int,int>* p1Map;
    Graph<string,int,int>* p2Map;
    int turn_parity;

    Graph<string,int,int>* p1ClaimableEdges;
    Graph<string,int,int>* p2ClaimableEdges;

    STLHashTable<string,bool> p1Cities;
    STLHashTable<string,bool> p2Cities;

    string bottomMessage;
    int winner_number;
};
