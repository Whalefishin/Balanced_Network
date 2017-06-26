#pragma once

#include <ostream>
#include <string>
#include <vector>
#include <utility>

#include "../goal.h"
#include "dataWriter.h"

class RailwayDataWriter : public DataWriter {
public:
    RailwayDataWriter(ostream& output, bool debug=false);

    /**
     * Writes the owner of a particular path on the board.  This is an integer
     * in the range of 0 through 2, where 0 indicates no owner.
     */
    void writeOwner(int owner);

    /**
     * Writes a point as an X-Y coordinate pair.
     */
    void writePoint(std::pair<int,int> point);

    /**
     * Writes a goal.
     */
    void writeGoal(Goal* goal);

    /**
     * Writes the status of a single player in the game.
     */
    void writePlayerState(int playerNumber, int score, int tracks, vector<Goal*> goals);

    /**
     * Writes a new game state message.
     */
    void writeGameMessage(std::string gameMessage);

    /**
     * Writes a new game graph.
     */
    void writeRailwayGraph(Graph<std::string,int,int>* graph);

    /**
     * Writes new vertex positions.
     */
    void writeVertexPositions(Dictionary<std::string, pair<int,int>>* positions);
};

