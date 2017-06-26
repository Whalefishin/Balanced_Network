#pragma once

#include <string>

#include "adts/graph.h"

#include "goal.h"
#include "gui/railwayDataReader.h"
#include "gui/railwayDataWriter.h"
#include "gui/subprocess.h"

class RailwayGUI {
public:
    /**
     * Creates a new Railways GUI.  This produces a GUI window for the user to
     * interact with.
     * @param backgroundImageFilename The name of the file to use as the
     *                                background image.
     * @param vertexPositions The positions of the vertices in the upcoming
     *                        graph.
     * @param debug true to put the GUI into debug mode.
     */
    RailwayGUI(
        std::string backgroundImageFilename,
        Dictionary<std::string,pair<int,int>>* vertexPositions,
        bool debug = false);

    /**
     * Destroys the Railways GUI.  This action closes the window if it is still
     * open.
     */
    ~RailwayGUI();

    /**
     * Updates the graph being displayed by the GUI.  This method copies the
     * contents of the provided graph into the display.  It must be called when
     * the map should be updated even if the same graph pointer is used.
     * @param graph The graph that should be used to update the GUI.  The caller
     *              retains ownership of the graph.
     */
    void updateRailwayMap(Graph<std::string,int,int>* graph);

    /**
     * Updates the message displayed in the GUI window.
     * @param message The new message to display.
     */
    void updateMessage(std::string message);

    /**
     * Updates the status of a player in the GUI window.
     * @param playerNumber The player to update: 1 or 2.
     * @param score The player's score.
     * @param tracks The number of tracks the player has.
     * @param goals The player's goals.
     */
    void updatePlayerStatus(
        int playerNumber, int score, int tracks, std::vector<Goal*> goals);

    /**
     * Retrieves the next move from the GUI.  This method will wait until the
     * user provides some input and then produce a return value that describes
     * that input; the result is always in the form of a list of strings for the
     * sake of simplicity.  The possible outcomes are:
     *   1. A vector containing just the string "pass", if the player has
     *      passed.
     *   2. A vector containing just the string "close", if the player has
     *      closed the window (meaning that the program should now stop and no
     *      further moves will be made).
     *   3. A vector containing the string "edge" followed by the names of the
     *      two vertices between which that edge has been placed, if the player
     *      has selected an edge.
     * If anything goes wrong during this process (e.g. getting a move from a
     * GUI that has been closed), an exception is thrown.
     * @return The next move in the GUI as described above.
     */
    std::vector<std::string> getNextMove();
private:
    const char** guiArgs;

    Subprocess* subprocess;
    RailwayDataReader* reader;
    RailwayDataWriter* writer;

    void stopGUI();
};

