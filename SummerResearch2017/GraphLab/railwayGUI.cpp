#include <stdexcept>

#include "railwayGUI.h"

using namespace std;

//const char* const GUI_ARGS[] = {"java","-jar","RailwaysGUI.jar",NULL};

RailwayGUI::RailwayGUI(
        std::string backgroundImageFilename,
        Dictionary<string,pair<int,int>>* vertexPositions,
        bool debug) {
    this->guiArgs = new const char*[5 + (debug?1:0)];
    int argIdx = 0;
    this->guiArgs[argIdx++] = "java";
    this->guiArgs[argIdx++] = "-jar";
    this->guiArgs[argIdx++] = "RailwayGUI.jar";
    if (debug) {
        this->guiArgs[argIdx++] = "-D";
    }
    this->guiArgs[argIdx++] = backgroundImageFilename.c_str();
    this->guiArgs[argIdx++] = NULL;

    this->subprocess = new Subprocess(this->guiArgs, true);
    this->reader = new RailwayDataReader(this->subprocess->stdout, debug);
    this->writer = new RailwayDataWriter(this->subprocess->stdin, debug);

    this->writer->writeVertexPositions(vertexPositions);
}

RailwayGUI::~RailwayGUI() {
    stopGUI();
    delete[] this->guiArgs;
}

void RailwayGUI::updateRailwayMap(Graph<string,int,int>* graph) {
    if (this->subprocess != NULL) {
        this->writer->writeRailwayGraph(graph);
    } else {
        throw runtime_error("The GUI is no longer functioning; cannot update map.");
    }
}

void RailwayGUI::updateMessage(string message) {
    if (this->subprocess != NULL) {
        this->writer->writeGameMessage(message);
    } else {
        throw runtime_error("The GUI is no longer functioning; cannot update message.");
    }
}

void RailwayGUI::updatePlayerStatus(
        int playerNumber, int score, int tracks, vector<Goal*> goals) {
    if (this->subprocess != NULL) {
        this->writer->writePlayerState(playerNumber, score, tracks, goals);
    } else {
        throw runtime_error("The GUI is no longer functioning; cannot update player.");
    }
}

vector<string> RailwayGUI::getNextMove() {
    if (this->subprocess != NULL) {
        vector<string> message = this->reader->readMessage();
        if (message.size() > 0 && message[0] == "close") {
            stopGUI();
        } else if (message.size() > 0 && message[0] == "error") {
            stopGUI();
            throw runtime_error("GUI error: " + message.at(1));
        }
        return message;
    } else {
        throw runtime_error("Attempted to get next move from GUI which is " \
                "no longer functioning.");
    }
}

void RailwayGUI::stopGUI() {
    if (this->subprocess != NULL) {
        this->subprocess->kill();
        delete this->subprocess;
        delete this->reader;
        delete this->writer;
        this->subprocess = NULL;
    }
}
