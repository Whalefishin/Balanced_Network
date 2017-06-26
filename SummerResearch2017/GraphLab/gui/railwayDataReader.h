#pragma once

#include <istream>
#include <string>

#include "dataReader.h"

class RailwayDataReader : public DataReader {
public:
    RailwayDataReader(std::istream& input, bool debug=false);

    /**
     * Reads an edge from the stream.  Edges from the UI are stored only as a
     * pair between a source and a target.
     */
    std::pair<std::string,std::string> readEdge();

    /**
     * Reads a message from the stream.
     *
     * Messages are returned in the form of a list of strings.  This is possible
     * because all messages received from the UI contain only string data.  This
     * is not an ideal design, but it is simple enough to be used with ease by
     * CS35 students and does not involve the usual mess one would have by
     * modeling a closed variant in C++ using polymorphism and casting.
     *
     * Possible messages are:
     *   ["close"]: indicates that the user closed the game window
     *   ["error", msg]: indicates that a GUI error has occurred and no further
     *      messages will be received
     *   ["edge", source, target]: indicates that an edge was selected
     *   ["pass"]: indicates that the pass button was pressed
     */
    std::vector<std::string> readMessage();
};

