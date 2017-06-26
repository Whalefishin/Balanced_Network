#include <istream>
#include <string>
#include <utility>

#include "formatError.h"
#include "railwayDataReader.h"

using namespace std;

RailwayDataReader::RailwayDataReader(std::istream& input, bool debug) : DataReader(input,debug) {
}

pair<string,string> RailwayDataReader::readEdge() {
    return readPair<string,string>(
        [&]() { return this->readString(); },
        [&]() { return this->readString(); }
    );
}

vector<string> RailwayDataReader::readMessage() {
    requireCharFromStream('M', "Expected message");
    string messageType = readString();
    vector<string> message = {messageType};
    if (messageType == "pass") {
        // Then there is nothing else to read.
    } else if (messageType == "edge") {
        pair<string,string> edge = readEdge();
        message.push_back(edge.first);
        message.push_back(edge.second);
    } else if (messageType == "error") {
        message.push_back(readString());
    } else if (messageType == "close") {
        // There's nothing else to read.
    }
    requireCharFromStream('$', "Expected message terminator");
    return message;
}
