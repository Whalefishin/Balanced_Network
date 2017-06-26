#include <functional>
#include <ostream>
#include <string>
#include <vector>
#include <utility>

#include "../adts/dictionary.h"
#include "../adts/edge.h"
#include "../adts/graph.h"
#include "../adts/list.h"
#include "dataWriter.h"

using namespace std;

DataWriter::DataWriter(ostream& output, bool debug) : output(output) {
    this->debug = debug;
    output.exceptions(std::ostream::failbit | std::ostream::badbit);
}

void DataWriter::writeCharToStream(char c) {
    output << c;
    if (debug) {
        cerr << "sent: " << c << endl;
    }
}

void DataWriter::writeStringToStream(string s) {
    output << s;
    if (debug) {
        cerr << "sent: " << s << endl;
    }
}

void DataWriter::flushStream() {
    output.flush();
}

void DataWriter::writeInt(int i) {
    writeCharToStream('i');
    writeStringToStream(to_string(i));
    writeCharToStream('$');
}

void DataWriter::writeString(string s) {
    writeCharToStream('s');
    writeInt(s.length());
    writeStringToStream(s);
    writeCharToStream('$');
}
