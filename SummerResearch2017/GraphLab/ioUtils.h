#pragma once

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "adts/dictionary.h"
#include "adts/edge.h"
#include "adts/graph.h"

/**
 * Reads a vertex position map from the specified file.
 * @param filename The file to read.
 * @return A pointer to a dictionary containing the vertex positions.  The
 *         caller becomes the owner of this memory.
 * @throws std::exception If an I/O error occurs.
 */
Dictionary<std::string, std::pair<int,int>>*
    readVertexPositions(std::string filename);

/**
 * Reads a Railway graph from the specified file.  The owner for all edges is
 * set to 0 (no owner).
 * @param filename The file to read.
 * @return A pointer to the graph which was loaded.  The caller becomes the
 *         owner of this memory.
 * @throws std::exception If an I/O error occurs.
 */
Graph<std::string, int, int>* readRailwayGraph(std::string filename);

