#include "goal.h"

Goal::Goal(std::string location1, std::string location2, std::string message, int points) {
    this->location1 = location1;
    this->location2 = location2;
    this->message = message;
    this->points = points;
}
