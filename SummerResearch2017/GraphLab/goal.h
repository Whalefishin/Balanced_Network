#pragma once

#include <string>

class Goal {
public:
    Goal(std::string location1, std::string location2, std::string message, int points);

    /** One of the locations in the objective. */
    std::string location1;
    /** The other location in the objective. */
    std::string location2;
    /** A message describing the state of the objective. */
    std::string message;
    /** the points that user can get from completeing the path */
    int points;

    // You can safely ignore the following code.  This eliminates some default
    // class routines, preventing you from using them accidentally.
    // Specifically, we are disabling the use of the copy constructor and the
    // copy assignment operator.  You can read more here:
    //   http://www.cplusplus.com/articles/y8hv0pDG/
    Goal(const Goal& goal) = delete;
    Goal& operator=(const Goal& goal) = delete;
};
