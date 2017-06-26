#pragma once

// Example of communication with a subprocess via stdin/stdout
// Original Author: Konstantin Tretyakov
// License: MIT
// https://gist.github.com/konstantint/d49ab683b978b3d74172

#include <ext/stdio_filebuf.h> // NB: Specific to libstdc++
#include <sys/wait.h>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <exception>

// Wrapping pipe in a class makes sure they are closed when we leave scope
class cpipe {
private:
    int fd[2];
public:
    cpipe();
    ~cpipe();
    const int read_fd() const;
    const int write_fd() const;
    void close();
};

