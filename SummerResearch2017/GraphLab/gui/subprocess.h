#pragma once

#include <ext/stdio_filebuf.h> // NB: Specific to libstdc++
#include <sys/wait.h>
#include <unistd.h>
#include <iostream>
#include <memory>
#include <exception>

#include "cpipe.h"

// Example of communication with a subprocess via stdin/stdout
// Original Author: Konstantin Tretyakov
// License: MIT
// https://gist.github.com/konstantint/d49ab683b978b3d74172

class Subprocess {
public:
    Subprocess(const char* const argv[], bool with_path = false, const char* const envp[] = 0);

    void sendEOF();
    int wait();
    void kill();

    int child_pid = -1;
    std::unique_ptr<__gnu_cxx::stdio_filebuf<char> > write_buf = NULL;
    std::unique_ptr<__gnu_cxx::stdio_filebuf<char> > read_buf = NULL;
    std::ostream stdin;
    std::istream stdout;

private:
    cpipe write_pipe;
    cpipe read_pipe;

    Subprocess(const Subprocess& other) = delete;
    Subprocess& operator=(const Subprocess& other) = delete;
};

