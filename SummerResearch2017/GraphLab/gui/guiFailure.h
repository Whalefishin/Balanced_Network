#pragma once

#include <stdexcept>
#include <string>

class GuiFailure : public std::runtime_error {
public:
    explicit GuiFailure(const char* message) : std::runtime_error(message) { }
    explicit GuiFailure(std::string& message) : std::runtime_error(message) { }
    virtual ~GuiFailure() throw (){}
};

