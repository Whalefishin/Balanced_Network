#pragma once

#include "guiFailure.h"

class FormatError : public GuiFailure {
public:
    explicit FormatError(const char* message) : GuiFailure(message) { }
    explicit FormatError(std::string& message) : GuiFailure(message) { }
    virtual ~FormatError() throw (){}
};

