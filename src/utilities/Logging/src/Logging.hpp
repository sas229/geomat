#ifndef LOGGING_H
#define LOGGING_H

#include <plog/Log.h>
#include <plog/Init.h>
#include <plog/Formatters/TxtFormatter.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <string>

namespace Logging {
    /**
     * @brief Function to initialise a console log using plog. 
     * Possible options include \"verbose\", \"debug\", \"info\", \"warning\",
     * \"error\", \"fatal\", or \"none\". Defaults to \"error\"."
     * 
     * @param level Log level to use. 
     */
    void initialise_log(std::string level);
}

#endif