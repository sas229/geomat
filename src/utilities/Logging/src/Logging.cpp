#include "Logging.hpp"

void Logging::initialise_log(std::string severity) {
    // If no logger found, initialise the logger.
    auto log = plog::get();
    if (log == NULL) {
        static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
        if (severity == "verbose") {
            plog::init(plog::verbose, &consoleAppender);
        } else if (severity == "debug") {
            plog::init(plog::debug, &consoleAppender);
        } else if (severity == "info") {
            plog::init(plog::info, &consoleAppender);
        } else if (severity == "warning") {
            plog::init(plog::warning, &consoleAppender);
        } else if (severity == "error") {
            plog::init(plog::error, &consoleAppender);
        } else if (severity == "fatal") {
            plog::init(plog::fatal, &consoleAppender);
        } else if (severity == "none") {
            plog::init(plog::none, &consoleAppender);
        } else {
            plog::init(plog::warning, &consoleAppender);
            PLOG_WARNING << "Unrecognised log severity set. Possible options include \"verbose\", \"debug\", \"info\", \"warning\", \"error\", \"fatal\", or \"none\". Defaulting to \"error\".";
            plog::get()->setMaxSeverity(plog::error);
        }
        std::cout << "Log initialised..." << std::endl;
    } else {
        // Otherwise adjust the severity.
        log->setMaxSeverity(plog::severityFromString(severity.c_str()));
    }
}