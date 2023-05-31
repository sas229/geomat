#include "Logging.hpp"

void Logging::initialise_log(std::string level) {
    static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
    if (level == "verbose") {
        plog::init(plog::verbose, &consoleAppender);
    } else if (level == "debug") {
        plog::init(plog::debug, &consoleAppender);
    } else if (level == "info") {
        plog::init(plog::info, &consoleAppender);
    } else if (level == "warning") {
        plog::init(plog::warning, &consoleAppender);
    } else if (level == "error") {
        plog::init(plog::error, &consoleAppender);
    } else if (level == "fatal") {
        plog::init(plog::fatal, &consoleAppender);
    } else if (level == "none") {
        plog::init(plog::none, &consoleAppender);
    } else {
        plog::init(plog::warning, &consoleAppender);
        PLOG_WARNING << "Unrecognised log level set. Possible options include \"verbose\", \"debug\", \"info\", \"warning\", \"error\", \"fatal\", or \"none\". Defaulting to \"error\".";
        plog::get()->setMaxSeverity(plog::error);
    }
}