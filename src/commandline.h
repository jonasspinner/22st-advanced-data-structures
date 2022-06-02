#ifndef ADVANCED_DATA_STRUCTURES_COMMANDLINE_H
#define ADVANCED_DATA_STRUCTURES_COMMANDLINE_H

#include <string>
#include <stdexcept>

namespace ads {
    struct ParseResult {
        enum class AlgoKind {
            bv, bp
        } algo{};
        std::string input;
        std::string output;
    };

    ParseResult parse(int argc, char *argv[]) {
        if (argc != 4) {
            throw std::invalid_argument("must contain 4 arguments");
        }
        ParseResult result;
        if (std::string_view(argv[1]) == "bv") {
            result.algo = ParseResult::AlgoKind::bv;
        } else if (std::string_view(argv[1]) == "bp") {
            result.algo = ParseResult::AlgoKind::bp;
        } else {
            throw std::invalid_argument("wrong algorithm name");
        }
        result.input = argv[2];
        result.output = argv[3];
        return result;
    }
}

#endif //ADVANCED_DATA_STRUCTURES_COMMANDLINE_H
