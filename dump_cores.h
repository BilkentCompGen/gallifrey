#pragma once

#ifndef _DUMP_CORES_H_
#define _DUMP_CORES_H_

#include <cstdint>
#include <string>

#include "asm.h"

namespace CLI {
class App;
}

// "dump-cores" subcommand: owns its parameters and registers itself with the
// parent app on construction (see gallifrey.cpp main). Copy/move deleted so the
// callback's captured `this` cannot dangle.
class DumpCoresProgram {
    std::string input;
    std::string output;
    syncmer_params syncmer;
    int lcp_level = 4;
    int lcmer_l = 0;
    uint32_t seed = kDefaultRandomSeed;

public:
    explicit DumpCoresProgram(CLI::App &parent);
    DumpCoresProgram(const DumpCoresProgram &) = delete;
    DumpCoresProgram &operator=(const DumpCoresProgram &) = delete;
};

#endif
