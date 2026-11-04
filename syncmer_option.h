#pragma once

#ifndef _SYNCMER_OPTION_H_
#define _SYNCMER_OPTION_H_

#include <CLI/CLI.hpp>
#include <string>

#include "asm.h"

// Registers --syncmer k/s[/t] on a subcommand, parsed into `out` with
// syncmer_params::parse. `out` starts at the default open syncmer (11/4/6),
// which the help shows.
inline CLI::Option *add_syncmer_option(CLI::App &app, syncmer_params &out) {
    out.apply_default_rule();
    return app
        .add_option_function<std::string>(
            "--syncmer",
            [&out](const std::string &text) {
                auto parsed = syncmer_params::parse(text);
                if (!parsed) throw CLI::ValidationError(parsed.error());
                out = *parsed;
            },
            "Syncmer selection: k-mer length k, s-mer length s and offset t "
            "of the smallest s-mer (open syncmer); without /t a closed "
            "syncmer (smallest s-mer at the first or last position). "
            "Lengths count homopolymer-compressed bases.")
        ->type_name("K/S[/T]")
        ->default_str(out.text());
}

#endif
