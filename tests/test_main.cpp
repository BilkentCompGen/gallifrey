// Entry point for the lasm_tests executable.
//
// The TEST_CASE bodies themselves live inline in the production translation
// units (asm2.cpp, lcmer.cpp) because they exercise TU-internal functions.
// The production `lasm` binary is compiled with DOCTEST_CONFIG_DISABLE so those
// bodies are stripped from it; this target compiles the same TUs without that
// define and provides doctest's main here.
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "ext/doctest.h"

// Provide the header-only lcptools implementation for the test binary, the same
// way lasm.cpp does for the production binary (the test target doesn't compile
// lasm.cpp).
#define LCPTOOLS_IMPL
#include "ext/lcptools/lcptools_ho.h"
