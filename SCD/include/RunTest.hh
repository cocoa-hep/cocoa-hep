enum test_t {
    kNONE = 0,
    kHADRON,
    kPI_ZERO,
    kJETS
};

/*
=============================
RunTest


Run tests to check code and physics.
Returns an exit code, so 0 is good, non-zero is bad.

=============================
*/

int RunTest( test_t type );
