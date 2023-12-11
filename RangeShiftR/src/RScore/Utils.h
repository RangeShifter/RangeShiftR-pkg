#ifndef UtilsH
#define UtilsH

#include <stdlib.h>
#include <iostream>
#include <sstream>

#include <cassert>
using namespace std;

const string Int2Str(const int x);
const string Float2Str(const float x);
const string Double2Str(const double x);

// Evaluate a lambda and assert we get the correct error 
void assert_error(const string& exptd_err_msg, void (*x)(void));

#endif // UtilsH
