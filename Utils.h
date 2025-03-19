#ifndef UtilsH
#define UtilsH

#include <stdlib.h>
#include <iostream>
#include <sstream>

#include <cassert>
using namespace std;

// Evaluate a lambda and assert we get the correct error 
void assert_error(const string& exptd_err_msg, void (*x)(void));

#endif // UtilsH
