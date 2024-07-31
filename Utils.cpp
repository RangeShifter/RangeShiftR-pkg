#include "Utils.h"

// Evaluate a lambda and assert we get the correct error 
void assert_error(const string& exptd_err_msg, void (*lambda)(void)) {
	string err_msg{ "No error.\n" };
	try { lambda(); } // evaluate
	catch (exception& e) { err_msg = e.what(); }
	assert(err_msg == exptd_err_msg);
}