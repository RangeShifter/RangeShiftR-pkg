#include "Utils.h"

const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
const string Float2Str(const float x) {
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
const string Double2Str(const double x) {
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}

// Evaluate a lambda and assert we get the correct error 
void assert_error(const string& exptd_err_msg, void (*lambda)(void)) {
	string err_msg{ "No error.\n" };
	try { lambda(); } // evaluate
	catch (exception& e) { err_msg = e.what(); }
	assert(err_msg == exptd_err_msg);
}