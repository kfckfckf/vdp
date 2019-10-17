#include "VDPApp.h"
int main(int argc, char *argv[])
{
	auto && app = VDPApp::GetInstance();
	int code = 0;
	if (!app.ParseCommand(argc, argv))
	{
		code = 1;
	}
	else if (!app.Run())
	{
		code = 2;
	}
	return code;
}