#include "VDPApp.h"
int main(int argc, char* argv[])
{
	auto& app = VDPApp::GetInstance();
	if (!app.ParseCommand(argc, argv))
	{
		return 1;
	}
	if (!app.Run())
	{
		return 2;
	}
	return 0;
}
