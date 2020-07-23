#include "ParameterManager.h"

ParameterManager::ParameterManager(void)
{
}

ParameterManager::~ParameterManager(void)
{
}

void ParameterManager::Clear(void)
{
	parameters.clear();
	optionalParameters.clear();
}

bool ParameterManager::IsProvided(const std::string& name) const
{
	if (auto p = GetParameter(name))
	{
		return p->IsProvided();
	}
	else
	{
		return false;
	}
}

bool ParameterManager::ParseCommand(int argc, char* argv[])
{
	processname = argv[0];
	auto slash = processname.find_last_of("\\/");
	if (slash != std::string::npos)
	{
		processname = processname.substr(slash + 1);
	}
	processname = processname.substr(0, processname.find_last_of("."));
	
	int pid = 0;
	for (int i = 1; i < argc; ++i)
	{
		if (argv[i][0] == '-')
		{
			if (argv[i][1] == '-') // keyword
			{
				std::string s = &argv[i][2];
				bool bfound = false;
				for (const auto& p : optionalParameters)
				{
					if (s == p.second->Keyword())
					{
						if (!ParseParameter(p.second, i + 1, argc, argv))
						{
							return false;
						}
						bfound = true;
						i += p.second->Size();
						break;
					}
				}
				if (!bfound)
				{
					std::cout << "Warning: unknown arguments.\n";
				}
			}
			else // shortcut
			{
				std::string s = &argv[i][1];
				bool bfound = false;
				for (const auto& p : optionalParameters)
				{
					if (s == p.second->Shortcut())
					{
						if (!ParseParameter(p.second, i + 1, argc, argv))
						{
							return false;
						}
						bfound = true;
						i += p.second->Size();
						break;
					}
				}
				if (!bfound)
				{
					std::cout << "Warning: unknown arguments.\n";
				}
			}
		}
		else
		{
			if (pid >= parameters.size())
			{
				std::cout << "Error: unknown parameter.\n";
				return false;
			}
			if (!ParseParameter(parameters[pid], i, argc, argv))
			{
				return false;
			}
			i += parameters[pid]->Size() ? (parameters[pid]->Size() - 1) : 0;
			++pid;
		}
	}
	for (const auto& p : parameters)
	{
		if (!p->IsProvided())
		{
			std::cout << "Error: not enough arguments.\n";
			return false;
		}
	}
	return true;
}

void ParameterManager::PrintUsage(void) const
{
	std::cout << "Usage: " << processname;
	for (const auto& p : parameters)
	{
		std::cout << " <" << p->Keyword() << ">";
	}
	for (const auto& p : optionalParameters)
	{
		std::cout << " [--" << p.second->Keyword();
		p.second->PrintUsage();
		std::cout << "]";
	}
	std::cout << std::endl;
}

bool ParameterManager::ParseParameter(std::shared_ptr<ParameterBase> p, int i, int argc, char* argv[])
{
	int size = p->Size();
	std::string strparameters;
	for (int j = 0; j < size; ++j)
	{
		if (i + j < argc)
		{
			strparameters += argv[i + j];
			strparameters += '\n';
		}
		else
		{
			std::cout << "Error: not enough parameter.\n";
			return false;
		}
	}
	if (!p->ReadValue(strparameters))
	{
		std::cout << "Error: wrong parameter.\n";
		return false;
	}
	p->SetProvided();
	return true;
}

std::shared_ptr<ParameterBase> ParameterManager::GetParameter(const std::string& name) const
{
	std::map<std::string, std::shared_ptr<ParameterBase>>::const_iterator iter;
	if ((iter = optionalParameters.find(name)) != optionalParameters.end())
	{
		return iter->second;
	}
	else
	{
		for (const auto& p : parameters)
		{
			if (p->Name() == name)
			{
				return p;
			}
		}
		return nullptr;
	}
}
