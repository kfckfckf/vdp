#pragma once
#include "ParameterT.h"
#include <map>
#include <memory>

class ParameterManager
{
public:
	ParameterManager(void);
	~ParameterManager(void);
	void Clear(void);
	template<typename T>
	void AddParameter(const std::string& name, const std::string& keyword, const std::string& shortcut, const T& defaultvalue, bool isoptional = true)
	{
		if (isoptional)
		{
			optionalParameters.insert({ name, std::make_shared<ParameterT<T>>(name, keyword, shortcut, defaultvalue) });
		}
		else
		{
			parameters.push_back(std::make_shared<ParameterT<T>>(name, keyword, shortcut, defaultvalue));
		}
	}
	template<typename T>
	auto GetParameter(const std::string& name) const
	{
		if (auto p = GetParameter(name))
		{
			return std::dynamic_pointer_cast<ParameterT<T>>(p)->Value();
		}
		else
		{
			std::cout << "Warning: no parameter " << name << ".\n";
			return T();
		}
	}

	template<typename T>
	bool SetParameter(const std::string& name, const T& value)
	{
		if (auto p = GetParameter(name))
		{
			std::dynamic_pointer_cast<ParameterT<T>>(p)->SetValue(value);
			return true;
		}
		else
		{
			return false;
		}
	}

	bool IsProvided(const std::string& name) const;
	bool ParseCommand(int argc, char* argv[]);
	void PrintUsage(void) const;
private:
	bool ParseParameter(std::shared_ptr<ParameterBase> p, int i, int argc, char* argv[]);
	std::shared_ptr<ParameterBase> GetParameter(const std::string& name) const;
	std::map<std::string, std::shared_ptr<ParameterBase>> optionalParameters;
	std::vector<std::shared_ptr<ParameterBase>> parameters;
	std::string processname;
};

