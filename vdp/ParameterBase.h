#pragma once
#include <string>

class ParameterBase
{
public:
	ParameterBase(const std::string& name, const std::string& keyword, const std::string& shortcut);
	virtual ~ParameterBase(void);
	const std::string& Name(void) const;
	const std::string& Keyword(void) const;
	const std::string& Shortcut(void) const;
	bool IsProvided(void) const;
	void SetProvided(bool b = true);
	virtual bool ReadValue(const std::string& value) = 0;
	virtual int Size(void) const = 0;
	virtual void PrintUsage(void) const = 0;
protected:
	std::string name;
	std::string keyword;
	std::string shortcut;
	bool isprovided = false;
};
