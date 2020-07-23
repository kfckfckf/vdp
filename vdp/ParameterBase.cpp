#include "ParameterBase.h"

ParameterBase::ParameterBase(const std::string& name, const std::string& keyword, const std::string& shortcut)
	:name(name), keyword(keyword), shortcut(shortcut)
{
}

ParameterBase::~ParameterBase(void)
{
}

const std::string& ParameterBase::Name(void) const
{
	return name;
}

const std::string& ParameterBase::Keyword(void) const
{
	return keyword;
}

const std::string& ParameterBase::Shortcut(void) const
{
	return shortcut;
}

bool ParameterBase::IsProvided(void) const
{
	return isprovided;
}

void ParameterBase::SetProvided(bool b)
{
	isprovided = b;
}
