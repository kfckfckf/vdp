#pragma once
#include "ParameterBase.h"
#include <vector>
#include <sstream>
#include <iostream>

template<typename T>
struct ParameterTag {};

template<typename T>
class ParameterT :public ParameterBase
{
public:
	typedef T DataType;

	ParameterT(const std::string& name, const std::string& keyword, const std::string& shortcut, const DataType& defaultvalue)
		:ParameterBase(name, keyword, shortcut), data(defaultvalue) {}

	void SetValue(const DataType& value) { data = value; }

	const DataType& Value(void) const { return data; }

	bool ReadValue(const std::string& value) override
	{
		return ReadValue(value, ParameterTag<DataType>());
	}

	int Size(void) const override { return Size(ParameterTag<T>()); }

	void PrintUsage(void) const override
	{
		PrintUsage(ParameterTag<DataType>());
	}

private:
	template<typename Type>
	bool ReadValue(const std::string& value, ParameterTag<Type>)
	{
		std::istringstream iss(value);
		iss >> data;
		return !iss.fail();
	}

	template<>
	bool ReadValue(const std::string&, ParameterTag<bool>)
	{
		data = true;
		return true;
	}

	template<>
	bool ReadValue(const std::string& value, ParameterTag<std::string>)
	{
		std::istringstream iss(value);
		std::getline(iss, data);
		return !iss.fail();
	}

	template<typename ValueType>
	bool ReadValue(const std::string& value, ParameterTag<std::vector<ValueType>>)
	{
		std::istringstream iss(value);
		for (size_t i = 0; i < data.size(); i++)
		{
			iss >> data[i];
			if (iss.fail())
			{
				return false;
			}
		}
		return true;
	}

	template<>
	bool ReadValue(const std::string& value, ParameterTag<std::vector<std::string>>)
	{
		std::istringstream iss(value);
		for (size_t i = 0; i < data.size(); i++)
		{
			std::getline(iss, data[i]);
			if (iss.fail())
			{
				return false;
			}
		}
		return true;
	}

	template<typename Type>
	int Size(ParameterTag<Type>) const
	{
		return 1;
	}

	template<>
	int Size(ParameterTag<bool>) const
	{
		return 0;
	}

	template<typename Type>
	int Size(ParameterTag<std::vector<Type>>) const
	{
		return data.size();
	}

	template<typename Type>
	void PrintUsage(ParameterTag<Type>) const
	{
		std::cout << " <arg>";
	}

	template<>
	void PrintUsage(ParameterTag<bool>) const
	{
	}

	template<typename Type>
	void PrintUsage(ParameterTag<std::vector<Type>>) const
	{
		for (size_t i = 0; i < data.size(); ++i)
		{
			std::cout << " <arg" << i + 1 << ">";
		}
	}

	DataType data;
};
