#pragma once

/* InputFileReader

	Store your parameters in an input file:
	- comments start with a !
	- empty lines are allowed
	- always write: parameter value
	- per parameter only one value

	Get parameters with .getValue<type>(parameter);
*/

#include <string>
#include <map>
#include <experimental/filesystem>

#include <glog/logging.h>


class InputFileReader
{
	std::experimental::filesystem::path filename_;
	std::map<std::string, std::string> params_;

public:
	// Constructors
	InputFileReader() = delete;
	InputFileReader(const std::experimental::filesystem::path &filename);
	InputFileReader(const InputFileReader &) = delete;

	// Destructor
	~InputFileReader();

	// Operators
	InputFileReader& operator=(const InputFileReader &) = delete;

	// Data-access
	template<typename T>
	T getValue(const std::string &param);

private:
	void deleteComments(std::string &str);
	bool getLine(std::ifstream& ifs, std::string& line);
	void errorMsg(const std::string & msg, const std::string& file, const std::string& function, const int& line) const;
};

// Data-access
template<typename T>
T InputFileReader::getValue(const std::string &param)
{
	LOG_ASSERT( param.length() != 0 ) << "No parameter name given";
	LOG_ASSERT( params_.find( param ) != params_.end() ) << "No parameter named \"" << param << "\"";

	std::stringstream pstr;
	pstr << params_.at(param);
	T retVal;
	pstr >> retVal;
	return retVal;
}