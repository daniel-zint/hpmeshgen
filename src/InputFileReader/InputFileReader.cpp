#include "InputFileReader.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <functional>
#include <cctype>

// trim from start (in place)
static inline void ltrim( std::string &s ) {
	s.erase( s.begin(), std::find_if( s.begin(), s.end(), []( int ch ) {
		return !std::isspace( ch );
									  } ) );
}

// trim from end (in place)
static inline void rtrim( std::string &s ) {
	s.erase( std::find_if( s.rbegin(), s.rend(), []( int ch ) {
		return !std::isspace( ch );
						   } ).base(), s.end() );
}

void InputFileReader::deleteComments(std::string &str)
{
	if (str.find('!') == std::string::npos)
	{
		// no comment in this string
		return;
	}

	// delte comment
	str.erase(str.find('!'));

	// delete spaces at the beginning of string
	ltrim( str );
	// delete spaces at the end of string
	rtrim( str );
}
bool InputFileReader::getLine(std::ifstream& ifs, std::string& line)
{
	if (!std::getline(ifs, line))
		return false;

	deleteComments(line);
	while (line.size() == 0)
	{
		if (!std::getline(ifs, line))
		{
			return false;
		}

		deleteComments(line);
	}
	return true;
}

InputFileReader::InputFileReader(const std::experimental::filesystem::path& filename) : filename_(filename)
{
	std::ifstream ifs(filename_);
	
	LOG_ASSERT( ifs.is_open() ) << "Could not load file \"" << filename_ << "\"";

	std::string line;

	while (getLine(ifs, line)) {
		std::istringstream iss(line);

		std::string name;
		std::string value;
		std::string rest;
		iss >> name >> value >>rest;

		LOG_ASSERT( value.size() != 0 ) << "Invalid name/value-pair\nname = " << name << "\nvalue = " << value;
		LOG_ASSERT( rest.size() == 0 ) << "Invalid name/value-pair. There is more than one value in a line!\nname = " << name << "\nvalue = " << value;
		
		if( params_.find( name ) != params_.end() ) {
			LOG( WARNING ) << "Parameter '" << name << "' was set twice. Overwriting current value.";
		}
		params_[name] = value;
	}

	ifs.close();
}

InputFileReader::~InputFileReader()
{

}
