#ifndef _INPUT_PARAMS_H
#define _INPUT_PARAMS_H

#include <map>
#include <string>
#include <stdexcept>
 
/// Exception of input parameters
class InputParamsError : public std::runtime_error {
public:
  InputParamsError(const std::string& msg) : std::runtime_error(msg) {};
};
 
/// Input parameters
class InputParams {
protected:
  std::map<std::string, std::string> values;

public:
  InputParams() {};
  
  void parse_command_line(int argc, char* argv[]); ///< Parses command line for parameters. Throws an exception if fails.
  
  const std::string& get_string(const char* key) const; ///< Return a parameters as a string. Throws an exception if fails.
  const int get_int32(const char* key) const; ///< Return a parameters as an integer. Throws an exception if fails.
  const double get_double(const char* key) const; ///< Return a parameters as a double. Throws an exception if fails.
  const bool get_bool(const char* key, const bool needs_definition = true) const; ///< Return a parameters as a bool. In a case needs_definition is false, it returns true if parameter is defined and false if parameter is not defined. Otherwise it requires explicit value 'true' or 'false'. In this case, it is not case-sensitive.
  
  void set_default(const char* key, const char* value); ///< Set a value of a parameter.
  void set_default(const char* key, const int value); ///< Set a value of a parameter.
  void set_default(const char* key, const double value); ///< Set a value of a parameter.
  void set_default(const char* key, const bool value); ///< Set a value of a parameter.
};

#endif
