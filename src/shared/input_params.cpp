#include <sstream>
#include <stdio.h>
#include "input_params.h"

using namespace std;

void InputParams::parse_command_line(int argc, char *argv[])
{
  int par_inx = 1;
  while (par_inx < argc) {
    //get parameter name
    const char* key = argv[par_inx];
    if (key[0] != '-') {
      stringstream s;
      s << "excpected a parameter but value \"" << key << "\" found";
      throw InputParamsError(s.str());
    }
    if (key[1] == '\0') {
      stringstream s;
      s << "excpected a parameter name but a noname parameter at index " << par_inx << " found";
      throw InputParamsError(s.str());
    }
    key++;

    //get value
    const char* val = NULL;
    if ((par_inx + 1) < argc)
      val = argv[par_inx + 1];

    //check value
    if (val != NULL && val[0] != '-') { //there is some value
      values[key] = val;
      par_inx++;
    } else {
      values[key] = ""; //no value, just boolean
    }
    par_inx++;
  }
}

const string& InputParams::get_string(const char* key) const {
  map<string, string>::const_iterator iter_value = values.find(key);
  if (iter_value != values.end())
    return iter_value->second;
  else {
    stringstream s;
    s << "parameter \"" << key << "\" not found";
    throw InputParamsError(s.str());
  }
}

const int InputParams::get_int32(const char* key) const {
  const string& str_value = get_string(key);
  int value = 0;
  if (sscanf(str_value.c_str(), "%d", &value) == 1 || sscanf(str_value.c_str(), "0x%x", &value) == 1 || sscanf(str_value.c_str(), "0x%X", &value) == 1)
    return value;
  else {
    stringstream s;
    s << "parameter \"" << key << "\", value \"" << str_value << "\" is not integer";
    throw InputParamsError(s.str());
  }
}

const double InputParams::get_double(const char* key) const {
  const string& str_value = get_string(key);
  double value = 0;
  if (sscanf(str_value.c_str(), "%lf", &value) == 1)
    return value;
  else {
    stringstream s;
    s << "parameter \"" << key << "\", value \"" << str_value << "\" is not a floating point";
    throw InputParamsError(s.str());
  }
}

const bool InputParams::get_bool(const char* key, const bool needs_definition) const {
  if (needs_definition) {
    const string& str_value = get_string(key);

    //convert to lower case
    stringstream str_low_value;
    string::const_iterator iter = str_value.begin();
    while (iter != str_value.end()) {
      str_low_value << (char)tolower(*iter);
      iter++;
    }

    //compare
    if (str_low_value.str().compare("true") == 0)
      return true;
    else if (str_low_value.str().compare("false") == 0)
      return false;
    else {
      stringstream s;
      s << "parameter \"" << key << "\", value \"" << str_low_value << "\" is not boolean";
      throw InputParamsError(s.str());
    }
  }
  else {
    try { 
      get_string(key);
      return true;
    }
    catch(InputParamsError&) { return false; }
  }
}

void InputParams::set_default(const char* key, const char* value) {
  values[key] = value;
}

void InputParams::set_default(const char* key, const int value) {
  stringstream s;
  s << value;
  values[key] = s.str();  
}

void InputParams::set_default(const char* key, const double value) {
  stringstream s;
  s << scientific << value;
  values[key] = s.str();  
}

void InputParams::set_default(const char* key, const bool value) {
  if (value)
    values[key] = "true";
  else
    values[key] = "false";
}
