//=======================================================
// include guard
#ifndef __PARAMETERLIST_TD_H_INCLUDED__
#define __PARAMETERLIST_TD_H_INCLUDED__

// dependencies
#include <map>
#include <string>
#include <vector>
#include <sstream>


class ParameterList : public std::map<std::string, std::string>
{
public:
    ParameterList() { init(); };
    void init();
    template <class T> bool convert(const std::string s, T& result);
    void ParseVector(const std::string& s, std::vector<double>& v);
    void splitAndFill(const std::string& s);
};

template <class T> bool ParameterList::convert(const std::string s, T& result)
{
    std::string val = (*this)[s];
    std::stringstream ss(val);
    ss >> result;
    return val.empty();
}

#endif // __PARAMETERLIST_TD_H_INCLUDED__
