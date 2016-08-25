/*
 * =====================================================================================
 *
 *       Filename:  parameterlist.h
 *
 *    Description: Header for class that imports the input file, parses
 *                 and converts to parameters for the model.
 *
 *        Version:  1.0
 *        Created:  08/24/2016 16:50:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Thomas McDonald (), mcdonald@jimmy.harvard.edu
 *   Organization:  DFCI
 *
 * =====================================================================================
 */

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


#endif // __PARAMETERLIST_TD_H_INCLUDED__
