/*!
 * class Configurator
 *
 * modified from old class Configurator
 *
 *
 * Author: Jian Zhou
 * Date: 11-21-2013
 *
 */
#ifndef PETSYS_CFG_H
#define PETSYS_CFG_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <map>

#include <petsys_log.h>
#include <petsys_sysparms.h>

namespace UCDPETSYS
{

class Configurator
{
public:
    Configurator();
    ~Configurator();

public:
    void processOptionFile(const char* filename, SYSTEMPARAMETERS& sysparms);
    template<typename T>
    static T string_to_T(const std::string& s);

private:
    enum OP_state {
        OP_OK = 0,
        OP_EMPTY,
        OP_COMMENT_ONLY,
        OP_BAD_FORMAT,
        OP_UNKNOWN
    };

    static const char va_delim = ',';
    static const char eq_delim = '=';
    static const char cm_delim = '#';
    static const char char_space = ' ';

    std::map<std::string, std::pair<std::string, int> > m_key_map;

    Configurator::OP_state isValidOption(const std::string& option,
                                         std::string& key,
                                         std::string& value);
    bool isValidKey(const std::string& key,
                    const std::string& sys_options = g_sysopts,
                    const char syskey_delim = Configurator::va_delim);

    template<typename T>
    bool isValidValue(const std::string& value,
                      const char va_delim,
                      const int num_of_values_required,
                      std::vector<T>& t_buff);

private:
    void setDefaultParms(SYSTEMPARAMETERS& sysparms);
    bool fillSysParms(const std::string& key,
                      const std::string& value,
                      SYSTEMPARAMETERS& sysparms);
    void eraseCharFromString(std::string& s, const char c);
    void eraseCharFromString(std::string& s, const std::string& clist);
    void separateString(const std::string& s,
                        std::vector<std::string>& s_buff,
                        const char delim = char_space);
};

}; // end of UCDPETSYS

#include "../src/details/petsys_cfg.impl"

#endif // petsys_cfg.h
