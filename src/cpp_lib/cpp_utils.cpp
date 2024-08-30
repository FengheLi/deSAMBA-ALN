/*
 * cpp_utils.cpp
 *
 *  Created on: 2021年10月27日
 *      Author: fenghe
 */

#include "cpp_utils.hpp"
#include <cstring>

void split_string(std::vector<std::string> &item_value, const char * split_line, const char *split_str){
	item_value.clear();
	std::string s(split_line);
	std::string separator(split_str);
	int t = 1, pre = 0;
	for(int item_idx = 1; item_idx < 10000000; item_idx++){
        t = s.find(separator, pre);
        if (t != -1)
        	item_value.push_back(s.substr(pre, t - pre));
        else
            break;
        pre = t + separator.size();
	}
    if (pre < s.size())
    	item_value.push_back(s.substr(pre));
}
