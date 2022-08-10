#include "Utils.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

namespace bcm3 {

void tokenize(std::string str, std::vector<std::string>& tokens, const char* delims)
{
	tokens.clear();

	if (str.empty()) {
		return;
	}

	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(delims, NULL, boost::keep_empty_tokens);
	if (*str.rbegin() == '\n') {
		str.resize(str.size() - 1);
	}
	tokenizer tok(str, sep);
	tokenizer::const_iterator tiend = tok.end();
	for(tokenizer::const_iterator ti = tok.begin(); ti != tiend; ++ti) {
		tokens.push_back(*ti);
	}
}

void fix_path(std::string& path)
{
	if (path.empty()) {
		return;
	}
	for (std::string::iterator iter = path.begin(); iter != path.end(); ++iter) {
		if (*iter == '\\') {
			*iter = '/';
		}
	}
	if (*path.rbegin() != '/') {
		path += "/";
	}
}

}
