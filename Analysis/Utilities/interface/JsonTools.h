#ifndef ICHiggsTauTau_Utilities_JsonTools_h
#define ICHiggsTauTau_Utilities_JsonTools_h

#include <vector>
#include <string>
#include "Utilities/interface/json.h"

namespace ic {
/**
 * Extracts a Json::Value from the specified input file
 *
 * @param file The full path to the input json file
 */
Json::Value ExtractJsonFromFile(std::string const& file);

/**
 * Extracts a Json::Value from the given input string
 *
 * @param str Must be a valid JSON object, including the enclosing curly
 * braces, e.g: `{ key : "value" }`.
 */
Json::Value ExtractJsonFromString(std::string const& str);

/**
 * Extracts a Json::Value from the given input string
 *
 * @param str should separate keys and required values by : and any
 * lists of values by ^ e.g: key:value1^value2.
 *
 * To pass a list consisting of only one value precede the value by ^
 * e.g.: key:^value1, for an empty list leave off the value arguments:
 * key:
 * 
 */
Json::Value ExtractJsonFromFlatString(std::string const& str);


/**
 * Updates the values in one json from the values in another
 *
 * All the values in *b* are entered into *a*, overwritting any already
 * defined in *a*.
 *
 * @param a The Json::Value to be updated \param b The Json::Value containing
 * the updates
 */
void UpdateJson(Json::Value& a, Json::Value const& b);

/**
 * Create a single merged Json::Value from a mixture of json files and
 * json-formatted strings
 *
 * Each string contained in `argv` must either be the name of a file
 * containing json data or a string of json-formatted text. A string will only
 * be assumed to indicate a file when it ends with the file extension `.json`.
 *
 * @note It is assumed that the arguments of this function are those passed to
 * a main function, and accordingly the first element of `argv` will be
 * ignored.
 *
 * @param argc Number of arguments contained in argv
 * @param argv An array of strings
 */
Json::Value MergedJson(int argc, char* argv[]);

/**
 * Create a single merged Json::Value from a mixture of json files and
 * json-formatted strings
 *
 * Each string contained in `vec` must either be the name of a file containing
 * json data or a string of json-formatted text. A string will only be assumed
 * to indicate a file when it ends with the file extension `.json`.
 *
 * @param vec Vector of input strings
 */
Json::Value MergedJson(std::vector<std::string> const& vec);
}
#endif
