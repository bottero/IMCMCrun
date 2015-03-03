#include <iostream>
#include <map>
#include <string>
#define SEPARATOR '='
#define COMMENTS_CHAR "#"

namespace configFileParser
{
  //---------------------------------------------------------------------------
  // The configuration::data is a simple map string (key, value) pairs.
  // The file is stored as a simple listing of those pairs, one per line.
  // The key is separated from the value by an equal sign '='.
  // Commentaries are made with a hash ("#").
  //
  // Example:
  //   # This is an example
  //   source.directory = /a/path/to/a/directory
  //   file.types = *.jpg;*.gif;*.png;*.pix;*.tif;*.bmp
  //
  // Notice that the configuration file format does not permit values to span
  // more than one line or [section]s.
  //   
  struct data: std::map <std::string, std::string>
  {
    // Here is a little convenience method...
    bool iskey( const std::string& s ) const {
      return count( s ) != 0;
    }
  };
  //---------------------------------------------------------------------------
  // The extraction operator reads configuration::data until EOF.
  // Invalid data is ignored.
  //
  std::istream& operator >> ( std::istream& ins, data& d )
  {
    std::string s, key, value;

    // For each (key, value) pair in the file
    while (std::getline( ins, s )) {
      
      std::string::size_type begin = s.find_first_not_of( " \f\t\v" ); // ~int -> First non blank character 

      // Skip blank lines
      if (begin == std::string::npos) continue; // If the first non blank character is the end of the line

      // Skip commentary
      if (std::string( COMMENTS_CHAR ).find( s[ begin ] ) != std::string::npos) continue;

      // Extract the key value
      std::string::size_type end = s.find(SEPARATOR, begin );
      key = s.substr( begin, end - begin );

      // (No leading or trailing whitespace allowed)
      key.erase( key.find_last_not_of( " \f\t\v" ) + 1 );

      // No blank keys allowed
      if (key.empty()) continue;

      // Extract the value (no leading or trailing whitespace allowed)
      begin = s.find_first_not_of( " \f\n\r\t\v", end + 1 );
      if (s.find( "#") != std::string::npos) 
        end = s.find( "#");
      else
        end = s.find_last_not_of(  "\f\n\r\t\v" ) + 1;        

      value = s.substr( begin, end - begin );
        
      // No leading whitespace allowed : striping them
      std::string::size_type pos = value.find_last_not_of(' ');
      if (pos != value.length()-1) {
        if ( pos == std::string::npos )
          pos = -1;
        value.erase(pos+1);
      }
      // Insert the properly extracted (key, value) pair into the map
      d[ key ] = value;
      }

    return ins;
  }

  //---------------------------------------------------------------------------
  // The insertion operator writes all configuration::data to stream.
  //
  std::ostream& operator << ( std::ostream& outs, const data& d )
  {
    data::const_iterator iter;
    for (iter = d.begin(); iter != d.end(); iter++)
      outs << iter->first << " = " << iter->second << std::endl;
    return outs;
  }
} // namespace configuration 
