//+++HDR+++
//======================================================================
//  This file is part of the SL software library.
//
//  Copyright (C) 1993-2014 by Enrico Gobbetti (gobbetti@crs4.it)
//  Copyright (C) 1996-2014 by CRS4 Visual Computing Group, Pula, Italy
//
//  For more information, visit the CRS4 Visual Computing Group 
//  web pages at http://www.crs4.it/vvr/.
//
//  This file may be used under the terms of the GNU General Public
//  License as published by the Free Software Foundation and appearing
//  in the file LICENSE included in the packaging of this file.
//
//  CRS4 reserves all rights not expressly granted herein.
//  
//  This file is provided AS IS with NO WARRANTY OF ANY KIND, 
//  INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS 
//  FOR A PARTICULAR PURPOSE.
//
//======================================================================
//---HDR---//

#include <sl/utility.hpp>

 
#ifdef _WIN32
// Avoid sprintf warnings
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#endif

namespace sl {

  /// The characters used for directory separation in path names
  std::string pathname_directory_separators() {
    return "/\\"; // FIXME: Includes Windows style backslash...
  }

  /// The characters used for extension separation in path names
  std::string pathname_extension_separators() {
    return "."; 
  }

  /// Is c a directory separator?
  bool pathname_is_directory_separator(const char& c) {
    return pathname_directory_separators().find(c) != std::string::npos;
  }

  /// Is s an absolute path?
  bool pathname_is_absolute(const std::string& s) {
    std::size_t len = s.length();
 
    if(len == 0) return false;
    if(pathname_is_directory_separator(s[0])) return true;

    // FIXME The following is for windows style: c:/...
    if((len == 2 && isalpha(s[0]) && s[1] == ':') ||
       (len > 2 && isalpha(s[0]) && s[1] == ':' && pathname_is_directory_separator(s[2]))) {
      return true;
    }
    return false;
  }

  /// The path name before last directory separator
  std::string pathname_directory(const std::string& s) {
    std::size_t pos = s.find_last_of(pathname_directory_separators());
    if (pos == std::string::npos) return "."; else return s.substr(0,pos);
  }

  /// The path name after last directory separator
  std::string pathname_base(const std::string& s) {
    std::size_t pos = s.find_last_of(pathname_directory_separators());
    if (pos == std::string::npos) return s; else return s.substr(pos+1);
  }

  /// The base name after the last extension separator
  std::string pathname_extension(const std::string& s) {
    std::string base = pathname_base(s);
    std::size_t pos = base.find_last_of(pathname_extension_separators());
    if (pos == std::string::npos) {
      return "";
    } else {
      return base.substr(pos+1);
    }
  }

  /// The path name before the last '.'
  std::string pathname_without_extension(const std::string& s) {
    std::size_t dpos = s.find_last_of(pathname_directory_separators());
    std::size_t epos = s.find_last_of(pathname_extension_separators());
    if (epos == std::string::npos) {
      // No extension separator
      return s;
    } else if ((dpos != std::string::npos) && dpos > epos) {
      // Extension separator in directory section
      return s;
    } else {
      return s.substr(0,epos);
    }
  }

  
  /**
   * Matches a string with a pattern. Special characters:
   *	    '*'  Matches any string, including the null string.
   *        '?'  Matches any single character.
   *        '\'  Escapes next special character
   */
  bool matches(const char* str, const char* regexp) {
    if (str[0] == '\0') {
      if (regexp[0] == '\0') {
        return true;
      } else if (regexp[0] == '*') {
        return matches(str, regexp+1);
      } else {
        return false;
      }
    } else if (regexp[0] == '\\') {
      return regexp[1] == str[0] && matches(str+1, regexp+2);
    } else if (regexp[0] == '*') {
      return matches(str+1, regexp) || matches(str, regexp+1);
    } else if (regexp[0] == '\0') {
      return false;
    } else if (regexp[0] == str[0] || regexp[0] == '?') {
      return matches(str+1, regexp+1);
    } else {
      return false;
    }
  }
  
  /**
   * Matches a string with a pattern. Special characters:
   *	    '*'  Matches any string, including the null string.
   *        '?'  Matches any single character.
   *        '\'  Escapes next special character
   */
  bool matches(const std::string& str, const std::string& regexp) {
    return matches(str.c_str(), regexp.c_str());
  }
  
  std::string human_readable_duration(const time_duration& t) {
    static const uint64_t ms_usec = uint64_t(1000);
    static const uint64_t s_usec  = uint64_t(1000)*ms_usec;
    static const uint64_t m_usec  = uint64_t(  60)*s_usec;
    static const uint64_t h_usec  = uint64_t(  60)*m_usec;
    
    std::string result;
    uint64_t usec = t.as_microseconds();
    
    uint32_t  h  = uint32_t(usec / h_usec);  usec -=  uint64_t(h) * h_usec;
    uint32_t  m  = uint32_t(usec / m_usec);  usec -=  uint64_t(m) * m_usec;
    uint32_t  s  = uint32_t(usec / s_usec);  usec -=  uint64_t(s) * s_usec;
    uint32_t  ms = uint32_t(usec / ms_usec); usec -=  uint64_t(ms) * ms_usec;
    
    if (h>0) {
      result = to_string(h) + "h" + to_string(m) + "m";
    } else if (m>0) {
      result = to_string(m) + "m" + to_string(s) + "s";
    } else if (s>0) {
      result = to_string(s) + "s";
    } else if (ms>0) {
      result = to_string(ms) + "ms";
    } else {
      result = to_string(usec) + "us";
    }
    return result;
  }

  std::string human_readable_quantity(const uint64_t q,
				      const char*    unit,
				      const uint64_t K) {
    const uint64_t M  = uint64_t(1000)*K;
    const uint64_t G  = uint64_t(1000)*M;

    char cstr[64]; 
    if (q>=G) {
      sprintf(cstr, "%.1fG%s", double(q)/double(G), unit);
    } else if (q>=M) {
      sprintf(cstr, "%.1fM%s", double(q)/double(M), unit);
    } else if (q>=K) {
      sprintf(cstr, "%.1fK%s", double(q)/double(K), unit);
    } else {
      sprintf(cstr, "%d%s", int(q), unit);
    }
    return std::string(cstr);
  }
  
  std::string human_readable_size(const uint64_t sz) {
    return human_readable_quantity(sz, "B", 1024);
  }

  std::string human_readable_percent(const double percent) {
    char cstr[64]; 
    if (percent>=10.0f) {
      sprintf(cstr, "%4.0f%%", percent);
    } else if (percent>=1.0f) {
      sprintf(cstr, "%4.1f%%", percent);
    } else {
      sprintf(cstr, "%4.2f%%", percent);
    }
    return std::string(cstr);
  }

  /// The i-th prime number (0<=i<=255)
  std::size_t small_prime(const unsigned char i) {
    static const std::size_t ptable[256] = {
      2,      3,      5,     7,    11,    13,    17,    19,    23,    29, 
      31,     37,     41,     43,   47,    53,    59,    61,    67,    71, 
      73,     79,     83,     89,   97,   101,   103,   107,   109,   113, 
      127,    131,    137,    139,  149,   151,   157,   163,   167,   173, 
      179,    181,    191,    193,  197,   199,   211,   223,   227,   229, 
      233,    239,    241,    251,  257,   263,   269,   271,   277,   281, 
      283,    293,    307,    311,  313,   317,   331,   337,   347,   349, 
      353,    359,    367,   373,   379,   383,   389,   397,   401,   409, 
      419,    421,    431,   433,   439,   443,   449,   457,   461,   463, 
      467,    479,    487,   491,   499,   503,   509,   521,   523,   541, 
      547,    557,    563,   569,   571,   577,   587,   593,   599,   601, 
      607,    613,    617,   619,   631,   641,   643,   647,   653,   659, 
      661,    673,    677,   683,   691,   701,   709,   719,   727,   733, 
      739,    743,    751,   757,   761,   769,   773,   787,   797,   809,
      811,    821,    823,   827,   829,   839,   853,   857,   859,   863, 
      877,    881,    883,   887,   907,   911,   919,   929,   937,   941, 
      947,    953,    967,   971,   977,   983,   991,   997,   1009,   1013, 
      1019,   1021,   1031,  1033,   1039,  1049,  1051,  1061,  1063,  1069, 
      1087,   1091,   1093,  1097,   1103,  1109,  1117,  1123,  1129,  1151, 
      1153,   1163,   1171,  1181,   1187,  1193,  1201,  1213,  1217,  1223, 
      1229,   1231,   1237,  1249,   1259,  1277,  1279,  1283,  1289,  1291, 
      1297,   1301,   1303,  1307,   1319,  1321,  1327,  1361,  1367,  1373, 
      1381,   1399,   1409,  1423,   1427,  1429,  1433,  1439,  1447,  1451, 
      1453,   1459,   1471,  1481,   1483,  1487,  
    };
    return ptable[i];
  }

}
