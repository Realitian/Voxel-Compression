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
#ifndef SL_ANY_EXTRACTOR_HPP
#define SL_ANY_EXTRACTOR_HPP

#include <sl/any.hpp>
#include <string> 
#include <iostream>
#include <iomanip>
#include <stack>
#include <list>
#include <map>

namespace sl {

  /// Test if a char is a blanks - move somewhere else...
  inline bool is_blank(int c) {
    return ((c == ' ' ) || 
	    (c == '\t') || 
	    (c == '\b') || 
	    (c == '\n') || 
	    (c == '\r'));
  }

  /// Advance in the stream until istr.peek() is not a blank
  inline void skip_blanks(std::istream& istr) {
    while (is_blank(istr.peek())) { istr.get(); };
  }

  /**
   *  This is an iostream manipulator similar to width() or endl
   *  See April 2000 C++ Report, "Effective Standard C++ Library"
   *  for full description of manipulators and how to implement them.
   *
   *  The manipulator escapes a string, surrounding it with
   *  a given char (if non null, default double quotes) and prefixing 
   *  the surround character as well as all characters in a given 
   *  set with an escape character (backslash by default).
   */
  class escape {
  protected:
    const char *input_string_;
    const char *escape_string_;
    char escape_char_;
    char surround_char_;
  public:
    explicit escape(const std::string &input_string,
		    char escape_char='\\',
		    char surround_char='"',
		    const char *chars_to_escape=NULL)
      :
      input_string_(input_string.c_str()),
      escape_string_(chars_to_escape),
      escape_char_(escape_char),
      surround_char_(surround_char)
    {
    }

    explicit escape(const char *input_string,
		    char escape_char='\\',
		    char surround_char='"',
		    const char *chars_to_escape=NULL)
      :
      input_string_(input_string),
      escape_string_(chars_to_escape),
      escape_char_(escape_char),
      surround_char_(surround_char)
    {
    }

    std::ostream &do_escape(std::ostream &os) const {
      if (surround_char_ != '\0') {
	os.put(surround_char_);
      }
      if (input_string_) {
	int charPos=0;
	while (input_string_[charPos]) {
	  bool escapeIt=false;
	  if (input_string_[charPos] == escape_char_ ||
	      input_string_[charPos] == surround_char_) {
	    escapeIt=true;
	  } else if(escape_string_ != NULL ) {
	    // Possible enhancement: use strchr here instead
	    // of the for loop.
	    for (int srchPos=0; escape_string_[srchPos]; srchPos++) {
	      if (input_string_[charPos] == escape_string_[srchPos]) {
		escapeIt=true;
		break;
	      }
	    }
	  }
	  if (escapeIt) {
	    os.put(escape_char_);
	  }
	  os.put(input_string_[charPos++]);
	}
      }
      if (surround_char_ != '\0') {
	os.put(surround_char_);
      }
      return os;
    }
    
  };

}

inline std::ostream& operator << (std::ostream &os, const sl::escape &e) {
  return e.do_escape(os);
}

namespace sl {
  /**
   *  This is an iostream manipulator similar to width() or endl
   *  See April 2000 C++ Report, "Effective Standard C++ Library"
   *  for full description of manipulators and how to implement them.
   *
   *  The manipulator unescapes a string, removing surround and
   *  escape characters.
   */
  class unescape {
  protected:
    char escape_char_;
    char surround_char_;
    std::string &output_string_;
  public:
    explicit unescape(std::string &output_string,
		      char escape_char='\\',
		      char surround_char='"'
		      )
    :
    escape_char_(escape_char),
    surround_char_(surround_char),
    output_string_(output_string)
    {
    }
    std::istream &do_unescape(std::istream &is) {
      output_string_="";
      char nextChar ='\0';
      
      // Use extraction operator to skip white spaces
      is >> nextChar;
      if (nextChar == surround_char_) {
	nextChar=is.get();
      } 	
      
      // Keep reading characters until there is an error, we
      // hit the end of file, or we hit the terminator.
      while (is.good() && !is.eof()
	     && nextChar != surround_char_) {
	// simply read each character one at a time,
	// and un-escape them.
	if (nextChar == escape_char_) {
	  nextChar=is.get();
	  if (is.eof()) {
	    break;
	  }
	}
	output_string_ += nextChar;
	nextChar=is.get();
      }
      return is;
    }
  };
}

inline std::istream& operator >> (std::istream &is, sl::unescape &e) {
  return e.do_unescape(is);
}

namespace sl {

  /**
   * Objects that read values from an input stream. 
   */
  class any_extractor {
  protected:
    std::stack<std::streampos>          spos_;
    std::stack<std::ios::fmtflags>  sflags_;
  public:

    inline any_extractor() {
    }

    virtual ~any_extractor() {
    }

    /**
     *  Read x from is if possible. If anything goes wrong, 
     *  the function returns with failbit set and x unchanged.
     */
    virtual void extract(std::istream& istr, any& x) {
      extract_pre(istr);
#if 0
      if (istr.ipfx(0)) {
#endif
	extract_do(istr, x);
#if 0
      }
#endif
      extract_post(istr);
    }

    virtual std::string type_name() const {
      return std::string("???");
    }

    virtual std::string printable_representation(const any& x) const {
      return x.to_string();
    }

  protected:

    virtual void extract_pre(std::istream& istr) {
      spos_.push(istr.tellg());
      sflags_.push(istr.flags());
    }

    virtual void extract_post(std::istream& istr) {
      if (!istr.fail()) {
	istr.flags(sflags_.top());
      } else {
	istr.seekg(spos_.top());
	istr.flags(sflags_.top());
	istr.clear(std::ios::failbit | istr.rdstate());
      }
      spos_.pop();
      sflags_.pop();
    }

    virtual void extract_do(std::istream& istr, any& x) = 0;

  };

  template <class T>
  class generic_extractor: public any_extractor {
  public:

    generic_extractor() {}
    virtual ~generic_extractor() {}
    
    virtual std::string type_name() const { return typeid(T).name(); }

  protected:

    virtual void extract_do(std::istream& istr, any& x) {
      T tmp = T();
      istr >> tmp;
      if (!istr.fail()) {
	x = any(tmp);
      }
    }
  };

  
  template <>
  class generic_extractor<void>: public any_extractor {
  public:

    generic_extractor() {}
    virtual ~generic_extractor() {}
    
    virtual std::string type_name() const { return typeid(void).name(); }

  protected:

    virtual void extract_do(std::istream&, any& x) {
      x = any(true);
    }
  };

  template <>
  class generic_extractor<bool>: public any_extractor {
  public:

    generic_extractor() {}
    virtual ~generic_extractor() {}
    
    virtual std::string type_name() const { return typeid(bool).name(); }

  protected:

    virtual void extract_do(std::istream& istr, any& x) {
      bool tmp = false;
      istr.setf(std::ios::boolalpha);
      istr >> tmp;
      if (!istr.fail()) {
	x = any(tmp);
      }
    }
  };
    
  template <>
  class generic_extractor<std::string>: public any_extractor {
  public:

    generic_extractor() {}
    virtual ~generic_extractor() {}
    
    virtual std::string type_name() const { return typeid(std::string).name(); }

  protected:

    virtual void extract_do(std::istream& istr, any& x) {
      skip_blanks(istr);
      int nextchar = istr.peek();
      std::string tmp;
      if (nextchar == '"' || nextchar == '\'') {
	unescape(tmp,'\\',nextchar).do_unescape(istr); 	//istr >> unescape(tmp,'\\',nextchar);
      } else {
	istr >> tmp;
      }
      if (!istr.fail()) {
	x = any(tmp);
      }
    }
  };

  class token_extractor: public any_extractor {
  public:
    token_extractor() {}
    virtual ~token_extractor() {}
    
    virtual std::string type_name() const { return "token"; }

  protected:

    virtual void extract_do(std::istream& istr, any& x) {
      skip_blanks(istr);
      int c = istr.peek();
      
      std::string tmp;
      while (c != EOF && !is_blank(c)) {
	tmp += char(c);
	istr.get(); c = istr.peek();
      }
      if (!istr.fail()) {
	x = any(tmp);
      }
    }
    
  };

  class enum_extractor: public any_extractor {
  protected:
    std::map<std::string, std::size_t> label_to_code_;
    std::map<std::size_t, std::string> code_to_label_;
  public:
    enum_extractor(const std::map<std::string, std::size_t>& e) 
      : label_to_code_(e) {
      build_code_to_label_map();
    }

    enum_extractor(const char* labels[],
		   std::size_t labels_count) {
      for (std::size_t i=0; i<labels_count; ++i) {
	label_to_code_[std::string(labels[i])] = i;
      }
      build_code_to_label_map();
    }

    enum_extractor(const char* l0,
		    const char* l1 = NULL,
		    const char* l2 = NULL,
		    const char* l3 = NULL,
		    const char* l4 = NULL,
		    const char* l5 = NULL,
		    const char* l6 = NULL,
		    const char* l7 = NULL,
		    const char* l8 = NULL,
		    const char* l9 = NULL,
		    const char* l10 = NULL) {
      if (l0) label_to_code_[std::string(l0)] = std::size_t(0);
      if (l1) label_to_code_[std::string(l1)] = std::size_t(1);
      if (l2) label_to_code_[std::string(l2)] = std::size_t(2);
      if (l3) label_to_code_[std::string(l3)] = std::size_t(3);
      if (l4) label_to_code_[std::string(l4)] = std::size_t(4);
      if (l5) label_to_code_[std::string(l5)] = std::size_t(5);
      if (l6) label_to_code_[std::string(l6)] = std::size_t(6);
      if (l7) label_to_code_[std::string(l7)] = std::size_t(7);
      if (l8) label_to_code_[std::string(l9)] = std::size_t(8);
      if (l9) label_to_code_[std::string(l9)] = std::size_t(9);
      if (l10) label_to_code_[std::string(l10)] = std::size_t(10);
      build_code_to_label_map();
    }
      
      
    virtual ~enum_extractor() {}
    
    virtual std::string type_name() const { 
      std::string result;
      for (std::map<std::size_t, std::string>::const_iterator it = code_to_label_.begin();
	   it != code_to_label_.end();
	   ++it) {
	if (!result.empty()) { result += '|'; }
	result += (*it).second;
      }
      return result;
    }

    virtual std::string printable_representation(const any& x) const {
      const std::size_t* code_ptr = any_cast<std::size_t>(&x);
      if (code_ptr) {
	std::map<std::size_t, std::string>::const_iterator it =
	  code_to_label_.find(*code_ptr);
	if (it == code_to_label_.end()) {
	  return "INVALID ENUM: [" + x.to_string() + "]";
	} else {
	  return (*it).second;
	}
      } else {
	return "INVALID ENUM: [" + x.to_string() + "]";
      }
    }

  protected:

    void build_code_to_label_map() {
      code_to_label_.clear();
      for (std::map<std::string, std::size_t>::const_iterator it = label_to_code_.begin();
	   it != label_to_code_.end();
	   ++it) {
	code_to_label_[(*it).second] = (*it).first;
      }
    }
    
    virtual void extract_do(std::istream& istr, any& x) {
      skip_blanks(istr);
      int c = istr.peek();
      
      std::string tmp;
      while (c != EOF && !is_blank(c)) {
	tmp += char(c);
	istr.get(); c = istr.peek();
      }
      if (!istr.fail()) {
	if (label_to_code_.find(tmp) != label_to_code_.end()) {
	  x = any(label_to_code_[tmp]);
	} else {
	  istr.clear(std::ios::failbit | istr.rdstate());
	}
      }
    }

  };

} // namespace sl

#endif
