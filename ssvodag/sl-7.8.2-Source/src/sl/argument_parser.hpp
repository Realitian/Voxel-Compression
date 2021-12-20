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
#ifndef SL_ARGUMENT_PARSER_HPP
#define SL_ARGUMENT_PARSER_HPP

#include <sl/any.hpp>
#include <sl/any_extractor.hpp>
#include <list> 
#include <string> 
#include <sstream>

namespace sl {

  class argument_record {
  protected:
    std::string     tag_string_;
    std::string     name_string_;
    any_extractor   *extractor_;
    any             *value_pointer_;
    any             default_value_;
    std::string     description_;
  protected:

    bool            last_operation_success_;
    std::string     last_error_string_;
    
  public:

    argument_record(const char      *tag_string,
		    const char      *name_string,
		    any_extractor   *extractor,
		    any             *value_pointer,
		    const char      *description);

    argument_record(const std::string& tag_string,
		    const std::string& name_string,
		    any_extractor      *extractor,
		    any                *value_pointer,
		    const std::string& description);

    ~argument_record();

    void reset_defaults();

    bool is_optional() const;

    bool is_tagged() const;

    bool is_literal() const;

    void extract_front_arguments(std::istream& istr);

    bool last_operation_success() const;

    const std::string& tag_string() const;

    const std::string& last_error_string() const;
    
    const std::string& name_string() const;

    std::string syntax_description() const;

    const std::string& description() const;

    std::string glossary_description() const;
  };

  
  class argument_parser {
  protected:
    std::string     synopsis_;
    std::string     glossary_;

    std::list<argument_record> argument_syntax_;
    bool            is_accepting_extra_arguments_;

  protected:
    bool            last_operation_success_;
    std::string     last_error_string_;
    std::list<std::string> last_extra_arguments_;
  
  public:
    
    argument_parser();

    argument_parser(const std::list<argument_record>& arg_syntax,
		    bool accept_extra = false);

    argument_parser(argument_record* arg_syntax,
		    std::size_t      arg_syntax_count,
		    bool accept_extra = false);

    virtual ~argument_parser();

  public:

    void clear_syntax();

    void append_syntax(const std::list<argument_record>& arg_syntax);

    void append_syntax(argument_record* arg_syntax,
		       std::size_t      arg_syntax_count);
    
    void accept_extra_arguments(bool b);
    
    bool is_accepting_extra_arguments() const;
    
  public:
    
    const std::string& synopsis() const;

    const std::string& glossary() const;
    
  public:

    void reset_defaults();
		    
    void parse(std::size_t argc,
	       const char* argv[]);

    void parse(std::istream& istr);

    void parse(const std::string& str);

    bool last_operation_success() const;

    const std::string& last_error_string() const;

    const std::list<std::string>& last_extra_arguments() const;

  protected:
    
    void append_last_error(const std::string& s);

    void append_last_error(const argument_record& a);

    void rebuild_synopsis();

    void rebuild_glossary();
  };

}

#endif


