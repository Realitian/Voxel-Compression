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
#include <sl/argument_parser.hpp>
#include <sstream> // for conversion to/from strings
#include <set>
#include <map>

namespace sl {
  
  argument_record::argument_record(const std::string& tag_string,
				   const std::string& name_string,
				   any_extractor      *extractor,
				   any                *value_pointer,
				   const std::string& description)
    :
    tag_string_(tag_string),
    name_string_(name_string),
    extractor_(extractor),
    value_pointer_(value_pointer),
    default_value_(value_pointer ? *value_pointer : any()),
    description_(description),
    last_operation_success_(true),
    last_error_string_(std::string("")) {

    SL_REQUIRE("Value pointer exists", value_pointer);
    SL_REQUIRE("Extractor exists", extractor);

    if (name_string_.empty()) {
      name_string_ = extractor_->type_name();
    }
  }

  argument_record::argument_record(const char      *tag_string,
				   const char      *name_string,
				   any_extractor   *extractor,
				   any             *value_pointer,
				   const char      *description)
    :
    tag_string_(std::string(tag_string ? tag_string : "")),
    name_string_(std::string(name_string ? name_string : "")),
    extractor_(extractor),
    value_pointer_(value_pointer),
    default_value_(value_pointer ? *value_pointer : any()),
    description_(std::string(description ? description : "")),
    last_operation_success_(true),
    last_error_string_(std::string("")) {

    SL_REQUIRE("Value pointer exists", value_pointer);
    SL_REQUIRE("Extractor exists", extractor);

    if (name_string_.empty()) {
      name_string_ = extractor_->type_name();
    }
  }

  argument_record::~argument_record() {
  }

  void argument_record::reset_defaults() {
    last_operation_success_ = true;
    last_error_string_ = "";
    if (is_literal()) {
      *value_pointer_ = any(); // Mark not found
    } else  if (!is_optional()) {
      *value_pointer_ = default_value_;
    }
  }
  
  bool argument_record::is_optional() const {
    return is_literal() || (!default_value_.empty());
  }

  bool argument_record::is_tagged() const {
    return !tag_string_.empty();
  }

  bool argument_record::is_literal() const {
    return dynamic_cast< generic_extractor<void>* >(extractor_) != NULL;
  }

  const std::string& argument_record::tag_string() const {
    return tag_string_;
  }
  
  void argument_record::extract_front_arguments(std::istream& istr) {
    last_operation_success_ = true;
    last_error_string_ = std::string("");

    if (is_literal()) {
      *value_pointer_ = any(true); // Mark non-empty value
    } else {
      extractor_->extract(istr, *value_pointer_);
      if (istr.fail()) {
	last_operation_success_ = false;
	last_error_string_ = std::string("Bad value for argument of ") +  syntax_description();
      }
    }
  }
  
  bool argument_record::last_operation_success() const {
    return last_operation_success_;
  }

  const std::string& argument_record::last_error_string() const {
    return last_error_string_;
  }

  const std::string& argument_record::name_string() const {
    return name_string_;
  }

  std::string argument_record::syntax_description() const {
    std::string result;

    if (is_optional()) result += std::string("[ ");
    if (is_tagged()) {
      result += tag_string_;
      if (!is_literal()) result += std::string(" ");
    }
    if (!is_literal()) result += std::string("<") + name_string() + std::string(">");
    if (is_optional()) result += std::string(" ]");
    return result;
  }
  
  const std::string& argument_record::description() const {
    return description_;
  }
  
  std::string argument_record::glossary_description() const {
    std::string result;
    if (!description().empty()) {
      result += syntax_description();
      if (!is_literal()) {
	if (is_optional()) {
	  result += std::string(" (default = ") + extractor_->printable_representation(default_value_) + std::string(")");
	}
      }
      result += std::string("\n");
      result += description();
      result += std::string("\n");
    }
    return result;
  }

  ////////////////////////////////////////////////////////////////////

  argument_parser::argument_parser(const std::list<argument_record>& arg_syntax,
				   bool accept_extra) {
    last_operation_success_ = true; 
    clear_syntax();
    append_syntax(arg_syntax);
    accept_extra_arguments(accept_extra);
  }

  argument_parser::argument_parser(argument_record* arg_syntax,
				   std::size_t      arg_syntax_count,
				   bool accept_extra) {
    last_operation_success_ = true; 
    clear_syntax();
    append_syntax(arg_syntax,
		  arg_syntax_count);
    accept_extra_arguments(accept_extra);
  }

  argument_parser::argument_parser() {
    last_operation_success_ = true; 
    clear_syntax();
  }

  argument_parser::~argument_parser() {
  }

  void argument_parser::clear_syntax() {
    argument_syntax_.clear();
    last_operation_success_ = true; 
    last_error_string_ = "";
    last_extra_arguments_.clear();
    synopsis_ = "";
    glossary_ = "";
  }

  void argument_parser::append_syntax(const std::list<argument_record>& arg_syntax) {
    for (std::list<argument_record>::const_iterator it = arg_syntax.begin();
	 it != arg_syntax.end();
	 ++it) {
      argument_syntax_.push_back(*it);
    }
    last_operation_success_ = true; 
    last_error_string_ = "";
    last_extra_arguments_.clear();
    rebuild_synopsis();
    rebuild_glossary();
  }
  
  void argument_parser::append_syntax(argument_record* arg_syntax,
				      std::size_t      arg_syntax_count) {
    std::list<argument_record> syntax_as_list;
    for (std::size_t i=0; i<arg_syntax_count; ++i) {
      syntax_as_list.push_back(arg_syntax[i]);
    }
    append_syntax(syntax_as_list);
  }

  void argument_parser::accept_extra_arguments(bool b) {
    is_accepting_extra_arguments_ = b;
  }

  bool argument_parser::is_accepting_extra_arguments() const {
    return is_accepting_extra_arguments_;
  }
  
  const std::string& argument_parser::synopsis() const {
    return synopsis_;
  }
  
  const std::string& argument_parser::glossary() const {
    return glossary_;
  }
  
  void argument_parser::reset_defaults() {
    last_operation_success_ = true;
    last_error_string_ = "";
    for (std::list<argument_record>::iterator it = argument_syntax_.begin();
	 it != argument_syntax_.end();
	 ++it) {
      (*it).reset_defaults();
      bool ok = (*it).last_operation_success();
      if (!ok) {
	last_operation_success_ = false;
	if (!last_error_string_.empty()) {
	  last_error_string_ += std::string("\n");
	}
	last_error_string_ += (*it).last_error_string();
      }
    }
  }
		    
  void argument_parser::parse(std::size_t argc,
			      const char* argv[]) {
    std::stringstream iostr;
    for (std::size_t i=1; i<argc; ++i) {
      if (i != 1) iostr << " ";
      std::string s(argv[i]);
      if (s.find(' ')  != std::string::npos ||
	  s.find('\t') != std::string::npos ||
	  s.find('\b') != std::string::npos ||
	  s.find('\n') != std::string::npos ||
	  s.find('\r') != std::string::npos ||
	  s.find('\'') != std::string::npos ||
	  s.find('"')  != std::string::npos) {
	iostr << escape(s);
      } else {
	iostr << s;
      }
    }
    parse(iostr);
  }

  void argument_parser::parse(const std::string& str) {
   std::stringstream iostr;
   iostr << str;
   parse(iostr);
  }

  bool argument_parser::last_operation_success() const {
    return last_operation_success_;
  }
  
  const std::string&  argument_parser::last_error_string() const {
    return last_error_string_;
  }

  const std::list<std::string>&  argument_parser::last_extra_arguments() const {
    return last_extra_arguments_;
  }
  
  void argument_parser::parse(std::istream& istr) {
    last_operation_success_ = true;
    last_error_string_ = "";
    last_extra_arguments_.clear();

    std::set<const argument_record*>        unprocessed_required_arguments;
    std::list<argument_record*>             untagged_arguments;
    std::map<std::string, argument_record*> tagged_arguments;

    for (std::list<argument_record>::iterator it = argument_syntax_.begin();
	 it != argument_syntax_.end();
	 ++it) {  
      if (!(*it).is_optional()) {
	unprocessed_required_arguments.insert(&(*it));
      }
      if ((*it).is_tagged()) {
	tagged_arguments[(*it).tag_string()] = &(*it);
      } else {
	untagged_arguments.push_back(&(*it));
      }
    }
    
    generic_extractor<std::string> strextr;
    token_extractor                tkextr;

    bool finished = false;
    bool aborted = false;

    while (!finished && !aborted) {
      // Skip blanks
      while (is_blank(istr.peek())) { istr.get(); }
      if (istr.peek() == EOF) istr.get();

      finished = istr.eof();
      aborted  = !finished && istr.fail();  // (!istr.ipfx(0));
      if (aborted) {
	last_operation_success_ = false;
	if (!last_error_string_.empty()) last_error_string_ += std::string("\n");
	last_error_string_ += "Error on input stream";
      }

      if (!finished && !aborted) {
	// Try a tagged argument
	std::streampos spos = istr.tellg();
	any tmp;
	tkextr.extract(istr, tmp);
	finished = istr.fail() && istr.eof();
	if (!finished) {
	  if (tagged_arguments.find(tmp.to_string()) != tagged_arguments.end()) {
	    // Handle tagged argument
	    argument_record* a = tagged_arguments[tmp.to_string()];
	    unprocessed_required_arguments.erase(a);
	    a->extract_front_arguments(istr);
	    append_last_error(*a);
	    aborted = !last_operation_success_;
	  } else {
	    // Backtrack and look for untagged or extra arguments
	    istr.seekg(spos);
	    istr.clear();

	    if (!untagged_arguments.empty()) {
	      argument_record* a = untagged_arguments.front();
	      untagged_arguments.pop_front();
	      unprocessed_required_arguments.erase(a);
	      a->extract_front_arguments(istr);
	      append_last_error(*a);
	      aborted = !last_operation_success_;
	    } else if (is_accepting_extra_arguments()) {
	      any tmp;
	      strextr.extract(istr, tmp);
	      if (!istr.fail()) {
		last_extra_arguments_.push_back(tmp.to_string());
	      } else {
		last_operation_success_ = false;
		if (!last_error_string_.empty()) last_error_string_ += std::string("\n");
		last_error_string_ += "Error while reading extra argument";
		aborted = !last_operation_success_;
	      }	 
	    } else {
	      last_operation_success_ = false;
	      if (!last_error_string_.empty()) last_error_string_ += std::string("\n");
	      last_error_string_ += "Unacceptable extra argument";
	      aborted = !last_operation_success_;
	    }   
	  }
	}
      }
    }

    // Check whether all required arguments where processed
    if (last_operation_success_) {
      for (std::set<const argument_record*>::iterator it = unprocessed_required_arguments.begin();
	   it != unprocessed_required_arguments.end();
	   ++it) {
	last_operation_success_ = false;
	if (!last_error_string_.empty()) last_error_string_ += std::string("\n");
	last_error_string_ += "Missing required argument " + (*it)->syntax_description();
      }
    }

  }
  
  void argument_parser::append_last_error(const argument_record& a) {
    if (!a.last_operation_success()) {
      last_operation_success_ = false;
      if (!last_error_string_.empty()) last_error_string_ += std::string("\n");
      last_error_string_ += a.last_error_string();
    }
  }

  void argument_parser::rebuild_synopsis() {
    synopsis_ = "";
    for (std::list<argument_record>::iterator it = argument_syntax_.begin();
	 it != argument_syntax_.end();
	 ++it) {  
      if (!synopsis_.empty()) {
	synopsis_ += std::string(" ");
      }
      synopsis_ += (*it).syntax_description();
    }
    if (is_accepting_extra_arguments()) {
      if (!synopsis_.empty()) {
	synopsis_ += std::string(" ");
      }
      synopsis_ += "[...]";
    }
  }

  void argument_parser::rebuild_glossary() {
    glossary_ = "";
    for (std::list<argument_record>::iterator it = argument_syntax_.begin();
	 it != argument_syntax_.end();
	 ++it) {  
      if (!glossary_.empty()) {
	glossary_ += std::string("\n");
      }
      glossary_ += (*it).glossary_description();
    }
  }

}
