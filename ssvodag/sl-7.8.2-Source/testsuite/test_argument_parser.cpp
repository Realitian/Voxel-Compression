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
/////// ALWAYS TEST IN DEBUG MODE
#if !(defined(__sgi) && !defined(__GNUC__))
#  undef NDEBUG
#endif
///////

#include <sl/tester.hpp>
#include <sl/argument_parser.hpp>

//--------------------

static std::size_t failed_test_count = 0;

void test1() {
  sl::any arg_x;
  sl::any arg_y(1.0f);
  sl::any arg_z;
  sl::any arg_w(true);
  sl::any arg_t;
  sl::any arg_s(std::string("NOSTRING"));
  sl::any arg_e(std::size_t(1));

  sl::argument_record argtable [] = 
  {
    //-------------------tag-----name-----extractor-------------------------------------var--------description
    sl::argument_record(NULL,    "x",     new sl::generic_extractor<float>,             &arg_x,    "untagged required argument"),
    sl::argument_record(NULL,    "y",     new sl::generic_extractor<float>,             &arg_y,    "untagged optional argument"),
    sl::argument_record("-z",    "bool",  new sl::generic_extractor<bool>,              &arg_z,    "tagged required boolean argument"),
    sl::argument_record("-w",    "bool",  new sl::generic_extractor<bool>,              &arg_w,    "tagged optional boolean argument"),
    sl::argument_record("-t",    NULL,    new sl::generic_extractor<void>,              &arg_t,    "tagged required literal argument"),
    sl::argument_record("-s",    "string",new sl::generic_extractor<std::string>,       &arg_s,    "tagged optional string argument"),
    sl::argument_record("-e",    NULL,    new sl::enum_extractor("red","green","blue"), &arg_e,    "tagged optional enum argument")
  };
  std::size_t narg = sizeof(argtable)/sizeof(sl::argument_record);

  std::string expected_synopsis =
    "<x> [ <y> ] -z <bool> [ -w <bool> ] [ -t ] [ -s <string> ] [ -e <red|green|blue> ] [...]";
  std::string expected_glossary =
    "<x>\n"
    "untagged required argument\n"
    "\n"
    "[ <y> ] (default = 1)\n"
    "untagged optional argument\n"
    "\n"
    "-z <bool>\n"
    "tagged required boolean argument\n"
    "\n"
    "[ -w <bool> ] (default = 1)\n"
    "tagged optional boolean argument\n"
    "\n"
    "[ -t ]\n"
    "tagged required literal argument\n"
    "\n"
    "[ -s <string> ] (default = NOSTRING)\n"
    "tagged optional string argument\n"
    "\n"
    "[ -e <red|green|blue> ] (default = green)\n"
    "tagged optional enum argument\n"
    ;

  sl::tester tester("argument_parser");

  sl::argument_parser scanarg;

  scanarg.append_syntax(argtable, narg);
  scanarg.accept_extra_arguments(true);
  
  tester.test("Append syntax",  scanarg.last_operation_success(), true);
  tester.test("Append syntax",  scanarg.last_error_string(), std::string());
  
#if 0
  tester.test("Synopsis", scanarg.synopsis(), expected_synopsis);
#endif
  tester.test("Glossary", scanarg.glossary(), expected_glossary);

  arg_s = (std::string("NOSTRING"));
  scanarg.reset_defaults();
  scanarg.parse(std::string("1.0 2.0 -z true -w false -t"));
  tester.test("Parse-1-success",  scanarg.last_operation_success(), true);
  tester.test("Parse-1-string",  scanarg.last_error_string(), std::string());
  tester.test("Parse-1-x",  arg_x.to_string(), std::string("1"));
  tester.test("Parse-1-y",  arg_y.to_string(), std::string("2"));
  tester.test("Parse-1-z",  arg_z.to_string(), std::string("1"));
  tester.test("Parse-1-w",  arg_w.to_string(), std::string("0"));
  tester.test("Parse-1-t",  arg_t.to_string(), std::string("1"));
  tester.test("Parse-1-s",  arg_s.to_string(), std::string("NOSTRING"));
  tester.test("Parse-1-e",  arg_e.to_string(), std::string("1"));
  tester.test("Parse-1-extra",  scanarg.last_extra_arguments().size(), std::size_t(0));

  arg_s = (std::string("NOSTRING"));
  scanarg.reset_defaults();
  scanarg.parse(std::string("1.0 -z true -t"));
  tester.test("Parse-2-success",  scanarg.last_operation_success(), true);
  tester.test("Parse-2-string",  scanarg.last_error_string(), std::string());
  tester.test("Parse-2-x",  arg_x.to_string(), std::string("1"));
  tester.test("Parse-2-y",  arg_y.to_string(), std::string("2"));
  tester.test("Parse-2-z",  arg_z.to_string(), std::string("1"));
  tester.test("Parse-2-w",  arg_w.to_string(), std::string("0"));
  tester.test("Parse-2-t",  arg_t.to_string(), std::string("1"));
  tester.test("Parse-2-s",  arg_s.to_string(), std::string("NOSTRING"));
  tester.test("Parse-2-e",  arg_e.to_string(), std::string("1"));
  tester.test("Parse-2-extra",  scanarg.last_extra_arguments().size(), std::size_t(0));

  arg_s = (std::string("NOSTRING"));
  scanarg.reset_defaults();
  scanarg.parse(std::string("1.0 -z true -s 'escaped string' -t"));
  tester.test("Parse-3-success",  scanarg.last_operation_success(), true);
  tester.test("Parse-3-string",  scanarg.last_error_string(), std::string());
  tester.test("Parse-3-x",  arg_x.to_string(), std::string("1"));
  tester.test("Parse-3-y",  arg_y.to_string(), std::string("2"));
  tester.test("Parse-3-z",  arg_z.to_string(), std::string("1"));
  tester.test("Parse-3-w",  arg_w.to_string(), std::string("0"));
  tester.test("Parse-3-t",  arg_t.to_string(), std::string("1"));
  tester.test("Parse-3-s",  arg_s.to_string(), std::string("escaped string"));
  tester.test("Parse-3-e",  arg_e.to_string(), std::string("1"));
  tester.test("Parse-3-extra",  scanarg.last_extra_arguments().size(), std::size_t(0));

  arg_s = (std::string("NOSTRING"));
  scanarg.reset_defaults();
  scanarg.parse(std::string("1.0 2.0 -z true -w false -t extra -e blue"));
  tester.test("Parse-4-success",  scanarg.last_operation_success(), true);
  tester.test("Parse-4-string",  scanarg.last_error_string(), std::string());
  tester.test("Parse-4-x",  arg_x.to_string(), std::string("1"));
  tester.test("Parse-4-y",  arg_y.to_string(), std::string("2"));
  tester.test("Parse-4-z",  arg_z.to_string(), std::string("1"));
  tester.test("Parse-4-w",  arg_w.to_string(), std::string("0"));
  tester.test("Parse-4-t",  arg_t.to_string(), std::string("1"));
  tester.test("Parse-4-s",  arg_s.to_string(), std::string("NOSTRING"));
  tester.test("Parse-3-e",  arg_e.to_string(), std::string("2"));
  tester.test("Parse-4-extra",  scanarg.last_extra_arguments().size(), std::size_t(1));
  
  failed_test_count += (std::size_t)(tester.failed_test_count());
}

int main() {

  test1();

  return (int)failed_test_count;

}


