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
#ifndef SL_ANY_HPP
#define SL_ANY_HPP

#include <sl/utility.hpp>
#include <algorithm> // for std::swap
#include <typeinfo>

namespace sl {
  
  namespace detail {
    
    class placeholder {
    public:
      virtual inline ~placeholder() {
      }
    public:
      virtual const std::type_info & type() const = 0;
      virtual placeholder * clone() const = 0;
      virtual std::string to_string() const = 0;
    };

    template<class T>
    class holder: public placeholder {
    public:
      inline holder(const T & value): held_(value) {
      }
    public:

      virtual const std::type_info &type() const {
	return typeid(T);
      }
      virtual placeholder *clone() const {
	return new holder(held_);
      }
      virtual std::string to_string() const {
	std::ostringstream s;
	s << held_;
	return s.str();
      }
    public:
      T held_;
    };
  }

  /**
   *  Generic union type. Based on the class of the same name described in 
   *  "Valued Conversions" by Kevlin Henney, C++ Report 12(7), July/August 2000)
   */
  class any {
  public:
    detail::placeholder * content_;
  public:
    
    inline any() 
      : content_(0) {
    }

    template<class T>
    inline any(const T & value)
      : content_(new detail::holder<T>(value)) {
    }

    inline any(const any & other)
      : content_(other.content_ ? other.content_->clone() : 0) {
    }

    inline ~any() {
      if (content_) delete content_;
      content_ = 0;
    }

  public:

    inline any& swap(any & rhs) {
      std::swap(content_, rhs.content_);
      return *this;
    }

    template<class T>
    inline any& operator=(const T& rhs) {
      if (content_) delete content_;
      content_ = new detail::holder<T>(rhs);
      return *this;
    }

    inline any& operator=(const any& rhs) {
      if (content_) delete content_;
      content_ = (rhs.content_ ? rhs.content_->clone() : 0);
      return *this;
    }

  public:

    inline bool empty() const {
      return content_ == 0;
    }

    inline const std::type_info & type() const {
      return content_ ? content_->type() : typeid(void);
    }

    inline std::string to_string() const {
      return empty() ? std::string() : (content_)->to_string();
    }

  }; // class any

  class bad_any_cast : public std::bad_cast {
  public:
    virtual const char * what() const throw() {
      return "sl::bad_any_cast: failed conversion";
    }
  };
  
  template<class T>
  inline T* any_cast(any * operand) {
    return 
      (operand && (operand->type() == typeid(T))) 
      ? &(static_cast<detail::holder<T> *>(operand->content_)->held_)
      : 0;
  }
  
  template<class T>
  inline const T* any_cast(const any * operand) {
    return 
      (operand && (operand->type() == typeid(T))) 
      ? &(static_cast<const detail::holder<T> *>(operand->content_)->held_)
      : 0;
  }

  template<class T>
  inline T any_cast(const any & operand) {
    const T * result = any_cast<T>(&operand);
    if(!result)
      throw bad_any_cast();
    return *result;
  }
    
}

#endif
