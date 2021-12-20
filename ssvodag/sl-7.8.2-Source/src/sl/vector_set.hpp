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
#ifndef SL_VECTOR_SET_HPP
#define SL_VECTOR_SET_HPP
#define SL_VERSION_VECTOR_SET_ 0x00010010

#include <algorithm>
#include <vector>
#include <utility>
#include <functional>
#include <sl/serializer.hpp>

namespace sl {

  /**
   * A STL set/multiset implemented on top of a (non-sorted) vector. Useful only for
   * very small sets.
   */
	template<class Key, bool bNoDuplicates=false, class EqualKey=std::equal_to<Key>, class Alloc=std::allocator<Key> >
	class vector_set_base {
  public:
    typedef vector_set_base<Key,bNoDuplicates,EqualKey,Alloc> this_t;
    typedef std::vector<Key,Alloc>		                        base_container_t;
    typedef typename base_container_t::allocator_type	        allocator_type;
    typedef typename base_container_t::size_type			        size_type;
    typedef typename base_container_t::difference_type	      difference_type;
    typedef typename base_container_t::reference			        reference;
    typedef typename base_container_t::const_reference	      const_reference;
    typedef typename base_container_t::value_type		          value_type;
    typedef Key						                                    key_type;
    typedef typename base_container_t::iterator			          iterator;
    typedef typename base_container_t::const_iterator	        const_iterator;
    typedef EqualKey						                              key_equal;
    typedef EqualKey						                              value_equal;
    
    typedef typename base_container_t::const_reverse_iterator const_reverse_iterator;
    typedef typename base_container_t::reverse_iterator	      reverse_iterator;
    
    typedef std::pair<iterator, iterator>                     Pairii_;
    typedef std::pair<const_iterator, const_iterator>         Paircc_;
    typedef std::pair<iterator, bool>                         Pairib_;

  protected:
    key_equal         key_equal_;
    base_container_t  storage_;

  public: // Construction & assignment
    
    explicit vector_set_base(const EqualKey& pred = EqualKey(),const Alloc& al = Alloc())
        : key_equal_(pred),storage_(al) {
    }

    vector_set_base(const this_t& x)
        : key_equal_(x.key_equal_),storage_(x.storage_) {
    }

    ~vector_set_base() {
    }
    
    this_t& operator=(const this_t& x) {
      (this->storage_).operator=(x.storage_);
      (this->key_equal_)= x.key_equal_;
      return *this;
    }

    this_t& operator=(const base_container_t& x){
      (this->storage_).operator=(x);
      return *this;
    }
    
    void reserve(size_type n) {
      (this->storage_).reserve(n);
    }
    
    Alloc get_allocator() const {
      return (this->storage_).get_allocator();
    }
    
  public: // Iterators
    
    iterator begin() {
      return (this->storage_).begin();
    }
    
    const_iterator begin() const {
      return (this->storage_).begin();
    }
    
    iterator end() {
      return (this->storage_).end();
    }
    
    const_iterator end() const {
      return (this->storage_).end();
    }
    
    reverse_iterator rbegin() {
      return (this->storage_).rbegin();
    }
    
    const_reverse_iterator rbegin() const {
      return (this->storage_).rbegin();
    }
    
    reverse_iterator rend() {
      return (this->storage_).rend();
    }
    
    const_reverse_iterator rend() const {
      return (this->storage_).rend();
    }

  public: // Serialization
    
    void store_to(output_serializer& s) const {
      s << (this->storage_);
    }
    
    void retrieve_from(input_serializer& s) {
      s >> (this->storage_);
    }

  public: // Size
    
    size_type size() const {
      return (this->storage_).size();
    }
    
    size_type max_size() const {
      return (this->storage_).max_size();
    }
    
    bool empty() const {
      return (this->storage_).empty();
    }

  public: // Element access
    
    const_reference at(size_type p) const {
      return (this->storage_).at(p);
    }
    
    reference at(size_type p) {
      return (this->storage_).at(p);
    }
    
    const_reference operator[](size_type p) const {
      return (this->storage_).operator[](p);
    }
		
    reference operator[](size_type p)	{
      return (this->storage_).operator[](p);
    }
    
    reference front() {
      return (this->storage_).front();
    }
    
    const_reference front() const {
      return (this->storage_).front();
    }
    
    reference back() {
      return (this->storage_).back();
    }
    
    const_reference back() const {
      return (this->storage_).back();
    }
    
  public: // Element insertion
    
    Pairib_ insert(const value_type& x) {
      if (bNoDuplicates) { 
        iterator p=find(x);
        if (p==end()) {
          return Pairib_((this->storage_).insert(p,x),true);
        }else{
          return Pairib_(p,false);
        }
      }else{
        return Pairib_((this->storage_).insert(end(),x),true);
      }
    }

    iterator insert(iterator /*it*/, const value_type& x) {
      return insert(x).first;
    }
    
#if (_MSC_VER >= 1300)
    template<class It>
    void insert(It first, It beyond) {
      size_type n= std::distance(first,beyond);
      reserve(size()+n);
      for( ;first!=beyond;++first){
        insert(*first);
      }
    }
#else
    void insert(const_iterator first, const_iterator beyond) {
      size_type n= std::distance(first,beyond);
      reserve(size()+n);
      for( ;first!=beyond;++first){
        insert(*first);
      }
    }
#endif

  public: // Element removal
    
    void pop_back() {
      (this->storage_).pop_back();
    }

    iterator erase(iterator p) {
      return (this->storage_).erase(p);
    }
    
    iterator erase(iterator first, iterator beyond) {
      return (this->storage_).erase(first,beyond);
    }

    size_type erase(const Key& key) {
      size_type n = 0;
      for (iterator it=begin();it!=end();) {
        if ((this->key_equal_)(key,*it)) {
          ++n;
          (this->storage_).erase(it);
          if (bNoDuplicates) break;
        } else {
          ++it;
        }
      }
      return n;
    }
    
    void clear() { 
      return (this->storage_).clear();
    }

  public: // comparison
    
    bool Eq_(const this_t& x) const {
      bool result = (size() == x.size());
      if (result) {
        const_iterator x_end = x.end();
        for (const_iterator it=begin(); result && it!=end(); ++it) {
          result = (x.find(*it) != x_end);
        }
      }
      return result;
    }

  public: // swapping
    
    void swap(this_t& x) {
      (this->storage_).swap(x.storage_);
      std::swap((this->key_equal_),x.key_equal_);
    }
    
    friend void swap(this_t& x, this_t& Y_) {
      x.swap(Y_);
    }

  public: // Comparison
    
    key_equal key_eq() const {
      return (this->key_equal_);
    }
    
    value_equal value_eq() const {
      return (this->key_eq());
    }

  public: // Search
    
    iterator find(const Key& k) {
      iterator it=begin();
      while (it!=end()) {
        if ((this->key_equal_)(k,*it)) break;
        ++it;
      }
      return it;
    }
    
    const_iterator find(const Key& k) const {
      const_iterator it=begin();
      while (it!=end()) {
        if ((this->key_equal_)(k,*it)) break;
        ++it;
      }
      return it;
    }
    
    size_type count(const Key& k) const {
      size_type n = 0;
      for (const_iterator it=begin();it!=end();++it) {
        if ((this->key_equal_)(k,*it)) {
          ++n;
          if (bNoDuplicates) break;
        }
      }
      return n;
    }
    
  public: // functions for use with direct std::vector access

    base_container_t& get_container() {
      return (this->storage_);
    }

  };


  template<class Key,bool bNoDuplicates,class EqualKey, class Alloc> inline
	bool operator==(const vector_set_base<Key, bNoDuplicates,EqualKey,Alloc>& x,
                  const vector_set_base<Key, bNoDuplicates,EqualKey,Alloc>& Y_)
	{return x.Eq_(Y_); }
  template<class Key,bool bNoDuplicates,class EqualKey, class Alloc> inline
	bool operator!=(const vector_set_base<Key, bNoDuplicates,EqualKey,Alloc>& x,
                  const vector_set_base<Key, bNoDuplicates,EqualKey,Alloc>& Y_)
	{return !(x == Y_); }

  /**
   * A std::set implemented as a non-sorted vector
   */
  template<class K, class EqualKey = std::equal_to<K>, class A = std::allocator<K> >
  class vector_set: public vector_set_base<K,true,EqualKey,A> {
  public:
    typedef vector_set_base<K,true,EqualKey,A> super_t;
    typedef vector_set<K,EqualKey,A>  this_t;
    
    typedef typename super_t::base_container_t		 base_container_t;
    typedef typename super_t::allocator_type	allocator_type;
    typedef typename super_t::size_type			size_type;
    typedef typename super_t::difference_type	difference_type;
    typedef typename super_t::reference			reference;
    typedef typename super_t::const_reference	const_reference;
    typedef typename super_t::value_type		value_type;
    typedef K						key_type;
    typedef typename super_t::iterator			iterator;
    typedef typename super_t::const_iterator	const_iterator;
    typedef EqualKey						key_equal;
    typedef EqualKey						value_equal;
    
    typedef typename super_t::const_reverse_iterator const_reverse_iterator;
    typedef typename super_t::reverse_iterator	reverse_iterator;
    
    typedef std::pair<iterator, iterator> Pairii_;
    typedef std::pair<const_iterator, const_iterator> Paircc_;
    typedef std::pair<iterator, bool> Pairib_;
  public:

    explicit vector_set(const EqualKey& pred = EqualKey(),
                               const A& al = A())
        : super_t(pred, al){}
#if (_MSC_VER >= 1300)
    template<class It>
    vector_set(It first, It beyond, 
                      const EqualKey& pred = EqualKey(),const A& al = A())
        :super_t(first, beyond, pred, al) {}
#else
    vector_set(const_iterator first, const_iterator beyond, 
                      const EqualKey& pred = EqualKey(),const A& al = A())
        :super_t(first, beyond, pred, al) {}
#endif
    vector_set(const this_t& x)
        : super_t(x)
    {}
    ~vector_set() {}
    
    this_t& operator=(const this_t& x) {
      (this->storage_).operator=(x.storage_);
      (this->key_equal_)= x.key_equal_;
      return *this;
    }
    this_t& operator=(const base_container_t& x){
      (this->storage_).operator=(x);
      return *this;
    }
    bool Eq_(const this_t& x) const{
      return (this->size() == x.size()
              && std::equal(this->begin(), this->end(), x.begin()));
    }
    void swap(this_t& x) {
      (this->storage_).swap(x.storage_);
      std::swap((this->key_equal_),x.key_equal_);
    }
    
    friend void swap(this_t& x, this_t& y) {
      x.swap(y);
    }
  };
  
  /**
   * A std::multiset implemented as a non-sorted vector
   */
	template<class K, class EqualKey = std::equal_to<K>, class A = std::allocator<K> >
	class vector_multiset: public vector_set_base<K,false,EqualKey,A> {
  public:
    typedef vector_set_base<K,true,EqualKey,A> super_t;
    typedef vector_multiset<K,EqualKey,A>  this_t;
    
    typedef typename super_t::base_container_t		 base_container_t;
    typedef typename super_t::allocator_type	allocator_type;
    typedef typename super_t::size_type			size_type;
    typedef typename super_t::difference_type	difference_type;
    typedef typename super_t::reference			reference;
    typedef typename super_t::const_reference	const_reference;
    typedef typename super_t::value_type		value_type;
    typedef K						key_type;
    typedef typename super_t::iterator			iterator;
    typedef typename super_t::const_iterator	const_iterator;
    typedef EqualKey						key_equal;
    typedef EqualKey						value_equal;
    
    typedef typename super_t::const_reverse_iterator const_reverse_iterator;
    typedef typename super_t::reverse_iterator	reverse_iterator;
    
    typedef std::pair<iterator, iterator> Pairii_;
    typedef std::pair<const_iterator, const_iterator> Paircc_;
    typedef std::pair<iterator, bool> Pairib_;
  public:

    explicit vector_multiset(const EqualKey& pred = EqualKey(),
                               const A& al = A())
        : super_t(pred, al){}
#if (_MSC_VER >= 1300)
    template<class It>
    vector_multiset(It first, It beyond, 
                      const EqualKey& pred = EqualKey(),const A& al = A())
        :super_t(first, beyond, pred, al) {}
#else
    vector_multiset(const_iterator first, const_iterator beyond, 
                      const EqualKey& pred = EqualKey(),const A& al = A())
        :super_t(first, beyond, pred, al) {}
#endif
    vector_multiset(const this_t& x)
        : super_t(x)
    {}
    ~vector_multiset() {}
    this_t& operator=(const this_t& x) {
      (this->storage_).operator=(x.storage_);
      (this->key_equal_)= x.key_equal_;
      return *this;
    }
    this_t& operator=(const base_container_t& x) {
      (this->storage_).operator=(x);
      return *this;
    }
    bool Eq_(const this_t& x) const{
       return (this->size() == x.size()
               && std::equal(this->begin(), this->end(), x.begin()));
    }
    void swap(this_t& x)
          {(this->storage_).swap(x.storage_);std::swap((this->key_equal_),x.key_equal_);}
    
    friend void swap(this_t& x, this_t& Y_)
          {x.swap(Y_); }
  };
}

#elif SL_VERSION_VECTOR_SET_ != 0x00010010
#error You have included two vector_set with different version numbers
#endif
