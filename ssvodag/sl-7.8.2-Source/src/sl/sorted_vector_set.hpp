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
/* STL-conforming "sorted vector" container
 *
 * Based on the following code:
 *
 * (C) 2002 Martin Holzherr (holzherr@infobrain.com). All rights reserved.
 *
 * Permission is granted to use, distribute and modify this code provided that:
 *   · this copyright notice appears,
 *   · 
 * The author welcomes any suggestions on the code or reportings of actual
 * use of the code. Please send your comments to holzherr@infobrain.com.
 *
 * The author makes NO WARRANTY or representation, either express or implied,
 * with respect to this code, its quality, accuracy, merchantability, or
 * fitness for a particular purpose.  This software is provided "AS IS", and
 * you, its user, assume the entire risk as to its quality and accuracy.
 *
 * Created:			November 19th, 2002
 * Last modified:	November 27th, 2002 								
 */

#ifndef SL_SORTED_VECTOR_SET_HPP
#define SL_SORTED_VECTOR_SET_HPP
#define SL_VERSION_SORTED_VECTOR_SET_ 0x00010010


#include <algorithm>
#include <vector>
#include <utility>
#include <functional>
#include <sl/serializer.hpp>

namespace sl {

  /**
   * A sorted vector conforming to std::set when bNoDuplicates=true and std::multiset when
   * bNoDuplicates=false
   */
  template<class K, bool bNoDuplicates= false,class Pr = std::less<K>, class A = std::allocator<K> >
  class sorted_vector_set_base {
  public:
    typedef sorted_vector_set_base<K,bNoDuplicates,Pr,A> this_t;
    typedef std::vector<K,A>		base_container_t;
    typedef typename base_container_t::allocator_type	allocator_type;
    typedef typename base_container_t::size_type			size_type;
    typedef typename base_container_t::difference_type	difference_type;
    typedef typename base_container_t::reference			reference;
    typedef typename base_container_t::const_reference	const_reference;
    typedef typename base_container_t::value_type		value_type;
    typedef K						key_type;
    typedef typename base_container_t::iterator			iterator;
    typedef typename base_container_t::const_iterator	const_iterator;
    typedef Pr						key_compare;
    typedef Pr						value_compare;
    
    typedef typename base_container_t::const_reverse_iterator const_reverse_iterator;
    typedef typename base_container_t::reverse_iterator	reverse_iterator;
    
    typedef std::pair<iterator, iterator> Pairii_;
    typedef std::pair<const_iterator, const_iterator> Paircc_;
    typedef std::pair<iterator, bool> Pairib_;
  protected:
    key_compare         key_compare_;
    base_container_t                storage_;

  public: // Construction & assignment
    
    explicit sorted_vector_set_base(const Pr& pred = Pr(),const A& al = A())
        : key_compare_(pred),storage_(al)
    {
    }
    
#if (_MSC_VER >= 1300)
    template<class It>
    sorted_vector_set_base(It first, It beyond, 
                  const Pr& pred = Pr(),const A& al = A())
        : key_compare_(pred),storage_(first,beyond,al)
    {
      stable_sort();
    }
#else
    sorted_vector_set_base(const_iterator first, const_iterator beyond, 
                  const Pr& pred = Pr(),const A& al = A())
        : key_compare_(pred),storage_(first,beyond,al)
    {
      stable_sort();
    }
#endif
    sorted_vector_set_base(const this_t& x)
        : key_compare_(x.key_compare_),storage_(x.storage_)
    {
    }

    ~sorted_vector_set_base() {
    }
    
    this_t& operator=(const this_t& x) {
      (this->storage_).operator=(x.storage_);
      (this->key_compare_)= x.key_compare_;
      return *this;
    }

    this_t& operator=(const base_container_t& x){
      (this->storage_).operator=(x);
      sort();return *this;
    }
    
    void reserve(size_type n) {
      (this->storage_).reserve(n);
    }
    
    A get_allocator() const {
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
    
  public: // Removal
    
    void pop_back()								{
      (this->storage_).pop_back();
    }

  public: // Assignment
    
    void assign(const_iterator first, const_iterator beyond) {					
      (this->storage_).assign(first,beyond);
    }
    
    void assign(size_type n, const K& x = K()) {
      (this->storage_).assign(n,x);
    }
    
  public: // Element insertion
    
    Pairib_ insert(const value_type& x) {
      if(bNoDuplicates){
        iterator p= lower_bound(x);
        if(p==end()||(this->key_compare_)(x,*p)){
          return Pairib_(InsertImpl_(p,x),true);
        }else{
          return Pairib_(p,false);
        }
      }else{
        iterator p= upper_bound(x);
        return Pairib_(InsertImpl_(p,x),true);
      }
    }

    iterator insert(iterator it, const value_type& x) {
      // it is the int
      if(it!=end() ){
        if(bNoDuplicates){
          if((this->key_compare_)(*it,x)){
            if((it+1)==end()||KeyCompare_Gt_(*(it+1),x)){//use hint
              return InsertImpl_(it+1,x);
            }else if(KeyCompare_Geq_(*(it+1),x)){
              return end();
            }
          }
        }else {
          if(	KeyCompare_Leq_(*it,x)
              &&((it+1)==end()||KeyCompare_Geq_(*(it+1),x))){
            return InsertImpl_(it+1,x);
          }
        }
      }
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
    
    iterator erase(iterator p) {
      return (this->storage_).erase(p);
    }
    
    iterator erase(iterator first, iterator beyond) {
      return (this->storage_).erase(first,beyond);
    }

    size_type erase(const K& key) {
      Pairii_ begEnd= equal_range(key);
      size_type n= std::distance(begEnd.first,begEnd.second);
      erase(begEnd.first,begEnd.second);
      return n;
    }
    
    void clear(){
      return (this->storage_).clear();
    }

  public: // comparison
    
    bool Eq_(const this_t& x) const {
      return (size() == x.size()
              && std::equal(begin(), end(), x.begin()));
    }
    
    bool Lt_(const this_t& x) const {
      return (std::lexicographical_compare(begin(), end(),
                                           x.begin(), x.end()));
    }

  public: // swapping
    
    void swap(this_t& x) {
      (this->storage_).swap(x.storage_);
      std::swap((this->key_compare_),x.key_compare_);
    }
    
    friend void swap(this_t& x, this_t& Y_) {
      x.swap(Y_);
    }

  public: // Comparison
    
    key_compare key_comp() const {
      return (this->key_compare_);
    }
    
    value_compare value_comp() const {
      return (this->key_comp());
    }

  public: // Search
    
    iterator find(const K& k) {
      iterator p = lower_bound(k);
      return (p==end()||(this->key_compare_)(k, *p))? end():p;
		}
    
    const_iterator find(const K& k) const {
      const_iterator p = lower_bound(k);
      return (p==end()||(this->key_compare_)(k,*p))?end():p;
    }
    
    size_type count(const K& k) const {
      Paircc_ Ans_ = equal_range(k);
      size_type n = std::distance(Ans_.first, Ans_.second);
      return (n);
    }
    
    iterator lower_bound(const K& k) {
      return std::lower_bound(begin(), end(), k, (this->key_compare_));
    }
    
    const_iterator lower_bound(const K& k) const {
      return std::lower_bound(begin(), end(), k, (this->key_compare_));
    }
    
    iterator upper_bound(const K& k) {
      return std::upper_bound(begin(), end(), k, (this->key_compare_));
    }
    
    const_iterator upper_bound(const K& k) const {
      return std::upper_bound(begin(), end(), k, (this->key_compare_));
    }
    
    Pairii_ equal_range(const K& k) {
      return std::equal_range(begin(), end(), k, (this->key_compare_));
    }
    
    Paircc_ equal_range(const K& k) const {
      return std::equal_range(begin(), end(), k, (this->key_compare_));
    }
    
  public: // functions for use with direct std::vector access

    base_container_t& get_container() {
      return (this->storage_);
    }

    //restore sorted order after low level access 
    void sort() {
      std::sort((this->storage_).begin(),(this->storage_).end(),(this->key_compare_));
      if( bNoDuplicates ){
        (this->storage_).erase(Unique_(),(this->storage_).end());
      }
    }
    
    //restore sorted order after low level access 
    void stable_sort() {
      std::stable_sort((this->storage_).begin(),(this->storage_).end(),(this->key_compare_));
      if( bNoDuplicates ){
        erase(Unique_(),end());
      }
    }
    
  protected:
    iterator Unique_() {
      iterator front_= (this->storage_).begin(),out_= (this->storage_).end(),end_=(this->storage_).end();
      bool bCopy_= false;
      for(iterator prev_; (prev_=front_)!=end_ && ++front_!=end_; ){
        if( (this->key_compare_)(*prev_,*front_)){
          if(bCopy_){
            *out_= *front_;
            out_++;
          }
        }else{
          if(!bCopy_){out_=front_;bCopy_=true;}
        }
      }
      return out_;
    }
    iterator InsertImpl_(iterator p,const value_type& x) {
      return (this->storage_).insert(p,x);
    }
    bool KeyCompare_Leq_(const K& ty0,const K& ty1) {
      return !(this->key_compare_)(ty1,ty0);
    }
    
    bool KeyCompare_Geq_(const K& ty0,const K& ty1) {
      return !(this->key_compare_)(ty0,ty1);
    }
    
    bool KeyCompare_Gt_(const K& ty0,const K& ty1) {
      return (this->key_compare_)(ty1,ty0);
    }
  };


  template<class K,bool bNoDuplicates,class Pr, class A> inline
	bool operator==(const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& x,
                  const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& Y_)
	{return x.Eq_(Y_); }
  template<class K,bool bNoDuplicates,class Pr, class A> inline
	bool operator!=(const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& x,
                  const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& Y_)
	{return !(x == Y_); }
  template<class K,bool bNoDuplicates,class Pr, class A> inline
	bool operator<(const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& x,
                 const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& Y_)
	{return x.Lt_(Y_);}
  template<class K,bool bNoDuplicates,class Pr,class A> inline
	bool operator>(const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& x,
                 const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& Y_)
	{return Y_ < x; }
  template<class K,bool bNoDuplicates,class Pr, class A> inline
	bool operator<=(const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& x,
                  const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& Y_)
	{return !(Y_ < x); }
  template<class K, bool bNoDuplicates,class Pr,class A> inline
	bool operator>=(const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& x,
                  const sorted_vector_set_base<K, bNoDuplicates,Pr,A>& Y_)
	{return (!(x < Y_)); }


  /**
   * A std::set implemented as a sorted vector
   */
  template<class K, class Pr = std::less<K>, class A = std::allocator<K> >
  class sorted_vector_set: public sorted_vector_set_base<K,true,Pr,A> {
  public:
    typedef sorted_vector_set_base<K,true,Pr,A> super_t;
    typedef sorted_vector_set<K,Pr,A>  this_t;
    
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
    typedef Pr						key_compare;
    typedef Pr						value_compare;
    
    typedef typename super_t::const_reverse_iterator const_reverse_iterator;
    typedef typename super_t::reverse_iterator	reverse_iterator;
    
    typedef std::pair<iterator, iterator> Pairii_;
    typedef std::pair<const_iterator, const_iterator> Paircc_;
    typedef std::pair<iterator, bool> Pairib_;
  public:

    explicit sorted_vector_set(const Pr& pred = Pr(),
                               const A& al = A())
        : super_t(pred, al){}
#if (_MSC_VER >= 1300)
    template<class It>
    sorted_vector_set(It first, It beyond, 
                      const Pr& pred = Pr(),const A& al = A())
        :super_t(first, beyond, pred, al) {}
#else
    sorted_vector_set(const_iterator first, const_iterator beyond, 
                      const Pr& pred = Pr(),const A& al = A())
        :super_t(first, beyond, pred, al){}
#endif
    sorted_vector_set(const this_t& x)
        : super_t(x)
    {}
    ~sorted_vector_set() {}
    
    this_t& operator=(const this_t& x) {
      (this->storage_).operator=(x.storage_);
      (this->key_compare_)= x.key_compare_;
      return *this;
    }
    this_t& operator=(const base_container_t& x){
      (this->storage_).operator=(x);
      this->sort();
      return *this;
    }
    bool Eq_(const this_t& x) const {
      return (this->size() == x.size()
              && std::equal(this->begin(), this->end(), x.begin()));
    }
    bool Lt_(const this_t& x) const {
      return (std::lexicographical_compare(this->begin(), this->end(),
                                           x.begin(), x.end()));
    }
    void swap(this_t& x) {
      (this->storage_).swap(x.storage_);
      std::swap((this->key_compare_),x.key_compare_);
    }
    
    friend void swap(this_t& x, this_t& Y_){
      x.swap(Y_);
    }
  };
  
  /**
   * A std::multiset implemented as a sorted vector
   */
  template<class K, class Pr = std::less<K>, class A = std::allocator<K> >
  class sorted_vector_multiset: public sorted_vector_set_base<K,false,Pr,A> {
  public:
    typedef sorted_vector_set_base<K,false,Pr,A> super_t;
    typedef sorted_vector_multiset<K,Pr,A>  this_t;
    
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
    typedef Pr						key_compare;
    typedef Pr						value_compare;
    
    typedef typename super_t::const_reverse_iterator const_reverse_iterator;
    typedef typename super_t::reverse_iterator	reverse_iterator;
    
    typedef std::pair<iterator, iterator> Pairii_;
    typedef std::pair<const_iterator, const_iterator> Paircc_;
    typedef std::pair<iterator, bool> Pairib_;
  public:

    explicit sorted_vector_multiset(const Pr& pred = Pr(),
                               const A& al = A())
        : super_t(pred, al){}
#if (_MSC_VER >= 1300)
    template<class It>
    sorted_vector_multiset(It first, It beyond, 
                      const Pr& pred = Pr(),const A& al = A())
        :super_t(first, beyond, pred, al) {}
#else
    sorted_vector_multiset(const_iterator first, const_iterator beyond, 
                      const Pr& pred = Pr(),const A& al = A())
        :super_t(first, beyond, pred, al) {}
#endif
    sorted_vector_multiset(const this_t& x)
        : super_t(x)
    {}
    ~sorted_vector_multiset() {}
    this_t& operator=(const this_t& x) {
      (this->storage_).operator=(x.storage_);
      (this->key_compare_)= x.key_compare_;
      return *this;
    }
    this_t& operator=(const base_container_t& x){
      (this->storage_).operator=(x);
      this->sort();
      return *this;
    }
    bool Eq_(const this_t& x) const {
      return (this->size() == x.size()
              && std::equal(this->begin(), this->end(), x.begin()));
    }
    bool Lt_(const this_t& x) const {
      return (std::lexicographical_compare(this->begin(),
                                           this->end(),
                                           x.begin(),
                                           x.end()));
    }
    void swap(this_t& x) {
      (this->storage_).swap(x.storage_);
      std::swap((this->key_compare_),x.key_compare_);
    }
    
    friend void swap(this_t& x, this_t& y) {
      x.swap(y);
    }
  };
}

#elif SL_VERSION_SORTED_VECTOR_SET_ != 0x00010010
#error You have included two sorted_vector_set with different version numbers
#endif
