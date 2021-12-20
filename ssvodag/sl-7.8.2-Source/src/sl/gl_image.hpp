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
#ifndef SL_GL_IMAGE_HPP
#define SL_GL_IMAGE_HPP

#include <sl/math.hpp>
#include <sl/numeric_traits.hpp>
#include <sl/utility.hpp>
#include <sl/fixed_size_array.hpp>

// Standard C++ includes
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <cmath>

#include <algorithm>
#include <cassert>

namespace sl {

  /**
   *  A simple in-core image class of up to
   *  3 dimensions using OpenGL standard 
   *  texture storage conventions. Tuned
   *  for small images.
   */
  template <class T=unsigned char>
  class gl_image {
  public:
    typedef gl_image<T> this_t;
    typedef T  value_type;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
  protected:

    sl::fixed_size_array<4,std::size_t> extent_; // c w h d
    value_type* data_;

  public:

    /// Construct by default (null image)
    inline gl_image() {
      extent_[0] = 0;
      extent_[1] = 0;
      extent_[2] = 0;
      extent_[3] = 0;
      data_ = 0;
    }
      
    /// Construct image of given size
    inline gl_image(std::size_t dimc,
		    std::size_t dimx, 
		    std::size_t dimy = 1,
		    std::size_t dimz = 1) {
      extent_[0] = dimc;
      extent_[1] = dimx;
      extent_[2] = dimy;
      extent_[3] = dimz;
      data_ = new value_type[dimc*dimx*dimy*dimz];
    }

    /// Destruct 
    inline ~gl_image() {
      if (data_) delete [] data_;
      data_ = 0;
    }
      
    /// Delete contents and resize to null
    inline this_t& wipe_out() {
      if (data_) delete [] data_;
      data_ = 0;
      extent_[0] = 0;
      extent_[1] = 0;
      extent_[2] = 0;
      extent_[3] = 0;
      return *this;
    }

    /// Construct as copy of other
    gl_image(const this_t& other) {
      extent_ = other.extent_;
      std::size_t other_size = other.size();
      data_ = other_size ? (new value_type[other_size]) : 0;
      for (std::size_t i=0; i<other_size; ++i) {
	data_[i] = other.data_[i];
      }
    }

    /// Construct as copy of other
    template<typename U> 
    gl_image(const gl_image<U>& other) {
      extent_ = other.extent_;
      std::size_t other_size = other.size();
      data_ = other_size ? (new value_type[other_size]) : 0;
      for (std::size_t i=0; i<other_size; ++i) {
	data_[i] = static_cast<value_type>(other.data_[i]);
      }
    }

    /// Assign to copy of other
    template<typename U> 
    this_t& assign(const gl_image<U>& other) {
      std::size_t other_size = other.size();
      if (size() != other_size) {
	// Needs a resize of data array
	wipe_out();
	data_ = other_size ? (new value_type[other_size]) : 0;
      }

      // Copy
      extent_ = other.extent_;
      for (std::size_t i=0; i<other_size; ++i) {
	data_[i] = static_cast<value_type>(other.data_[i]);
      }

      return *this;
    }

    /// Assign to copy of data stored in C array
    template<typename U> 
    this_t& assign(const U* other_data,
		   std::size_t dimc,
		   std::size_t dimx,
		   std::size_t dimy = 1,
		   std::size_t dimz = 1) {
      std::size_t other_size = dimc*dimx*dimy*dimz;
      assert((other_size == 0) || other_data);
 
      if (size() != other_size) {
	// Needs a resize of data array
	wipe_out();
	data_ = other_size ? (new value_type[other_size]) : 0;
      }
      // Copy
      extent_[0] = dimc;
      extent_[1] = dimx;
      extent_[2] = dimy;
      extent_[3] = dimz;
      for (std::size_t i=0; i<other_size; ++i) {
	data_[i] = static_cast<value_type>(other_data[i]);
      }

      return *this;
    }

    /// Assign to copy of data stored in C array
    template<typename U> 
    gl_image(const U* other_data,
	     std::size_t dimc,
	     std::size_t dimx,
	     std::size_t dimy = 1,
	     std::size_t dimz = 1) {
      extent_[0] = 0;
      extent_[1] = 0;
      extent_[2] = 0;
      extent_[3] = 0;
      data_ = 0;
      assign(other_data, dimc, dimx, dimy, dimz);
    }

    /// Assign to copy of other
    this_t& operator = (const this_t& other) {
      return assign(other);
    }

    /// Assign to copy of other
    template<typename U> 
    this_t& operator = (const gl_image<U>& other) {
      return assign(other);
    }

    /// Swap fields of an image (use it carefully!)
    void swap(this_t& other) {
      std::swap(extent_[0],other.extent_[0]);
      std::swap(extent_[1],other.extent_[1]);
      std::swap(extent_[2],other.extent_[2]);
      std::swap(extent_[3],other.extent_[3]);
      std::swap(data_, other.data_);
    }
    
  public: // Access

    inline std::size_t extent(std::size_t i) const {
      assert(i<4);
      return extent_[i];
    }

    inline std::size_t channels() const { return extent_[0]; }

    inline std::size_t width() const { return extent_[1]; }

    inline std::size_t height() const { return extent_[2]; }

    inline std::size_t depth() const { return extent_[3]; }

    inline std::size_t size() const { return extent_[0]*extent_[1]*extent_[2]*extent_[3]; }
      
    inline bool is_empty() const { return size() == 0; }
      
    inline std::size_t offset() const {
      return 0;
    }

    inline std::size_t offset(int x) const {
      return extent_[0]*x;
    }

    inline std::size_t offset(int x, int y) const {
      return extent_[0]*(extent_[1]*y + x );
    }

    inline std::size_t offset(int x, int y, int z) const {
      return extent_[0]*(extent_[1]*(extent_[2]*z + y) + x );
    }

    const value_type* to_pointer() const {
      return data_;
    }

    const value_type* to_pointer(int x) const {
      return &(data_[offset(x)]);
    }

    const value_type* to_pointer(int x, int y) const {
      return &(data_[offset(x, y)]);
    }
            
    const value_type* to_pointer(int x, int y, int z) const {
      return &(data_[offset(x, y, z)]);
    }

    value_type* to_pointer() {
      return data_;
    }

    value_type* to_pointer(int x) {
      return &(data_[offset(x)]);
    }

    value_type* to_pointer(int x, int y) {
      return &(data_[offset(x, y)]);
    }
            
    value_type* to_pointer(int x, int y, int z) {
      return &(data_[offset(x, y, z)]);
    }
            
    const value_type& at(int c, int x) const {
      assert(!is_empty());
      assert(c>= 0 && c<int(channels()));
      assert(x>= 0 && x<int(width()));
      return *(to_pointer(x)+c);
    }

    const value_type& at(int c, int x, int y) const {
      assert(!is_empty());
      assert(c>= 0 && c<int(channels()));
      assert(x>= 0 && x<int(width()));
      assert(y>= 0 && y<int(height()));
      return *(to_pointer(x, y)+c);
    }

    const value_type& at(int c, int x, int y, int z) const {
      assert(!is_empty());
      assert(c>= 0 && c<int(channels()));
      assert(x>= 0 && x<int(width()));
      assert(y>= 0 && y<int(height()));
      assert(z>= 0 && z<int(depth()));
      return *(to_pointer(x, y, z)+c);
    }

    value_type& at(int c, int x) {
      assert(!is_empty());
      assert(c>= 0 && c<int(channels()));
      assert(x>= 0 && x<int(width()));
      return *(to_pointer(x)+c);
    }

    value_type& at(int c, int x, int y) {
      assert(!is_empty());
      assert(c>= 0 && c<int(channels()));
      assert(x>= 0 && x<int(width()));
      assert(y>= 0 && y<int(height()));
      return *(to_pointer(x, y)+c);
    }

    value_type& at(int c, int x, int y, int z) {
      assert(!is_empty());
      assert(c>= 0 && c<int(channels()));
      assert(x>= 0 && x<int(width()));
      assert(y>= 0 && y<int(height()));
      assert(z>= 0 && z<int(depth()));
      return *(to_pointer(x, y, z)+c);
    }

    const value_type& operator()(int c, int x) const {
      return at(c,x);
    }

    const value_type& operator()(int c, int x, int y) const {
      return at(c,x,y);
    }

    const value_type& operator()(int c, int x, int y, int z) const {
      return at(c,x,y,z);
    }

    value_type& operator()(int c, int x) {
      return at(c,x);
    }

    value_type& operator()(int c, int x, int y) {
      return at(c,x,y);
    }

    value_type& operator()(int c, int x, int y, int z) {
      return at(c,x,y,z);
    }

    const_iterator begin() const {
      return data_;
    }

    const_iterator end() const {
      return data_+size();
    }

    iterator begin() {
      return data_;
    }

    iterator end() {
      return data_+size();
    }

    const value_type& front() const {
      assert(!is_empty());
      return data_[0];
    }

    const value_type& back() const {
      assert(!is_empty());
      return data_[size()-1];
    }

    value_type& front() {
      assert(!is_empty());
      return data_[0];
    }

    value_type& back() {
      assert(!is_empty());
      return data_[size()-1];
    }

  public: // Init

    /// Fill sequentially all image with value c
    this_t& fill(const value_type& c) {
      std::size_t sz = size();
      for (std::size_t i=0; i<sz; ++i) {
	data_[i] = c;
      }
      return *this;
    }

    /// Fill sequentially all image with values c0 c1
    this_t& fill(const value_type& c0, 
		 const value_type& c1) {
      if (!is_empty()) {
	value_type *ptr, *ptr_end = end()-1;
	for (ptr=begin(); ptr<ptr_end; ) { 
	  *(ptr++)=c0; 
	  *(ptr++)=c1; 
	}
	if (ptr!=ptr_end+1) *(ptr++)=c0;
      }
      return *this;
    }

    /// Fill sequentially all image with values c0 c1 c2
    this_t& fill(const value_type& c0, 
		 const value_type& c1,
		 const value_type& c2) {
      if (!is_empty()) {
	value_type *ptr, *ptr_end = end()-2;
	for (ptr=begin(); ptr<ptr_end; ) { 
	  *(ptr++)=c0; 
	  *(ptr++)=c1; 
	  *(ptr++)=c2;
	}
	ptr_end+=2;
	switch (ptr_end-ptr) {
	case 2: *(--ptr_end)=c1;
	case 1: *(--ptr_end)=c0;
	}
      }
      return *this;
    }

    /// Fill sequentially all image with values c0 c1 c2 c3
    this_t& fill(const value_type& c0, 
		 const value_type& c1,
		 const value_type& c2,
		 const value_type& c3) {
      if (!is_empty()) {
	value_type *ptr, *ptr_end = end()-3;
	for (ptr=begin(); ptr<ptr_end; ) { 
	  *(ptr++)=c0; 
	  *(ptr++)=c1; 
	  *(ptr++)=c2;
	  *(ptr++)=c3;
	}
	ptr_end+=3;
	switch (ptr_end-ptr) {
	case 3: *(--ptr_end)=c2;
	case 2: *(--ptr_end)=c1;
	case 1: *(--ptr_end)=c0;
	}
      }
      return *this;
    }

  public: // Alpha

    /// Fill channel c with value v
    this_t& fill_channel(int c,
			 const value_type& v) {
      assert(c>= 0 && c<int(channels()));

      std::size_t nc=channels();

      value_type *ptr, *ptr_end = end();
      for (ptr=begin(); ptr<ptr_end; ptr+=nc) {
	ptr[c] = v;
      }

      return *this;
    }

  public: // Alpha

    /// rescale channel c between min_out and max_out
    void channel_rescale(std::size_t c, 
			 value_type min_out, 
			 value_type max_out,
			 value_type in_zero = scalar_math<value_type>::zero(),
			 value_type in_one = scalar_math<value_type>::finite_upper_bound()) {
      assert(c<channels());
      std::size_t nc=channels()-1;
      const float zero = float(in_zero);
      const float one  = float(in_one);
      value_type *ptr, *ptr_end = end();
      for (ptr=begin(); ptr<ptr_end; ptr+=nc+1) {
	float v =  ((one-float(ptr[c]))*min_out+(float(ptr[c])-zero)*max_out)/(one-zero);
	ptr[c] = static_cast<value_type>(median(v,zero,one)/*+0.4*/);
      }
    }

    /// Set last channel to zero if average of others is less than v, to max instead
    this_t& set_alpha_from_black(value_type v=10,
				 value_type zero = scalar_math<value_type>::zero(),
				 value_type one  = scalar_math<value_type>::finite_upper_bound()) {
      assert(channels()>1);
      std::size_t nc=channels()-1;

      SL_SUMTYPENAME(value_type) thr = nc * v;

      value_type *ptr, *ptr_end = end();
      for (ptr=begin(); ptr<ptr_end; ptr+=nc+1) {
	SL_SUMTYPENAME(value_type) sum = zero;
	for (std::size_t i=0; i<nc; ++i) {
	  sum += ptr[i];
	}
	ptr[nc] = (sum<=thr) ? zero : one;
      }
      return *this;
    }

    /// Set last channel to zero if average of others is more than v, to max instead
    this_t& set_alpha_from_white(value_type v=250,
				 value_type zero = scalar_math<value_type>::zero(),
				 value_type one  = scalar_math<value_type>::finite_upper_bound()) {
      assert(channels()>1);
      std::size_t nc=channels()-1;

      SL_SUMTYPENAME(value_type) thr = nc * v;

      value_type *ptr, *ptr_end = end();
      for (ptr=begin(); ptr<ptr_end; ptr+=nc+1) {
	SL_SUMTYPENAME(value_type) sum = zero;
	for (std::size_t i=0; i<nc; ++i) {
	  sum += ptr[i];
	}
	ptr[nc] = (sum>=thr) ? zero : one;
      }
      return *this;
    }

    /// Set last channel to average of others 
    this_t& set_alpha_from_value() {
      assert(channels()>1);
      std::size_t nc=channels()-1;

      value_type zero = scalar_math<value_type>::zero();
	
      value_type *ptr, *ptr_end = end();
      for (ptr=begin(); ptr<ptr_end; ptr+=nc+1) {
	SL_SUMTYPENAME(value_type) sum = zero;
	for (std::size_t i=0; i<nc; ++i) {
	  sum += ptr[i];
	}
	ptr[nc] = value_type(sum/nc);
      }
      return *this;
    }

  public: // Crop

    /// Reshape image to have given extents
    this_t& reshape(std::size_t dimc,
		    std::size_t dimx, 
		    std::size_t dimy = 1,
		    std::size_t dimz = 1) {
      std::size_t new_sz = dimc*dimx*dimy*dimz;
      if (new_sz != size()) {
	wipe_out();
	data_ = new_sz ? (new value_type[new_sz]) : 0;
      }
      extent_[0] = dimc;
      extent_[1] = dimx;
      extent_[2] = dimy;
      extent_[3] = dimz;

      return *this;
    }

    this_t& assign_crop(const this_t& src,
			int c0, int x0, int y0, int z0,
			int c1, int x1, int y1, int z1) {
      assert(c0>= 0 && c0<src.channels());
      assert(x0>= 0 && x0<src.width());
      assert(y0>= 0 && y0<src.height());
      assert(z0>= 0 && z0<src.depth());
      assert(c1>= c0 && c1<src.channels());
      assert(x1>= x0 && x1<src.width());
      assert(y1>= y0 && y1<src.height());
      assert(z1>= z0 && z1<src.depth());

      const std::size_t dx=x1-x0+1;
      const std::size_t dy=y1-y0+1;
      const std::size_t dz=z1-z0+1;
      const std::size_t dc=c1-c0+1;

      reshape(dc,dx,dy,dz);

      const value_type *psrc = src.to_pointer(c0,x0,y0,z0);
      value_type *pdest = to_pointer();
	
      // FIXME optimize most common cases
      for (unsigned int z=0; z<dz; ++z) {
	for (unsigned int y=0; y<dy; ++y) {
	  for (unsigned int x=0; x<dx; ++x) {
	    std::memcpy(pdest,psrc,dc*sizeof(value_type)); // works only for POD
	    pdest+=dc;
	    psrc+=extent_[0]; // channels
	  }
	  psrc+=extent_[0]*(extent_[1]-dx); // skip to crop begin
	}
	psrc+=extent_[0]*extent_[1]*(extent_[2]-dy); // skip to crop begin
      }

      return *this;
    }

    this_t& assign_crop(const this_t& src,
			int x0, int y0, int z0,
			int x1, int y1, int z1) {
      int c0= 0; int c1 = int(src.channels())-1;
      assert(c0>= 0 && c0<src.channels());
      assert(x0>= 0 && x0<src.width());
      assert(y0>= 0 && y0<src.height());
      assert(z0>= 0 && z0<src.depth());
      assert(c1>= c0 && c1<src.channels());
      assert(x1>= x0 && x1<src.width());
      assert(y1>= y0 && y1<src.height());
      assert(z1>= z0 && z1<src.depth());
      assign_crop(src, 
		  c0, x0, y0, z0,
		  c1, x1, y1, z1);

      return *this;
    }
 
    this_t& assign_crop(const this_t& src,
			int x0, int y0,
			int x1, int y1) {
      int c0= 0; int c1 = int(src.channels())-1;
      int z0= 0; int z1 = int(src.depth())-1;
      assert(c0>= 0 && c0<src.channels());
      assert(x0>= 0 && x0<src.width());
      assert(y0>= 0 && y0<src.height());
      assert(z0>= 0 && z0<src.depth());
      assert(c1>= c0 && c1<src.channels());
      assert(x1>= x0 && x1<src.width());
      assert(y1>= y0 && y1<src.height());
      assert(z1>= z0 && z1<src.depth());
      assign_crop(src, 
		  c0, x0, y0, z0,
		  c1, x1, y1, z1);

      return *this;
    }
      
    this_t& assign_crop(const this_t& src,
			int x0, 
			int x1) {
      int c0= 0; int c1 = int(src.channels())-1;
      int z0= 0; int z1 = int(src.depth())-1;
      int y0= 0; int y1 = int(src.height())-1;
      assert(c0>= 0 && c0<src.channels());
      assert(x0>= 0 && x0<src.width());
      assert(y0>= 0 && y0<src.height());
      assert(z0>= 0 && z0<src.depth());
      assert(c1>= c0 && c1<src.channels());
      assert(x1>= x0 && x1<src.width());
      assert(y1>= y0 && y1<src.height());
      assert(z1>= z0 && z1<src.depth());
      assign_crop(src, 
		  c0, x0, y0, z0,
		  c1, x1, y1, z1);

      return *this;
    }

    this_t& crop(int c0, int x0, int y0, int z0,
		 int c1, int x1, int y1, int z1) {
      this_t tmp;
      tmp.assign_crop(c0, x0, y0, z0,
		      c1, x1, y1, z1);
      tmp.swap(*this);

      return *this;
    }

    this_t& crop(int x0, int y0, int z0,
		 int x1, int y1, int z1) {
      this_t tmp;
      tmp.assign_crop(x0, y0, z0,
		      x1, y1, z1);
      tmp.swap(*this);

      return *this;
    }
        
    this_t& crop(int x0, int y0,
		 int x1, int y1) {
      this_t tmp;
      tmp.assign_crop(x0, y0,
		      x1, y1);
      tmp.swap(*this);

      return *this;
    }

    this_t& crop(int x0,
		 int x1) {
      this_t tmp;
      tmp.assign_crop(x0,
		      x1);
      tmp.swap(*this);

      return *this;
    }

  };
} // namespace sl

#endif
