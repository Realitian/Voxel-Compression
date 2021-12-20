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
#ifndef SL_HISTOGRAM_HPP
#define SL_HISTOGRAM_HPP

#include <sl/config.hpp>
#include <sl/std_serializer.hpp>
#include <sl/cstdint.hpp>
#include <sl/math.hpp>
#include <vector>
#include <algorithm>
#include <cassert>

namespace sl {

  template <class T>
  class histogram {
  public:
    typedef T             value_t;
    typedef histogram<T>  this_t;
  protected:
    std::vector <value_t> bin_count_x_; //! Counters for bins.
    std::vector <value_t> bin_sum_x_; 	//! Sum of value in bin
    std::vector <value_t> bin_range_; 	//! Range for bins.
    value_t lower_bound_; 	        //! Minimum value.
    value_t upper_bound_;	        //! Maximum value.
    std::size_t  bin_count_;	//! Number of vaild intervals betwen lower and upper bound
    
    value_t min_inserted_sample_; 	        //! Minimum value.
    value_t max_inserted_sample_;	        //! Maximum value.

    /// incrementally updated values
    value_t sample_count_;	//! Number of accumulated samples.
    value_t sum_x_;	//! Average.
    value_t sum_x2_; 	//! Root mean square.

  protected:
    
    /** 
     * Returns the index of the bin which contains a given value.
     */
    inline std::size_t bin_index(value_t x) const {
      typename std::vector<value_t>::const_iterator it = std::lower_bound(bin_range_.begin(),
									  bin_range_.end(),
									  x); // This is a form of binary search
      assert(it!=bin_range_.begin());
      assert(it!=bin_range_.end());
      assert((*it)>=val);
      int pos = it-bin_range_.begin();
      assert(pos >=1);
      pos -= 1; 
      assert (bin_range_[pos] < val);
      assert (val <= bin_range_[pos+1]);
      return std::size_t(pos);
    }
    
  public:

    histogram() {
      clear();
    }

  public:

    void store_to(output_serializer& s) const {
      s <<
	bin_count_x_ << bin_sum_x_ << bin_range_ <<
	lower_bound_ << upper_bound_ <<
	sl::uint32_t(bin_count_) <<
	min_inserted_sample_ <<
	max_inserted_sample_ <<
	sample_count_ <<
	sum_x_ <<
	sum_x2_;
    }
    
    void retrieve_from(input_serializer& s) {
      sl::uint32_t n;
      s >>
	bin_count_x_ >> bin_sum_x_ >> bin_range_ >>
	lower_bound_ >> upper_bound_ >>
	n >>
	min_inserted_sample_ >>
	max_inserted_sample_ >>
	sample_count_ >>
	sum_x_ >>
	sum_x2_;
      bin_count_ = n;
    }
    
  public:
    
    void clear() {
      bin_count_x_.clear();
      bin_sum_x_.clear();
      bin_range_.clear();
      bin_count_=0;
      lower_bound_=0;
      upper_bound_=1;
      
      reset();
    }

    void reset() {
      std::fill(bin_count_x_.begin(), bin_count_x_.end(), 0);
      std::fill(bin_sum_x_.begin(), bin_sum_x_.end(), 0);
      sample_count_=0;
      sum_x_=0;
      sum_x2_=0;
      min_inserted_sample_ = sl::scalar_math<value_t>::finite_upper_bound();
      max_inserted_sample_ = sl::scalar_math<value_t>::finite_lower_bound();
    }
    
    /** 
     * Init ranges to with a (i/n)^gamma distribution
     */
    void set_range(value_t l, value_t u, std::size_t n, value_t gamma=1.0) {
      clear();
	
      lower_bound_=l;
      upper_bound_= u;
      bin_count_=n;

      bin_count_x_.resize(n+2);
      std::fill(bin_count_x_.begin(),bin_count_x_.end(),0);
      bin_sum_x_.resize(n+2);
      std::fill(bin_sum_x_.begin(),bin_sum_x_.end(),0);
      bin_range_.resize(n+3);
	
      bin_range_[0]   = sl::scalar_math<value_t>::finite_lower_bound(); 
      bin_range_[n+2] = sl::scalar_math<value_t>::finite_upper_bound(); 

      double delta=(upper_bound_-lower_bound_);
      if (gamma==1.0) {
	for(std::size_t i=0; i<=n; ++i) {
	  bin_range_[i+1] = lower_bound_ + delta*value_t(i)/value_t(n);
	}
      } else{
	for (std::size_t i=0; i<=n; ++i) {
	  bin_range_[i+1] = lower_bound_ + delta*std::pow(value_t(i)/value_t(n),gamma);
	}
      }
    }

    inline value_t lower_bound() const {return lower_bound_; }
    inline value_t upper_bound() const {return upper_bound_; }
    
    inline value_t min_inserted_sample() const {return min_inserted_sample_; }
    inline value_t max_inserted_sample() const {return max_inserted_sample_; }
      
    inline void insert(value_t x, value_t multiplicity=value_t(1.0)) {
      std::size_t pos=bin_index(x);
      if (x<min_inserted_sample_) min_inserted_sample_=x;
      if (x>max_inserted_sample_) max_inserted_sample_=x;
      if (pos>=0 && pos<=bin_count_+1) { 
	bin_count_x_[pos]+=multiplicity;
	bin_sum_x_[pos]+=x*multiplicity;
	sample_count_ += multiplicity;
	sum_x_  += x*multiplicity;
	sum_x2_ += (x*x)*multiplicity;
      }
    }

    inline void insert(const this_t& other) {
      // Insert into this all the elements of other
      if (!other.is_empty()) {
	for (std::size_t i=0; i<=other.bin_count_+1; ++i) {
	  // We assume that each bin contains the same value
	  // with multiplicity
	  value_t ci = other.bin_count_x_[i];
	  if (ci) {
	    value_t si = other.bin_sum_x_[i];
	    insert(si/ci, ci);
	  }
	}
	// Insert extrema with zero multiplicity, just to
	// update bounds
	insert(other.min_inserted_sample_, 0);
	insert(other.max_inserted_sample_, 0);
      }
    }
		       
    inline std::size_t bin_count() const { return bin_count_; }
    inline bool        is_empty() const { return bin_count_==0; }
    
    inline value_t bin_lower_bound(int index) const {return bin_range_[index]; }
    inline value_t bin_upper_bound(int index) const {return bin_range_[index+1]; }

    inline std::size_t sample_count() const {
      return sample_count_;
    }
    
    inline value_t max_sample_count_in_bin() const {
      return *(std::max_element(bin_count_x_.begin(),bin_count_x_.end()));
    }

    inline value_t sample_count_in_bin(int idx) const {
      return bin_count_x_[idx];
    }
            
    inline value_t sample_count_in_range(value_t l, value_t u) const {
      std::size_t first_bin= bin_index(l);
      std::size_t last_bin = bin_index(u);
      value_t result=0;
      for (std::size_t i=first_bin; i<=last_bin; ++i) {
	result += bin_count_x_[i];
      }
      return result;
    }

    inline value_t sample_count_in_bin_of(value_t v) const {
      return bin_count_x_[bin_index(v)];
    }
    
    inline value_t bin_width_of(value_t v) {
      std::size_t pos = bin_index(v);
      return bin_range_[pos+1]-bin_range_[pos];
    }
    
    value_t percentile(value_t frac) const {
      assert(frac >= 0 && frac <= 1);
      value_t result = 0;
      if (sample_count_>0) {
	const std::size_t N = bin_count_x_.size();

	value_t total_sum = 0;
	for (std::size_t i=0; i<N; ++i) {
	  total_sum+= bin_count_x_[i];
	}
	total_sum*=frac;
	
	value_t partial_sum=0;
	for (std::size_t i=0; i<N; ++i) {
	  partial_sum+= bin_count_x_[i];
	  if (partial_sum>=total_sum) {
	    result = std::min(max_inserted_sample_, bin_range_[i+1]);
	    break;
	  }
	}
      }	
      return result;
    }
        
    //! Returns the average of the data.
    inline value_t average() const {
      return sum_x_/sample_count_;
    }
	
    //! Returns the Root Mean Square of the data.
    inline value_t RMS() const {
      return value_t(std::sqrt(double(sum_x2_)/double(sample_count_)));
    }
	
    //! Returns the variance of the data.
    inline value_t variance() const {
      double avg = double(average());
      return value_t(std::abs(double(sum_x2_)/double(sample_count_)-avg*avg));
    }

    //! Returns the standard deviation of the data.
    inline value_t standard_deviation() const {
      return std::sqrt(variance());
    }
  };

  typedef histogram<float>  histogramf;
  typedef histogram<double> histogramd;
  
} // namespace sl

#endif
