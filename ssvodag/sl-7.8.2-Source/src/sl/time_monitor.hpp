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

#ifndef SL_TIME_MONITOR_HPP 
#define SL_TIME_MONITOR_HPP

#include <sl/time_duration.hpp>
#include <sl/clock.hpp>
#include <map>
#include <vector>

namespace sl {
  
  /**
   *  Objects that monitor cpu and real time
   */
  class time_monitor {
  public:
    typedef std::pair<time_duration, time_duration> rt_cpu_duration_t;
    
  protected:
    real_time_clock rt_clock_;
    cpu_time_clock  cpu_clock_;
    
    std::map<std::size_t, std::string> tag_to_string_;
    std::vector< std::pair<std::size_t, rt_cpu_duration_t> > stack_;
    std::map<std::size_t, rt_cpu_duration_t> tag_to_rt_cpu_duration_;

    rt_cpu_duration_t level_0_duration_;

    rt_cpu_duration_t half_resolution_;
    
  public:

    time_monitor();

    ~time_monitor();
    
    void insert_category(std::size_t tag, const std::string& description);

    const std::string& category_tag(std::size_t tag) const;

    std::size_t category_count() const;
    
    const std::map<std::size_t, std::string>& category_map() const;
    
    void clear();

    void restart();
    
    void begin(std::size_t tag);

    void end();

    const time_duration& real_time_duration(std::size_t tag) const;

    const time_duration& cpu_time_duration(std::size_t tag) const;

    const time_duration& real_time_duration() const;
    
    const time_duration& cpu_time_duration() const;
    
  };

} // namespace sl

extern std::ostream& operator<< (std::ostream& os, const sl::time_monitor& t);
 
#endif
