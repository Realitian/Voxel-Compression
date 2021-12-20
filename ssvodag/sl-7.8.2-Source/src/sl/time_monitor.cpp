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
#include <sl/time_monitor.hpp>
#include <sl/utility.hpp>
#include <iomanip>

namespace sl {

  time_monitor::time_monitor() {
    half_resolution_.first  = rt_clock_.resolution();
    half_resolution_.second = cpu_clock_.resolution();
    for (std::size_t i=0; i<10; ++i) {
      half_resolution_.first = std::min(half_resolution_.first, rt_clock_.resolution());
      half_resolution_.second = std::min(half_resolution_.second, cpu_clock_.resolution());
    }
    half_resolution_.first  /= 2;
    half_resolution_.second /= 2;
  }
  
  time_monitor::~time_monitor() {}
  
  void time_monitor::insert_category(std::size_t tag, const std::string& description) {
    tag_to_string_[tag] = description;
  }
  
  const std::string& time_monitor::category_tag(std::size_t tag) const {
    return const_cast<time_monitor*>(this)->tag_to_string_[tag];
  }
  
  std::size_t time_monitor::category_count() const {
    return tag_to_string_.size();
  }
  
  const std::map<std::size_t, std::string>& time_monitor::category_map() const {
    return tag_to_string_;
  }
  
  void time_monitor::clear() {
    tag_to_string_.clear();
    restart();
  }
  
  void time_monitor::restart() {
    cpu_clock_.restart();
    rt_clock_.restart();
    tag_to_rt_cpu_duration_.clear();
    stack_.clear();
    level_0_duration_ = rt_cpu_duration_t();
  }
  
  void time_monitor::begin(std::size_t tag) {
    std::pair<std::size_t, rt_cpu_duration_t> tag_start;
    tag_start.first = tag;
    tag_start.second.first  = rt_clock_.elapsed();
    tag_start.second.second = cpu_clock_.elapsed();
    stack_.push_back(tag_start);
  }
  
  void time_monitor::end() {
    std::pair<std::size_t, rt_cpu_duration_t> tag_start = stack_.back();
    stack_.pop_back();
    rt_cpu_duration_t& old_duration = tag_to_rt_cpu_duration_[tag_start.first];
    rt_cpu_duration_t  delta_duration = std::make_pair(rt_clock_.elapsed()-tag_start.second.first+half_resolution_.first,
						       cpu_clock_.elapsed()-tag_start.second.second+half_resolution_.second);
    delta_duration.second = std::min(delta_duration.first,delta_duration.second); // enforce cpu time <= real time 
    old_duration.first  += delta_duration.first;
    old_duration.second += delta_duration.second;
    if (stack_.empty()) {
      level_0_duration_.first += delta_duration.first;
      level_0_duration_.second += delta_duration.second;
    }
  }
  
  const time_duration& time_monitor::real_time_duration(std::size_t tag) const {
    return const_cast<time_monitor*>(this)->tag_to_rt_cpu_duration_[tag].first;
  }
  
  const time_duration& time_monitor::cpu_time_duration(std::size_t tag) const {
    return const_cast<time_monitor*>(this)->tag_to_rt_cpu_duration_[tag].second;
  }
  
  const time_duration& time_monitor::real_time_duration() const {
    return level_0_duration_.first;
  }
  
  const time_duration& time_monitor::cpu_time_duration() const {
    return level_0_duration_.second;
  }

} // namespace sl

std::ostream& operator<< (std::ostream& os, const sl::time_monitor& t) {
  os << std::left << std::setw(36) << "TAG" << std::right << std::setw(12) << "REALTIME (s)" << std::setw(12) << "CPUTIME (s)" << std::setw(12) << "%CPU" << std::endl;
  sl::time_duration t_rt_all = t.real_time_duration();
  sl::time_duration t_cpu_all = t.cpu_time_duration();
  
  for (std::map<std::size_t, std::string>::const_iterator it = t.category_map().begin();
       it != t.category_map().end();
       ++it) {
    sl::time_duration t_rt  = t.real_time_duration(it->first);
    sl::time_duration t_cpu = t.cpu_time_duration(it->first);
    double cpu_percent = sl::median(0.0, 100.0,
				    100.0*double(t_cpu.as_microseconds())/double(t_rt.as_microseconds()+1));
    os << std::left << std::setw(36) << it->second << std::right << std::setw(12) << t_rt.as_seconds() << std::setw(12) << t_cpu.as_seconds() << std::setw(12) << (int)cpu_percent << std::endl;
  }
  double cpu_percent_all = sl::median(0.0, 100.0,
				      100.0*double(t_cpu_all.as_microseconds())/double(t_rt_all.as_microseconds()+1));
  os << std::left << std::setw(36) << "TOTAL" << std::right << std::setw(12) << t_rt_all.as_seconds() << std::setw(12) << t_cpu_all.as_seconds() << std::setw(12) << (int)cpu_percent_all << std::endl;
  
  return os;
}

