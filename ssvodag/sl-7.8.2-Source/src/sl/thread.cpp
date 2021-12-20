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
#include <sl/config.hpp>
#if SL_HAVE_THREADS

#include <sl/thread.hpp>

#if HAVE_WIN32_THREADS
#  include <sl/thread-windows.cpp>
#elif HAVE_POSIX_THREADS
#  include <sl/thread-posix.cpp>
#else
#  error "Unsupported native thread library"
#endif

namespace sl {

  //--------------------------------------------------------------------
  // Mutex
  //--------------------------------------------------------------------
  
  mutex::mutex(RecursionMode mode) {
	  impl_ = new detail::mutex(mode == Recursive ? true : false);
  }

  mutex::~mutex() {
    delete impl_; impl_ = 0;
  }

  void mutex::lock() {
    impl_->lock();
  }

  bool mutex::try_lock() {
    return impl_->try_lock();
  }

  void mutex::unlock() {
    impl_->unlock();
  }

  //--------------------------------------------------------------------
  // Wait condition
  //--------------------------------------------------------------------

  wait_condition::wait_condition() {
    impl_ = new detail::wait_condition();
  }

  wait_condition::~wait_condition() {
    delete impl_; impl_ = 0;
  }

  void wait_condition::wait(mutex &m) {
    impl_->wait(*m.impl_);
  }
		
  void wait_condition::notify_one() {
    impl_->notify_one();
  }

  void wait_condition::notify_all() {
    impl_->notify_all();
  }

  //--------------------------------------------------------------------
  // Thread
  //--------------------------------------------------------------------

  thread::thread(): impl_(0), stack_size_(0), priority_(InheritPriority), is_running_(false), is_finished_(false) {
    impl_ = new detail::thread(this);
  }
  
  thread::~thread() {
    delete impl_; impl_ = 0;
  }

  bool thread::is_finished() const {
    mutex_.lock();
    bool result = is_finished_;
    mutex_.unlock();
    return result;
  }
  
  bool thread::is_running() const {
    mutex_.lock();
    bool result = is_running_;
    mutex_.unlock();
    return result;
  }
  
  thread::Priority thread::priority() const {
    mutex_.lock();
    Priority result = priority_;
    mutex_.unlock();
    return result;
  }
  
  void thread::set_priority(Priority p) {
    mutex_.lock();
    priority_ = p;
    mutex_.unlock();
  }
  
  std::size_t thread::stack_size() const {
    mutex_.lock();
    std::size_t result = stack_size_;
    mutex_.unlock();
    return result;
  }
  
  void thread::set_stack_size(std::size_t s) {
    mutex_.lock();
    stack_size_ = s;
    mutex_.unlock();
  }

  bool thread::wait() {
    return impl_->wait();
  }

  void thread::start() {
    mutex_.lock();
    bool is_running = is_running_;
    Priority priority = priority_;
    std::size_t stack_size = stack_size_;
    mutex_.unlock();
    if (!is_running) {
      impl_->start(priority, stack_size);
    }
  }

  void thread::usleep(unsigned long usecs) {
    detail::thread::usleep(usecs);
  }

  void thread::yield_current_thread() {
    detail::thread::yield_current_thread();
  }

  std::size_t thread::hardware_concurrency() {
    return detail::thread::hardware_concurrency();
  }
} // namespace sl

#endif // SL_HAVE_THREADS
