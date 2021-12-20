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
#include <sl/thread.hpp>
#include <pthread.h>
#include <signal.h>
#include <sched.h>
#include <unistd.h>

namespace sl {

  namespace detail {

    //--------------------------------------------------------------------
    // Mutex
    //--------------------------------------------------------------------

    class mutex {
      SL_DISABLE_COPY(mutex);
    protected:
      friend class wait_condition;

      pthread_mutexattr_t attr_;
      pthread_mutex_t handle_;
      bool is_recursive_;
    public:
      
      inline mutex(bool is_recursive = false): is_recursive_(is_recursive) {
	pthread_mutexattr_init(&attr_);
	if (is_recursive_) {
	  pthread_mutexattr_settype(&attr_, PTHREAD_MUTEX_RECURSIVE);
	} else {
	  pthread_mutexattr_settype(&attr_, PTHREAD_MUTEX_NORMAL);
	}
	pthread_mutex_init(&handle_, &attr_);
      }

      inline ~mutex() {
	pthread_mutex_destroy(&handle_);
	pthread_mutexattr_destroy(&attr_);
      }

      inline void lock() {
	pthread_mutex_lock(&handle_);
      }

      inline  bool try_lock() {
	return (pthread_mutex_trylock(&handle_) == 0) ? true : false;
      }
      
      inline void unlock() {
	pthread_mutex_unlock(&handle_);
      }
    }; // class mutex

    //--------------------------------------------------------------------
    // Wait condition
    //--------------------------------------------------------------------

    class wait_condition {
      SL_DISABLE_COPY(wait_condition);
    protected:     
      pthread_cond_t handle_;

    public:
      inline wait_condition() {
	pthread_cond_init(&handle_, NULL);
      }

      inline ~wait_condition() {
	pthread_cond_destroy(&handle_);
      }

      inline void wait(mutex &m) {
	pthread_cond_wait(&handle_, &(m.handle_));
      }
      
      inline void notify_one() {
	pthread_cond_signal(&handle_);
      }

      inline void notify_all() {
	pthread_cond_broadcast(&handle_);
      }
    }; // class wait_condition

    //--------------------------------------------------------------------
    // Thread
    //--------------------------------------------------------------------

    class thread {
      SL_DISABLE_COPY(thread);
    protected:
      sl::thread* owner_;
      pthread_t handle_; ///< Thread handle.
      
      static void * _sys_thread_wrapper(void* self) {
	thread* t = (thread*)self;
	t->run();
		  return 0;
      }
      
    public:
      
      inline thread(sl::thread* owner) :
	owner_(owner),
	handle_(0) {
      }
      
      inline ~thread() {
	if (owner_->is_running()) {
	  // Destroying a still running thread!
	  std::terminate();
	}
      }
      
      inline bool wait() {
	if (owner_->is_running()) {
	  pthread_join(handle_, NULL);
	  return !owner_->is_running();
	} else {
	  return true;
	}
      }

      inline void start(sl::thread::Priority priority,
			std::size_t stack_size) {
	owner_->mutex_.lock();
	owner_->is_running_ = true;
	owner_->is_finished_ = false;
	
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE); // For portability

	// Stack size
	if (stack_size!= 0) pthread_attr_setstacksize (&attr, stack_size);

	// Priority
	if (priority == sl::thread::InheritPriority) {
#ifndef ANDROID
      pthread_attr_setinheritsched(&attr, PTHREAD_INHERIT_SCHED);
#endif// ANDROID
	} else {
	  int sched_policy;
	  if (pthread_attr_getschedpolicy(&attr, &sched_policy) != 0) {
	    // failed to get the scheduling policy, don't bother
	    // setting the priority
	  } else {
	    int prio_min = sched_get_priority_min(sched_policy);
	    int prio_max = sched_get_priority_max(sched_policy);
	    if (prio_min == -1 || prio_max == -1) {
	      // failed to get the scheduling parameters, don't
	      // bother setting the priority
	    } else {
	      // crudely scale our priority enum values to the prio_min/prio_max
	      float t =
		(float(sl::thread::TimeCriticalPriority)-float(sl::thread::IdlePriority)) /
		(float(priority) - float(sl::thread::IdlePriority));
	      int prio = int(float(prio_min) + float(prio_max - prio_min) * t + 0.4f);
	      if (prio<prio_min) prio = prio_min;
	      if (prio>prio_max) prio = prio_max;
	    	      
	      sched_param sp;
	      sp.sched_priority = prio;
#ifndef ANDROID
	      pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);
#endif//ANDROID
	      pthread_attr_setschedparam(&attr, &sp);
	    }
	  }
	}

	// Creation!
	if (pthread_create(&handle_, &attr, _sys_thread_wrapper, (void *) this) != 0) {
	  handle_ = 0;
	  owner_->is_running_ = false;
	}
	pthread_attr_destroy(&attr);
	owner_->mutex_.unlock();
      }

      inline void run() {
	try {
	  owner_->run();
	} catch(...) {
	  // Uncaught exceptions will terminate the application (default behavior
	  // according to the C++0x draft)
	  std::terminate();
	}

	// The thread is no longer executing
	owner_->mutex_.lock();
	owner_->is_running_ = false;
	owner_->is_finished_ = true;
	owner_->mutex_.unlock();
      }

      static inline void usleep(unsigned long usecs) {
	::usleep(usecs);
      }

      static inline void yield_current_thread() {
	::sched_yield();
      }

      static inline std::size_t hardware_concurrency() {
#if defined(_SC_NPROCESSORS_ONLN)
	return (std::size_t) sysconf(_SC_NPROCESSORS_ONLN);
#elif defined(_SC_NPROC_ONLN)
	return (std::size_t) sysconf(_SC_NPROC_ONLN);
#else
	// The standard requires this function to return zero if the number of
	// hardware cores could not be determined.
	return 0;
#endif
      }
      
    }; // class thread
  } // namespace detail
} // namespace sl
