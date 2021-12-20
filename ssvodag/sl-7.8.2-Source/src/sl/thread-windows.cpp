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
#include <windows.h>

namespace sl {

  namespace detail {

    //--------------------------------------------------------------------
    // Mutex
    //--------------------------------------------------------------------

    // NOTE: The Win32 implementation is based on critical section for
    // recursive mutexes and on the  Interlocked* mechanism for non-recursive
    // ones. The latter method is the one implemented in win32-pthreads.

    class mutex {
      SL_DISABLE_COPY(mutex);
    protected:
      friend class wait_condition;

      bool is_recursive_;
      union {
	struct {
	  CRITICAL_SECTION cs_;
	} as_recursive_;
	struct {
	  LONG lock_idx_;             /* Provides exclusive access to mutex state
				       via the Interlocked* mechanism.
				       0: unlocked/free.
				       1: locked - no other waiters.
				       -1: locked - with possible other waiters.
				      */
	  HANDLE event_;              /* Mutex release notification to waiting
					 threads. */
	} as_nonrecursive_;
      } primitive_;

    public:
      
      inline mutex(bool is_recursive = false): is_recursive_(is_recursive) {
	if (is_recursive_) {
	  InitializeCriticalSection(&primitive_.as_recursive_.cs_);
	} else {
	  primitive_.as_nonrecursive_.lock_idx_ = 0;
	  primitive_.as_nonrecursive_.event_ = CreateEvent (NULL, 0,    /* manual reset = No */
							   0,          /* initial state = not signaled */
							   NULL);      /* event name */
	}
      }

      inline ~mutex() {
	// FIXME: Check whether locked by someone???
	if (is_recursive_) {
	  DeleteCriticalSection(&primitive_.as_recursive_.cs_);
	} else {
	  CloseHandle(primitive_.as_nonrecursive_.event_);
	}
      }

      inline void lock() {
	if (is_recursive_) {
	  EnterCriticalSection(&primitive_.as_recursive_.cs_);
	} else {
	  if ((LONG)InterlockedExchange((LPLONG)&(primitive_.as_nonrecursive_.lock_idx_), (LONG)1) != (LONG)0) {
	    while ((LONG)InterlockedExchange((LPLONG)&(primitive_.as_nonrecursive_.lock_idx_), (LONG)-1) != (LONG)0) {
	      if (WAIT_OBJECT_0 != WaitForSingleObject(primitive_.as_nonrecursive_.event_, INFINITE)) {
		// FIXME result = EINVAL;
		// FIXME break;
	      }
	    }
	  }
	}
      }

      inline bool try_lock() {
	if (is_recursive_) {
	  return TryEnterCriticalSection(&primitive_.as_recursive_.cs_) ? true : false;
	} else {
	  if ((LONG)0 == (LONG) InterlockedCompareExchange((LPLONG) &primitive_.as_nonrecursive_.lock_idx_,
						     (LONG) 1,
						     (LONG) 0)) {
	    // Lock acquired
	    return true;
	  } else {
	    // Lock busy
	    return false;
	  }
	}
      }
      
      inline void unlock() {
	if (is_recursive_) {
	  LeaveCriticalSection(&primitive_.as_recursive_.cs_);
	} else {
	  LONG idx = InterlockedExchange((LPLONG) &primitive_.as_nonrecursive_.lock_idx_, (LONG) 0);
	  if (idx != 0) {
	    if (idx < 0) {
	      /*
	       * Someone may be waiting on that mutex.
	       */
	      if (SetEvent (primitive_.as_nonrecursive_.event_) == 0) {
		//result = EINVAL;
	      }
	    }
	  } else {
	    /*
	     * Was not locked (so can't be owned by us).
	     */
	    // result = EPERM;
	  }
	}
      }
    }; // class mutex

    //--------------------------------------------------------------------
    // Wait condition
    //--------------------------------------------------------------------

    // NOTE: The Win32 implementation is based on the corresponding
    // implementation in GLFW and TinyThreads, in turn based on a
    // description by Douglas C. Schmidt and Irfan Pyarali:
    // http://www.cs.wustl.edu/~schmidt/win32-cv-1.html

    class wait_condition {
      SL_DISABLE_COPY(wait_condition);
    protected:
      enum { CONDITION_EVENT_ONE= 0, CONDITION_EVENT_ALL=1};
      
      HANDLE           events_[2];          ///< Signal and broadcast event HANDLEs.
      unsigned int     waiters_count_;      ///< Count of the number of waiters.
      CRITICAL_SECTION waiters_count_lock_; ///< Serialize access to waiters_count_.
    public:
      inline wait_condition() {
	waiters_count_ = 0;
	events_[CONDITION_EVENT_ONE] = CreateEvent(NULL, FALSE, FALSE, NULL);
	events_[CONDITION_EVENT_ALL] = CreateEvent(NULL, TRUE, FALSE, NULL);
	InitializeCriticalSection(&waiters_count_lock_);
      }

      inline ~wait_condition() {
	CloseHandle(events_[CONDITION_EVENT_ONE]);
	CloseHandle(events_[CONDITION_EVENT_ALL]);
	DeleteCriticalSection(&waiters_count_lock_);
      }

      inline void wait(mutex &m) {
	// Increment number of waiters
	EnterCriticalSection(&waiters_count_lock_);
	++ waiters_count_;
	LeaveCriticalSection(&waiters_count_lock_);
	
	// It's ok to release the mutex here since Win32 manual-reset events maintain
	// state when used with SetEvent()
	m.unlock(); // FIXME LeaveCriticalSection(&m.handle_);
	
	// Wait for either event to become signaled due to notify_one() or
	// notify_all() being called
	int result = WaitForMultipleObjects(2, events_, FALSE, INFINITE);
	
	// Check if we are the last waiter
	EnterCriticalSection(&waiters_count_lock_);
	-- waiters_count_;
	bool is_last_waiter =
	  (result == (WAIT_OBJECT_0 + CONDITION_EVENT_ALL)) &&
	  (waiters_count_ == 0);
	LeaveCriticalSection(&waiters_count_lock_);
	
	// If we are the last waiter to be notified to stop waiting, reset the event
	if(is_last_waiter) {
	  ResetEvent(events_[CONDITION_EVENT_ALL]);
	}
	// Reacquire the mutex
	m.lock(); // FIXME EnterCriticalSection(&m.handle_);
      }
      
      inline void notify_one() {
	// Are there any waiters?
	EnterCriticalSection(&waiters_count_lock_);
	bool is_having_waiters = (waiters_count_ > 0);
	LeaveCriticalSection(&waiters_count_lock_);

	// If we have any waiting threads, send them a signal
	if(is_having_waiters) {
	  SetEvent(events_[CONDITION_EVENT_ONE]);
	}
      }

      inline void notify_all() {
	// Are there any waiters?
	EnterCriticalSection(&waiters_count_lock_);
	bool is_having_waiters = (waiters_count_ > 0);
	LeaveCriticalSection(&waiters_count_lock_);

	// If we have any waiting threads, send them a signal
	if(is_having_waiters) {
	  SetEvent(events_[CONDITION_EVENT_ALL]);
	}
      }
    }; // class wait_condition
    
    //--------------------------------------------------------------------
    // Thread
    //--------------------------------------------------------------------

    class thread {
      SL_DISABLE_COPY(thread);
    protected:
      sl::thread* owner_;
      HANDLE handle_; ///< Thread handle.
      DWORD  win32_thread_id_;
      int    win32_priority_;
      
      static DWORD WINAPI _sys_thread_wrapper(LPVOID self) {
	thread* t = (thread*)self;
	t->run();
	return (0);
      }
      
    public:
      
      inline thread(sl::thread* owner) :
	owner_(owner),
	handle_(0), win32_thread_id_(0) {
      }
      
      inline ~thread() {
	if (owner_->is_running()) {
	  // Destroying a still running thread!
	  std::terminate();
	}
      }
      
      inline bool wait() {
	if (owner_->is_running()) {
	  WaitForSingleObject(handle_, INFINITE);
	  return !owner_->is_running(); // FIXME
	} else {
	  return true;
	}
      }

      inline void start(sl::thread::Priority priority,
			std::size_t stack_size) {
	owner_->mutex_.lock();
	owner_->is_running_ = true;
	owner_->is_finished_ = false;

	switch (priority) {
	case sl::thread::IdlePriority:         win32_priority_ = THREAD_PRIORITY_IDLE; break;
	case sl::thread::LowestPriority:       win32_priority_ = THREAD_PRIORITY_LOWEST; break;
	case sl::thread::LowPriority:          win32_priority_ = THREAD_PRIORITY_BELOW_NORMAL; break;
	case sl::thread::NormalPriority:       win32_priority_ = THREAD_PRIORITY_NORMAL; break;
	case sl::thread::HighPriority:         win32_priority_ = THREAD_PRIORITY_ABOVE_NORMAL; break;
	case sl::thread::HighestPriority:      win32_priority_ = THREAD_PRIORITY_HIGHEST; break;
	case sl::thread::TimeCriticalPriority: win32_priority_ = THREAD_PRIORITY_TIME_CRITICAL; break;
	case sl::thread::InheritPriority:      win32_priority_ = GetThreadPriority(GetCurrentThread()); break;
	};
	
	handle_ = CreateThread(0, stack_size, _sys_thread_wrapper, (LPVOID)this, 0, &win32_thread_id_);
	if (handle_ == 0) {
	  owner_->is_running_ = false;
	}
	owner_->mutex_.unlock();
      }
      
      inline void run() {
	SetThreadPriority(handle_, win32_priority_);
	
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
	// FIXME can we do better than that?
	unsigned long msecs = (unsigned long)(usecs/1000);
	if (msecs<1) msecs = 1;
	Sleep(msecs); 
      }

      static inline void yield_current_thread() {
	Sleep(0);
      }

      static inline std::size_t hardware_concurrency() {
	SYSTEM_INFO si;
	GetSystemInfo(&si);
	return (int) si.dwNumberOfProcessors;
      }
      
    }; // class thread
  } // namespace detail
} // namespace sl
