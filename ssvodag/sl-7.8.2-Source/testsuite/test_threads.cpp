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
#include <sl/tester.hpp>

#if !SL_HAVE_THREADS
int main() {
  std::cerr << "No threads on this platform." << std::endl;
  return 0;
}
#else

////////////////////////////////////////////////////////////////////

#include <sl/thread.hpp>
#include <list>
//--------------------

static std::size_t failed_test_count = 0;

// Mutex + global count variable
sl::mutex global_mutex;
int       global_count;

// Condition variable
sl::wait_condition global_condition;

class test_id_thread: public sl::thread {
public:
  test_id_thread(){}
  virtual void run() {
    std::cout << " My thread id is " << "???" << std::endl;
  }
};

class test_lock_thread: public sl::thread {
public:
  test_lock_thread(){}
  virtual void run() {
    for (int i = 0; i < 1000; ++i) {
      sl::mutex_locker lock(&global_mutex);
      ++global_count;
    }
  }
};

class test_condition_thread_1: public sl::thread {
public:
  test_condition_thread_1(){}
  virtual void run() {
    sl::mutex_locker lock(&global_mutex);
    --global_count;
    global_condition.notify_all();
  }
};

class test_condition_thread_2: public sl::thread {
public:
  test_condition_thread_2(){}
  virtual void run() {
    std::cout << " Waiting..." << std::flush;
    sl::mutex_locker lock(&global_mutex);
    while (global_count > 0) {
      std::cout << "." << std::flush;
      global_condition.wait(global_mutex);
    }
    std::cout << "." << std::endl;
  }
};

class test_yield_thread: public sl::thread {
public:
  test_yield_thread(){}
  virtual void run() {
    yield_current_thread();
  }
};


// This is the main program (i.e. the main thread)
int main() {
  sl::tester tester("Threading subsystem");

  std::cerr << std::endl << "PART I: Hardware queries" << std::endl;
  tester.unchecked_test("Number of processor cores", sl::thread::hardware_concurrency());
  
  // Test 2: thread IDs
  std::cerr << std::endl << "PART II: Thread IDs" << std::endl;
  {
    // Show the main thread ID
    // std::cout << " Main thread id is " << this_thread::get_id() << "." << std::endl;

    // Start a bunch of child threads - only run a single thread at a time
    test_id_thread t1;
    t1.start();
    t1.wait();
    test_id_thread t2;
    t2.start();
    t2.wait();
    test_id_thread t3;
    t3.start();
    t3.wait();
  }

  // Test 3: mutex locking
  std::cerr << std::endl << "PART III: Mutex locking (100 threads x 1000 iterations)" << std::endl;
  {
    // Clear the global counter.
    global_count = 0;

    // Start a bunch of child threads
    std::list<sl::thread*> thread_list;
    for(int i = 0; i < 100; ++ i) {
      sl::thread* t = new test_lock_thread;
      t->start();
      thread_list.push_back(t);
    }

    // Wait for the threads to finish
    for(std::list<sl::thread*>::iterator it = thread_list.begin();
	it != thread_list.end(); ++it) {
      sl::thread* t = *it;
      t->wait();
      delete t;
    }

    // Check the global count
    tester.test("Parallel inc", global_count, 100*1000);
  }

  // Test 4: condition variable
  std::cout << std::endl << "PART IV: Condition variable (40 + 1 threads)" << std::endl;
  {
    // Set the global counter to the number of threads to run.
    global_count = 40;

    // Start the waiting thread (it will wait for gCount to reach zero).
    test_condition_thread_2 t1;
    t1.start();

    // Start a bunch of child threads (these will decrease gCount by 1 when they
    // finish)
    std::list<sl::thread*> thread_list;
    for(int i = 0; i < 40; ++ i) {
      sl::thread* t = new test_condition_thread_1;
      t->start();
      thread_list.push_back(t);
    }

    // Wait for the waiting thread to finish
    t1.wait();

   // Wait for the other threads to finish
    for(std::list<sl::thread*>::iterator it = thread_list.begin();
	it != thread_list.end(); ++it) {
      sl::thread* t = *it;
      t->wait();
      delete t;
    }
    tester.test("Parallel dec", global_count, 0);
  }

  // Test 5: yield
  std::cout << std::endl << "PART V: Yield (40 + 1 threads)" << std::endl;
  {
    // Start a bunch of child threads
    std::list<sl::thread*> thread_list;
    for(int i = 0; i < 40; ++ i) {
      sl::thread* t = new test_yield_thread;
      t->start();
      thread_list.push_back(t);
    }

    // Yield...
    sl::thread::yield_current_thread();

    // Wait for the threads to finish
    for(std::list<sl::thread*>::iterator it = thread_list.begin();
	it != thread_list.end(); ++it) {
      sl::thread* t = *it;
      t->wait();
      delete t;
    }
  }

  failed_test_count += tester.failed_test_count();
  return failed_test_count;
}

#endif // SL_HAVE_THREADS
