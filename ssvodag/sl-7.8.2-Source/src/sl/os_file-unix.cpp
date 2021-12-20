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
/// --------- Bgn Large files on 32 bit machines
#undef _FILE_OFFSET_BITS
#undef _LARGEFILE_SOURCE
#undef _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
/// --------- End Large files on 32 bit machines
#
#include <sl/os_file.hpp>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>

#ifndef O_LARGEFILE
#define O_LARGEFILE 0
#endif

namespace sl {
  
  os_file::file_handle_t os_file::file_open(const char* file_name,
                                            data_access_flags_t   access_flags,
                                            file_creation_flags_t creation_flags) {
    int desired_access;
    mode_t share_mode;
    switch (access_flags) {
    case OS_READ_ONLY :
      desired_access = O_RDONLY | O_LARGEFILE;
      share_mode = S_IRUSR;
      break;
    case OS_WRITE_ONLY:
      desired_access = O_WRONLY | O_LARGEFILE;
      share_mode = S_IWUSR;
      break;
    case OS_READ_WRITE:
      desired_access = O_RDWR | O_LARGEFILE;
      share_mode = S_IRUSR | S_IWUSR;
      break;
    default:
      // std::cerr << "os_file::open: wrong flag opening " << file_name << std::endl;
      return -1;
      break;
    }
    switch (creation_flags) {
    case OS_OPEN_EXISTING:
      break;
    case OS_OPEN_CREATE_IF_NOT_PRESENT:
      desired_access |= O_CREAT;
      break;
    }
    return ::open( file_name, desired_access, share_mode); // -1 if an error occurred
  }

  void os_file::file_close(file_handle_t file_descriptor) { 
    ::close( (int)file_descriptor );
  }

  void os_file::file_delete(const char* file_name) { 
    ::unlink( file_name );
  }

  void os_file::file_resize(file_handle_t fd, uint64_t len) {
    uint64_t sz = file_size((int)fd);
    if (len<sz) {
      int result=::ftruncate((int)fd,len);
      SL_USEVAR(result);
    } else if (len>sz) {
      char zero='\0';
      ::lseek((int)fd, len-1, SEEK_SET);
      int result=::write((int)fd, &zero, 1);
      SL_USEVAR(result);
    }
  }
  
  uint64_t os_file::file_size(file_handle_t file_descriptor) { 
    off_t result = ::lseek((int)file_descriptor, 0L, SEEK_END);	      
    if ( result == (off_t)-1 ) {
      // std::cerr << "os_file::file_length failed\n";
    }
    return result;
  }

  size_t os_file::memory_page_size() { 
    return getpagesize();
  }

  void* os_file::memory_map(size_t len, data_access_flags_t flags, file_handle_t file_descriptor, uint64_t offset) { 
    //set proper protection flag
    int prot = 0;    
    switch ( flags ) {
    case OS_READ_ONLY :
      prot = PROT_READ;
      break;
    case OS_WRITE_ONLY:
      prot = PROT_WRITE;
      break;
    case OS_READ_WRITE:
      prot = PROT_WRITE | PROT_READ;
      break;
    default:
      // std::cerr << "os_file::mmap: wrong flag\n";
      return 0;
      break;
    }
    
    // extend file for writing
    if (flags != OS_READ_ONLY && file_size(file_descriptor)<len) {
      file_resize(file_descriptor, len);
    }

    // memory map the file
    errno = 0;
    void *data = ::mmap(0, len, prot, MAP_SHARED, (int)file_descriptor, offset);
      
    if (data == MAP_FAILED) {
      std::cerr << "os_file::mmap failed: [" << data << "] " << strerror( errno ) << std::endl;
      std::cerr << "file_descriptor = " << (int)file_descriptor << std::endl;
      std::cerr << "len = " << len << std::endl;
      std::cerr << "offset = " << offset << std::endl;
      data = 0;
    }
    return data;
  }

  void os_file::memory_unmap(void* addr, size_t len) { 
    errno = 0;
    int ret = ::munmap(addr, len);
    if (ret == -1) {
      std::cerr << "os_file::munmap failed: [" << ret << "] " << strerror( errno ) << std::endl;
      std::cerr << "len = " << len << std::endl;
    }      
  }

  void os_file::memory_advise_willneed(const void* addr, size_t len) { 
    ::madvise( (void*)addr, len, MADV_WILLNEED );
  }

  void os_file::memory_advise_dontneed(const void* addr, size_t len) { 
    ::madvise( (void*)addr, len, MADV_DONTNEED );
  }

  void os_file::memory_incore(const void* addr, size_t len, unsigned char *vec) {
#if defined(__APPLE__) && defined(__MACH__)
    ::mincore( (void*)addr, len, (char*)vec );
#else
    ::mincore( (void*)addr, len, vec );
#endif
  }

} // os namespace
