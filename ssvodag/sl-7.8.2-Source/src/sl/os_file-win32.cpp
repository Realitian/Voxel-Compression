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
#include <sl/os_file.hpp>
#include <iostream>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#include <stdio.h>
#include <windows.h>

namespace sl {

  os_file::file_handle_t os_file::file_open(const char* file_name,
                                            data_access_flags_t   access_flags,
                                            file_creation_flags_t creation_flags) {
    // get file handle
    DWORD desired_access;
    DWORD share_mode;
    DWORD creation_disposition;
    switch ( access_flags ) {
    case OS_READ_ONLY:
      desired_access = GENERIC_READ;
      share_mode = FILE_SHARE_READ;
      break;
    case OS_WRITE_ONLY:
      desired_access = GENERIC_WRITE;
      share_mode = FILE_SHARE_WRITE;
      break;
    case OS_READ_WRITE:
      desired_access = GENERIC_READ | GENERIC_WRITE;
      share_mode = FILE_SHARE_READ | FILE_SHARE_WRITE;
      break;
    default:
      // std::cerr << "os_file::open: wrong flag opening " << file_name << std::endl;
      return -1;
      break;
    }
    switch (creation_flags) {
    case OS_OPEN_EXISTING:
      creation_disposition = OPEN_EXISTING;
    break;
    case OS_OPEN_CREATE_IF_NOT_PRESENT:
      creation_disposition = OPEN_ALWAYS;
      break;
    }
    // CreateFile take as input an array of short: convert file_name
    int str_len = strlen( file_name );
    char* str = new char [ str_len+1 ];
    for(int i = 0; i <= str_len; i++ )
      str[ i ] = file_name[ i ];

    // open the file and get an handle to it.
    file_handle_t file_descriptor = (file_handle_t) ::CreateFile( str, desired_access, 
                                                                share_mode, NULL, creation_disposition,
                                                                FILE_ATTRIBUTE_NORMAL | FILE_FLAG_RANDOM_ACCESS, NULL );

    delete[] str;
    return file_descriptor;
  }
 
  void os_file::file_close(file_handle_t file_descriptor) {
    // close file
    CloseHandle((HANDLE)file_descriptor);
  }

  void os_file::file_delete(const char* file_name) { 
    // DeleteFile take as input an array of short: convert file_name
    int str_len = strlen( file_name );
    char* str = new char [ str_len + 1];
    for(int i = 0; i <= str_len; i++ )
      str[ i ] = file_name[ i ];

    ::DeleteFile( str );
    delete[] str;
  }

  void os_file::file_resize(file_handle_t fd, uint64_t len) {
    LARGE_INTEGER qlen;
    qlen.QuadPart = len;
    SetFilePointerEx(((HANDLE)fd), qlen, NULL, FILE_BEGIN);
    SetEndOfFile(((HANDLE)fd));
  }

  uint64_t os_file::file_size(file_handle_t file_descriptor) { 
    LARGE_INTEGER qsz;
    GetFileSizeEx((HANDLE)file_descriptor, &qsz);
    return qsz.QuadPart;
  }
  
  size_t os_file::memory_page_size() { 
    SYSTEM_INFO sys_info;
    GetSystemInfo( &sys_info );
    return sys_info.dwAllocationGranularity;
  }

  void* os_file::memory_map(size_t len, data_access_flags_t prot, file_handle_t file_descriptor, uint64_t offset) { 
    // set flags
    DWORD protect;
    DWORD desired_access;
    switch( prot ) {
    case OS_READ_ONLY :
      protect = PAGE_READONLY;
      desired_access = FILE_MAP_READ;
      break;
    case OS_WRITE_ONLY:
      protect = PAGE_READWRITE;
      desired_access = FILE_MAP_WRITE;
      break;
    case OS_READ_WRITE:
      protect = PAGE_READWRITE;
      desired_access = FILE_MAP_ALL_ACCESS;
      break;
    default:
      SL_TRACE_OUT(-1) << "os_file::mmap: wrong flag\n";
      return 0;
      break;
    }

    // get handle for memory mapping
    // FIXME Hi/Lo    
    LARGE_INTEGER qlen;
    qlen.QuadPart = len+offset;
    int mapping = (int) CreateFileMapping((HANDLE)file_descriptor, NULL, 
					  protect, qlen.HighPart, qlen.LowPart, NULL);
    if ( mapping == -1 ) {
      SL_TRACE_OUT(-1) << "os_file::mmap: CreateFileMapping() failed\n";
      return 0;
    }

    // memory map offset is assumed to be a double word int
    // offset must be page aligned.
    void* data;
    // FIXME Hi/Lo
    LARGE_INTEGER qoffset;
    qoffset.QuadPart = offset;
    if (!(data = MapViewOfFile((HANDLE)mapping, desired_access, qoffset.HighPart, qoffset.LowPart, len))) {
      SL_TRACE_OUT(-1) << "os_file::mmap: MapViewOfFile() failed\n";
      CloseHandle((HANDLE)mapping);
      return 0;
    }

    CloseHandle((HANDLE)mapping);

    return data;
  }

  void os_file::memory_unmap(void* addr,  size_t /* len */ ) { 
    // len parameter is ignored on win32.
    UnmapViewOfFile( addr );
  }


  void os_file::memory_advise_willneed(const void* addr, size_t len) { 
    // madvise function not implemented on windows
  }

  void os_file::memory_advise_dontneed(const void* addr, size_t len) { 
    // madvise function not implemented on windows
  }

  void os_file::memory_incore(const void* addr, size_t len,unsigned char *vec) {
    // mincore function not implemented on windows  
    int num_bytes = len / memory_page_size();
    for( int i = 0; i < num_bytes; i++ ) {
      vec[ i ] = 1;
    }
  }

} // os namespace
