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
#ifndef OS_FILE_HPP
#define OS_FILE_HPP

#include <sl/cstdint.hpp>

namespace sl {
  
  /**
   * Operating system independent functions for file access and
   * memory mapping.
   */
  class os_file {
  public: 
    typedef enum data_access_flags {
      OS_READ_ONLY = 0, OS_WRITE_ONLY = 1, OS_READ_WRITE = 2
    } data_access_flags_t;

    typedef enum file_creation_flags {
      OS_OPEN_EXISTING = 0, OS_OPEN_CREATE_IF_NOT_PRESENT = 1
    } file_creation_flags_t;

    typedef int64_t file_handle_t;
    
  public:

    /**
     * Open file file_name. flags can be OS_RDONLY, OS_WRONLY, OS_RDWR.
     * Return a file_descriptor (HANDLE under win32).
     * Return -1 if fail.
     */
    static file_handle_t file_open(const char* file_name,
                                   data_access_flags_t   access_flags,
                                   file_creation_flags_t creation_flags);

    /**
     * Close the file identified by file_descriptor, returned by
     * open function.
     */
    static void file_close(file_handle_t file_descriptor);

    /**
     * Delete file file_name. 
     */
    static void file_delete(const char* file_name);

    /**
     * Return file length in bytes of the file identified 
     * by file_descriptor, returned by open function.
     */ 
    static uint64_t file_size(file_handle_t file_descriptor);

    /**
     * Truncate or extend file to specific length
     */ 
    static void file_resize(file_handle_t file_descriptor, uint64_t len);

    /**
     * Return operating system memory page size (granularity).
     */
    static size_t memory_page_size();

    /**
     * Memory map a region of the file (file_descriptor) , 
     * starting from offset for len bytes; flags can be
     * OS_RDONLY, OS_WRONLY, OS_RDWR.
     * Offset must be page_size aligned.
     * Return a void* to the first mmapped byte.
     */
    static void* memory_map(size_t len, data_access_flags_t flags, file_handle_t file_descriptor, uint64_t offset);

    /**
     * Unmap a pointer obtained by mmap.
     */
    static void memory_unmap(void* addr, size_t len);

    /*
     * Advise an area of memory. Implemented only in unix.
     * Addr is the pointer to the area, len is its length.
     */
    static void memory_advise_willneed(const void* addr, size_t len);

    /*
     * Advise an area of memory. Implemented only in unix.
     * Addr is the pointer to the area, len is its length.
     */
    static void memory_advise_dontneed(const void* addr, size_t len);

    /**
     * Check if an area of memory is incore.
     * Addr is the pointer to the area, len is its length.
     * On win32 return always true
     */
    static void memory_incore(const void* addr, size_t len, unsigned char *vec);

  }; // class os_file


} // namespace sl

#endif // OS_FILE_HPP
