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
#ifndef SL_GL_GL_HPP
#define SL_GL_GL_HPP

#include <sl/config.hpp>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN /* somewhate limit Win32 pollution */
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>

// Overloads - this is needed for templates

#define SL_WRAP_GL_PROC_1(fname, code, tp) \
  static inline fname(tp arg1) { fname##code(arg1); } 
#define SL_WRAP_GL_PROC_2(fname, code, tp) \
  static inline fname(tp arg1, tp arg2) { fname##code(arg1,arg2); } 
#define SL_WRAP_GL_PROC_3(fname, code, tp) \
  static inline fname(tp arg1, tp arg2, tp arg3) { fname##code(arg1,arg2,arg3); } 
#define SL_WRAP_GL_PROC_4(fname, code, tp) \
  static inline fname(tp arg1, tp arg2, tp arg3, tp arg4) { fname##code(arg1,arg2,arg3,arg4); }  

SL_WRAP_GL_PROC_2(glVertex, 2d, GLdouble)
SL_WRAP_GL_PROC_2(glVertex, 2f, GLfloat)
SL_WRAP_GL_PROC_2(glVertex, 2i, GLint)
SL_WRAP_GL_PROC_2(glVertex, 2s, GLshort)
  
SL_WRAP_GL_PROC_3(glVertex, 3d, GLdouble)
SL_WRAP_GL_PROC_3(glVertex, 3f, GLfloat)
SL_WRAP_GL_PROC_3(glVertex, 3i, GLint)
SL_WRAP_GL_PROC_3(glVertex, 3s, GLshort)

SL_WRAP_GL_PROC_4(glVertex, 4d, GLdouble)
SL_WRAP_GL_PROC_4(glVertex, 4f, GLfloat)
SL_WRAP_GL_PROC_4(glVertex, 4i, GLint)
SL_WRAP_GL_PROC_4(glVertex, 4s, GLshort)

SL_WRAP_GL_PROC_1(glVertex2, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glVertex2, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glVertex2, iv, const GLint*)
SL_WRAP_GL_PROC_1(glVertex2, sv, const GLshort*)

SL_WRAP_GL_PROC_1(glVertex3, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glVertex3, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glVertex3, iv, const GLint*)
SL_WRAP_GL_PROC_1(glVertex3, sv, const GLshort*)

SL_WRAP_GL_PROC_1(glVertex4, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glVertex4, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glVertex4, iv, const GLint*)
SL_WRAP_GL_PROC_1(glVertex4, sv, const GLshort*)

SL_WRAP_GL_PROC_3(glNormal, 3d, GLdouble)
SL_WRAP_GL_PROC_3(glNormal, 3f, GLfloat)
SL_WRAP_GL_PROC_3(glNormal, 3i, GLint)
SL_WRAP_GL_PROC_3(glNormal, 3s, GLshort)

SL_WRAP_GL_PROC_1(glNormal3, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glNormal3, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glNormal3, iv, const GLint*)
SL_WRAP_GL_PROC_1(glNormal3, sv, const GLshort*)

SL_WRAP_GL_PROC_1(glIndex, d, GLdouble)
SL_WRAP_GL_PROC_1(glIndex, f, GLfloat)
SL_WRAP_GL_PROC_1(glIndex, i, GLint)
SL_WRAP_GL_PROC_1(glIndex, s, GLshort)
SL_WRAP_GL_PROC_1(glIndex, ub, GLubyte) /* 1.1 */

SL_WRAP_GL_PROC_1(glIndex, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glIndex, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glIndex, iv, const GLint*)
SL_WRAP_GL_PROC_1(glIndex, sv, const GLshort*)
SL_WRAP_GL_PROC_1(glIndex, ubv, const GLubyte*) /* 1.1 */

SL_WRAP_GL_PROC_3(glColor, 3d, GLdouble)
SL_WRAP_GL_PROC_3(glColor, 3f, GLfloat)
SL_WRAP_GL_PROC_3(glColor, 3i, GLint)
SL_WRAP_GL_PROC_3(glColor, 3s, GLshort)
SL_WRAP_GL_PROC_3(glColor, 3ub, GLubyte) 
SL_WRAP_GL_PROC_3(glColor, 3ui, GLuint) 
SL_WRAP_GL_PROC_3(glColor, 3us, GLushort) 

SL_WRAP_GL_PROC_4(glColor, 4d, GLdouble)
SL_WRAP_GL_PROC_4(glColor, 4f, GLfloat)
SL_WRAP_GL_PROC_4(glColor, 4i, GLint)
SL_WRAP_GL_PROC_4(glColor, 4s, GLshort)
SL_WRAP_GL_PROC_4(glColor, 4ub, GLubyte) 
SL_WRAP_GL_PROC_4(glColor, 4ui, GLuint) 
SL_WRAP_GL_PROC_4(glColor, 4us, GLushort) 

SL_WRAP_GL_PROC_1(glColor, 3d, const GLdouble*)
SL_WRAP_GL_PROC_1(glColor, 3f, const GLfloat*)
SL_WRAP_GL_PROC_1(glColor, 3i, const GLint*)
SL_WRAP_GL_PROC_1(glColor, 3s, const GLshort*)
SL_WRAP_GL_PROC_1(glColor, 3ub, const GLubyte*) 
SL_WRAP_GL_PROC_1(glColor, 3ui, const GLuint*) 
SL_WRAP_GL_PROC_1(glColor, 3us, const GLushort*) 

SL_WRAP_GL_PROC_1(glColor, 4d, const GLdouble*)
SL_WRAP_GL_PROC_1(glColor, 4f, const GLfloat*)
SL_WRAP_GL_PROC_1(glColor, 4i, const GLint*)
SL_WRAP_GL_PROC_1(glColor, 4s, const GLshort*)
SL_WRAP_GL_PROC_1(glColor, 4ub, const GLubyte*) 
SL_WRAP_GL_PROC_1(glColor, 4ui, const GLuint*) 
SL_WRAP_GL_PROC_1(glColor, 4us, const GLushort*) 

SL_WRAP_GL_PROC_1(glTexCoord, 1d, GLdouble)
SL_WRAP_GL_PROC_1(glTexCoord, 1f, GLfloat)
SL_WRAP_GL_PROC_1(glTexCoord, 1i, GLint)
SL_WRAP_GL_PROC_1(glTexCoord, 1s, GLshort)

SL_WRAP_GL_PROC_2(glTexCoord, 2d, GLdouble)
SL_WRAP_GL_PROC_2(glTexCoord, 2f, GLfloat)
SL_WRAP_GL_PROC_2(glTexCoord, 2i, GLint)
SL_WRAP_GL_PROC_2(glTexCoord, 2s, GLshort)

SL_WRAP_GL_PROC_3(glTexCoord, 3d, GLdouble)
SL_WRAP_GL_PROC_3(glTexCoord, 3f, GLfloat)
SL_WRAP_GL_PROC_3(glTexCoord, 3i, GLint)
SL_WRAP_GL_PROC_3(glTexCoord, 3s, GLshort)

SL_WRAP_GL_PROC_4(glTexCoord, 4d, GLdouble)
SL_WRAP_GL_PROC_4(glTexCoord, 4f, GLfloat)
SL_WRAP_GL_PROC_4(glTexCoord, 4i, GLint)
SL_WRAP_GL_PROC_4(glTexCoord, 4s, GLshort)

SL_WRAP_GL_PROC_1(glTexCoord1, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glTexCoord1, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glTexCoord1, iv, const GLint*)
SL_WRAP_GL_PROC_1(glTexCoord1, sv, const GLshort*)

SL_WRAP_GL_PROC_1(glTexCoord2, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glTexCoord2, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glTexCoord2, iv, const GLint*)
SL_WRAP_GL_PROC_1(glTexCoord2, sv, const GLshort*)

SL_WRAP_GL_PROC_1(glTexCoord3, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glTexCoord3, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glTexCoord3, iv, const GLint*)
SL_WRAP_GL_PROC_1(glTexCoord3, sv, const GLshort*)

SL_WRAP_GL_PROC_1(glTexCoord4, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glTexCoord4, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glTexCoord4, iv, const GLint*)
SL_WRAP_GL_PROC_1(glTexCoord4, sv, const GLshort*)

SL_WRAP_GL_PROC_2(glRasterPos, 2d, GLdouble)
SL_WRAP_GL_PROC_2(glRasterPos, 2f, GLfloat)
SL_WRAP_GL_PROC_2(glRasterPos, 2i, GLint)
SL_WRAP_GL_PROC_2(glRasterPos, 2s, GLshort)

SL_WRAP_GL_PROC_3(glRasterPos, 3d, GLdouble)
SL_WRAP_GL_PROC_3(glRasterPos, 3f, GLfloat)
SL_WRAP_GL_PROC_3(glRasterPos, 3i, GLint)
SL_WRAP_GL_PROC_3(glRasterPos, 3s, GLshort)

SL_WRAP_GL_PROC_4(glRasterPos, 4d, GLdouble)
SL_WRAP_GL_PROC_4(glRasterPos, 4f, GLfloat)
SL_WRAP_GL_PROC_4(glRasterPos, 4i, GLint)
SL_WRAP_GL_PROC_4(glRasterPos, 4s, GLshort)

SL_WRAP_GL_PROC_1(glRasterPos2, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glRasterPos2, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glRasterPos2, iv, const GLint*)
SL_WRAP_GL_PROC_1(glRasterPos2, sv, const GLshort*)

SL_WRAP_GL_PROC_1(glRasterPos3, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glRasterPos3, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glRasterPos3, iv, const GLint*)
SL_WRAP_GL_PROC_1(glRasterPos3, sv, const GLshort*)

SL_WRAP_GL_PROC_1(glRasterPos4, dv, const GLdouble*)
SL_WRAP_GL_PROC_1(glRasterPos4, fv, const GLfloat*)
SL_WRAP_GL_PROC_1(glRasterPos4, iv, const GLint*)
SL_WRAP_GL_PROC_1(glRasterPos4, sv, const GLshort*)

SL_WRAP_GL_PROC_4(glRect, d, GLdouble)
SL_WRAP_GL_PROC_4(glRect, f, GLfloat)
SL_WRAP_GL_PROC_4(glRect, i, GLint)
SL_WRAP_GL_PROC_4(glRect, s, GLshort)

SL_WRAP_GL_PROC_2(glRect, dv, const GLdouble*)
SL_WRAP_GL_PROC_2(glRect, fv, const GLfloat*)
SL_WRAP_GL_PROC_2(glRect, iv, const GLint*)
SL_WRAP_GL_PROC_2(glRect, sv, const GLshort*)

// Extension support

extern bool glIsExtensionSupported(const char* extension);

#ifdef _WIN32
#define SL_DECLARE_GLEXT(tp, fn) tp fn
#define SL_GET_GLEXT(tp, fn) fn = wglGetProcAddress(fn)
#else
#define SL_DECLARE_GLEXT(tp, fn)
#define SL_GET_GLEXT(tp, fn)
#endif

#endif



