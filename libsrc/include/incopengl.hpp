#ifndef INCOPENGL_HPP___
#define INCOPENGL_HPP___
#define GL_GLEXT_PROTOTYPES

#include <mydefs.hpp>

#  if defined(TOGL_AGL) || defined(TOGL_AGL_CLASSIC) || defined(TOGL_NSOPENGL)
#    include <OpenGL/gl.h>
#    include <OpenGL/glu.h>
#  else
#    include <GL/gl.h>
#    include <GL/glu.h>
#  endif


#ifdef TOGL_X11
// parallel
#define GLX_GLXEXT_PROTOTYPES
#include <GL/glx.h>
#include <GL/glxext.h>
#endif

#ifdef WIN32
// part of OpenGL 1.2, but not in Microsoft's OpenGL 1.1 header:
// GL version should be checked at runtime
#define GL_CLAMP_TO_EDGE                  0x812F
#define GL_ARRAY_BUFFER                   0x8892
#define GL_ELEMENT_ARRAY_BUFFER           0x8893
#define GL_STATIC_DRAW                    0x88E4
typedef ptrdiff_t GLintptr;
typedef ptrdiff_t GLsizeiptr;
extern void (*glBindBuffer) (GLenum a, GLuint b);
extern void (*glDeleteBuffers) (GLsizei a, const GLuint *b);
extern void (*glGenBuffers) (GLsizei a, GLuint *b);
extern void (*glBufferData) (GLenum a, GLsizeiptr b, const GLvoid *c, GLenum d);
extern void (*glBufferSubData) (GLenum a, GLintptr b, GLsizeiptr c, const GLvoid *d);
#endif
DLL_HEADER void LoadOpenGLFunctionPointers();
#endif // INCOPENGL_HPP___
