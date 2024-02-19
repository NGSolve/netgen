#ifndef INCOPENGL_HPP___
#define INCOPENGL_HPP___
#define GL_GLEXT_PROTOTYPES

#include <mystdlib.h>
#include <mydefs.hpp>

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#  ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#define GL_DO_NOT_WARN_IF_MULTI_GL_VERSION_HEADERS_INCLUDED
#    include <OpenGL/gl3.h>
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
#define GL_FRAMEBUFFER_COMPLETE           0x8CD5
#define GL_FRAMEBUFFER                    0x8D40
#define GL_RENDERBUFFER                   0x8D41
#define GL_DEPTH_ATTACHMENT               0x8D00
#define GL_COLOR_ATTACHMENT0              0x8CE0
typedef ptrdiff_t GLintptr;
typedef ptrdiff_t GLsizeiptr;
extern void (*glBindBuffer) (GLenum a, GLuint b);
extern void (*glDeleteBuffers) (GLsizei a, const GLuint *b);
extern void (*glGenBuffers) (GLsizei a, GLuint *b);
extern void (*glBufferData) (GLenum a, GLsizeiptr b, const GLvoid *c, GLenum d);
extern void (*glBufferSubData) (GLenum a, GLintptr b, GLsizeiptr c, const GLvoid *d);

extern GLenum (*glCheckFramebufferStatus) (GLenum target);
extern void (*glBindFramebuffer) (GLenum target, GLuint framebuffer);
extern void (*glBindRenderbuffer) (GLenum target, GLuint renderbuffer);
extern void (*glDeleteFramebuffers) (GLsizei n, const GLuint *framebuffers);
extern void (*glDeleteRenderbuffers) (GLsizei n, const GLuint *renderbuffers);
extern void (*glGenFramebuffers) (GLsizei n, GLuint *framebuffers);
extern void (*glGenRenderbuffers) (GLsizei n, GLuint *renderbuffers);
extern void (*glRenderbufferStorage) (GLenum target, GLenum internalformat, GLsizei width, GLsizei height);
extern void (*glFramebufferRenderbuffer) (GLenum target, GLenum attachment, GLenum renderbuffertarget, GLuint renderbuffer);
#endif
DLL_HEADER void LoadOpenGLFunctionPointers();
#endif // INCOPENGL_HPP___
