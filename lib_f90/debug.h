
/* 
       DEBUG.H 
       Fortran debugging output tools 
*/

#define TRACE_MACROS 1
#define DEBUG 0

#if TRACE_MACROS==1

#define MACRO(token) print *, "MACRO: ", token

#else

#define MACRO(token)

#endif

#if DEBUG==1

/** Show a message */
#define MSG(message)  print *, "F90: ", message
/** Show a variable */
#define VAR(variable) print *, "F90: ", "variable: ", variable
/** Show a vector */
#define VEC(variable) print *, "F90: ", "variable:", achar(10), \
                               variable, achar(10), "  --"

#else

#define MSG(message)
#define VAR(variable)
#define VEC(variable)

#endif
