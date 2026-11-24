/*
 * diffuse_core.h - C interface to the DISCUS suite (libdiffuse_core)
 *
 * Implemented in python/diffuse_capi.f90 (Fortran, bind(C)). Used by the
 * Python package through ctypes; usable from C, C++, Julia and others.
 * See migration/WP2_Binding_Decision.md (option B).
 *
 * Conventions
 *   - strings in : char array + length; no terminating zero needed
 *   - strings out: buffers of the sizes given below, blank padded, no zero
 *   - arrays     : pointer + size; the caller allocates
 *   - sizes      : int64_t; flags and codes: int (0 = ok, -1 = error)
 *   - one suite per process; calls must not overlap (no thread safety)
 *   - output goes to the process's stdout/stderr (capture it if needed)
 *   - library mode: an error never ends the process; the first error since
 *     dc_clear_error is kept for dc_last_error
 */
#ifndef DIFFUSE_CORE_H
#define DIFFUSE_CORE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DC_API_VERSION      1
#define DC_ERR_CLASS_LEN    4     /* "COMM", "I/O ", "FORT" ...            */
#define DC_ERR_TEXT_LEN     45    /* standard message                      */
#define DC_ERR_MSGS_LEN     560   /* 7 extra lines of 80 characters        */
#define DC_MAX_LINE         1024  /* longest command line                  */

int     dc_api_version(void);                 /* DC_API_VERSION of the library */

void    dc_init(void);                        /* once, before anything else    */
void    dc_reset(void);                       /* the suite command "reset"     */
void    dc_flush_output(void);                /* flush Fortran's stdout/stderr */

/* Run command lines separated by '\n' (a macro held in memory). Returns 0
 * when they ran (check dc_last_error for an error), else the length of the
 * longest line, which exceeds DC_MAX_LINE (nothing was run).               */
int64_t dc_run(const char *text, int64_t n);

/* Value of a DISCUS expression; 0 on error (see dc_last_error). About 12
 * significant digits: as many as DISCUS's expression evaluator keeps.       */
double  dc_value(const char *expr, int64_t n);
/* Set a user variable (declared with "variable real, name"). Built-in
 * arrays such as r[] are set with a command: dc_run("r[2] = 2.5").          */
void    dc_set_value(const char *var, int64_t n, double value);

void    dc_clear_error(void);
void    dc_last_error(int *num, int *type,
                      char cls[DC_ERR_CLASS_LEN], char text[DC_ERR_TEXT_LEN],
                      char msgs[DC_ERR_MSGS_LEN]);

void    dc_seed(int value);                   /* 0: seed from the clock        */
int64_t dc_random_size(void);
int     dc_get_random_state(int *state, int64_t n);   /* n = dc_random_size() */
int     dc_set_random_state(const int *state, int64_t n);

int     dc_data_sets(void);                   /* KUPLOT data sets loaded       */
int64_t dc_data_length(int ik);               /* points of set ik (1-based), -1 */
int     dc_get_data(int ik, float *x, float *y, int64_t n);
int     dc_get_data_2d(int ik, float *x, float *y, float *z, int64_t n,
                       int64_t *nx, int64_t *ny);

#ifdef __cplusplus
}
#endif

#endif /* DIFFUSE_CORE_H */
