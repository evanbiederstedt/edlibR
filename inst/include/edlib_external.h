#ifndef EDLIB_EXTERNAL_H
#define EDLIB_EXTERNAL_H

#include "edlib.h"
#include <R_ext/Rdynload.h>

// use snake_case for function names to differentiate exported functions
// Camel case is used for the internal source file functions

inline EdlibAlignConfig edlib_new_align_config(int k, EdlibAlignMode mode, EdlibAlignTask task,
                                               const EdlibEqualityPair* additional_equalities,
                                               int additional_equalities_length) {
    using fun_t = EdlibAlignConfig(*)(int, EdlibAlignMode, EdlibAlignTask,
                                      const EdlibEqualityPair*, int);
    static fun_t fun = (fun_t)R_GetCCallable("edlibR","edlibNewAlignConfig");
    return fun(k, mode, task, additional_equalities, additional_equalities_length);
}

inline EdlibAlignConfig edlib_default_align_config(void) {
    using fun_t = EdlibAlignConfig(*)(void);
    static fun_t fun = (fun_t)R_GetCCallable("edlibR","edlibDefaultAlignConfig");
    return fun();
}

inline void edlib_free_align_result(EdlibAlignResult result) {
    using fun_t = void(*)(EdlibAlignResult);
    static fun_t fun = (fun_t)R_GetCCallable("edlibR","edlibFreeAlignResult");
    fun(result);
}

inline EdlibAlignResult edlib_align(const char* query, int query_length,
                                    const char* target, int target_length,
                                    EdlibAlignConfig config) {
    using fun_t = EdlibAlignResult(*)(const char*,int,const char*,int,EdlibAlignConfig);
    static fun_t fun = (fun_t)R_GetCCallable("edlibR","edlibAlign");
    return fun(query, query_length, target, target_length, config);
}

inline char* edlib_alignment_to_cigar(const unsigned char* alignment, int alignment_length,
                                      EdlibCigarFormat cigar_format) {
    using fun_t = char*(*)(const unsigned char*,int,EdlibCigarFormat);
    static fun_t fun = (fun_t)R_GetCCallable("edlibR","edlibAlignmentToCigar");
    return fun(alignment, alignment_length, cigar_format);
}

#endif // EDLIB_EXTERNAL_H
