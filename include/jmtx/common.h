//
// Created by jan on 21.11.2023.
//

#ifndef JMTX_COMMON_H
#define JMTX_COMMON_H

#include <stdlib.h>
#include <stdint.h>

//  Use standard C if possible
#if __STDC_VERSION__ > 201710L
//    #if __has_c_attribute(nodiscard)
    #define JMTX_NODISCARD_FUNCTION [[nodiscard]]
//    #endif
#endif


//  Compiler specific
#ifdef __GNUC__
    #ifndef JMTX_NODISCARD_FUNCTION
            #define JMTX_NODISCARD_FUNCTION __attribute__(nodiscard)
    #endif

    #define JMTX_EXTERNAL_FUNCTION __attribute__((visibility("default")))
    #define JMTX_INTERNAL_FUNCTION __attribute__((visibility("hidden")))
#endif



//  Fallback
#ifndef JMTX_NODISCARD_FUNCTION
    #define JMTX_NODISCARD_FUNCTION
#endif


#ifndef JMTX_EXTERNAL_FUNCTION
    #define JMTX_EXTERNAL_FUNCTION
#endif

#ifndef JMTX_INTERNAL_FUNCTION
    #define JMTX_INTERNAL_FUNCTION
#endif

#endif //JMTX_COMMON_H
