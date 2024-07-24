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
            #define JMTX_NODISCARD_FUNCTION __attribute__((warn_unused_result))
    #endif

    #ifndef JMTX_HOT_FUNCTION
        #define JMTX_HOT_FUNCTION __attribute__((hot))
    #endif
#endif



//  Fallback
#ifndef JMTX_NODISCARD_FUNCTION
    #define JMTX_NODISCARD_FUNCTION
#endif

#ifndef JMTX_HOT_FUNCTION
    #define JMTX_HOT_FUNCTION
#endif

#ifndef JMTX_EXTERNAL_FUNCTION
    #define JMTX_EXTERNAL_FUNCTION
#endif

#ifndef JMTX_INTERNAL_FUNCTION
    #define JMTX_INTERNAL_FUNCTION
#endif

#ifndef _MSC_BUILD
    #define JMTX_ARRAY_ATTRIB(x) x
#else//_MSC_BUILD
    #define JMTX_ARRAY_ATTRIB(x)
#endif//_MSC_BUILD

#endif //JMTX_COMMON_H
