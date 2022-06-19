//
// Created by jan on 18.6.2022.
//

#ifndef MTXLIB_ERRORS_H
#define MTXLIB_ERRORS_H

#include <stdint.h>
#include <stdlib.h>
#include <errno.h>


/**
 * Compilation flags:
 *  - MTX_ERROR_MESSAGES: turns on error messages which can be intercepted
 */


#ifdef MTX_ERROR_MESSAGES

struct ErrorCallStack
{
    char* name;
    const char** callers;
    const char** files;
    int* lines;
    const char** callees;
    uint capacity;
    uint size;
};

_Thread_local extern struct ErrorCallStack error_thread_call_stack;
void error_enter_function(struct ErrorCallStack* stack, const char* caller, const char* file, int line, const char* callee);
void error_leave_function(struct ErrorCallStack* stack);
#define CALL_FUNCTION(exp) (error_enter_function(&error_thread_call_stack, __func__, __FILE__, __LINE__, #exp), exp)
#define LEAVE_FUNCTION() error_leave_function(&error_thread_call_stack)

struct ErrorCallStack error_thread_begin(const char* name, ...);
void error_thread_end(struct ErrorCallStack* stack);
#define THREAD_BEGIN(fmt, ...) error_thread_call_stack = error_thread_begin(fmt __VA_OPT__(,) __VA_ARGS__)
#define THREAD_END error_thread_end(&error_thread_call_stack)

void report_error_message(struct ErrorCallStack* stack, const char* fmt, ...);
#define REPORT_ERROR_MESSAGE(msg, ...) report_error_message(&error_thread_call_stack, "Error in function %s in %s - %d:\n\t" msg "\n", __func__, __FILE__, __LINE__ __VA_OPT__(,) __VA_ARGS__)
#define MALLOC_FAILED(size) REPORT_ERROR_MESSAGE("Could not allocate %zu bytes using malloc, reason: %s", (size_t)(size), strerror(errno))
#define REALLOC_FAILED(size) REPORT_ERROR_MESSAGE("Could not allocate %zu bytes using realloc, reason: %s", (size_t)(size), strerror(errno))
#define CALLOC_FAILED(size) REPORT_ERROR_MESSAGE("Could not allocate %zu bytes using calloc, reason: %s", (size_t)(size), strerror(errno))
void set_error_hook(void(*new_hook)(const char* error_message, void* param), void* param);



#else
#define CALL_FUNCTION(exp) exp
#define LEAVE_FUNCTION() (void)0
#define REPORT_ERROR_MESSAGE(msg, ...) (void)0
#define MALLOC_FAILED(size) (void)0
#define REALLOC_FAILED(size) (void)0
#define CALLOC_FAILED(size) (void)0
#define THREAD_BEGIN(fmt, ...) (void)0
#define THREAD_END (void)0
#endif



#endif //MTXLIB_ERRORS_H
