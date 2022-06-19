//
// Created by jan on 18.6.2022.
//

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "errors.h"

_Thread_local struct ErrorCallStack error_thread_call_stack;

void error_enter_function(struct ErrorCallStack* stack, const char* caller, const char* file, int line, const char* callee)
{
    if (stack->capacity == stack->size)
    {
        const uint new_capacity = stack->capacity + 64;

        const char** new_caller_ptr = realloc(stack->callers, sizeof(char*) * new_capacity);
        if (!new_caller_ptr)
        {
            perror("Could not reallocate memory for caller array");
            exit(EXIT_FAILURE);
        }
        memset(new_caller_ptr + stack->size, 0, 64 * sizeof(char*));
        stack->callers = new_caller_ptr;

        const char** new_callee_ptr = realloc(stack->callees, sizeof(char*) * new_capacity);
        if (!new_callee_ptr)
        {
            perror("Could not reallocate memory for callee array");
            exit(EXIT_FAILURE);
        }
        memset(new_callee_ptr + stack->size, 0, 64 * sizeof(char*));
        stack->callees = new_callee_ptr;

        const char** new_files_ptr = realloc(stack->files, sizeof(char*) * new_capacity);
        if (!new_files_ptr)
        {
            perror("Could not reallocate memory for files array");
            exit(EXIT_FAILURE);
        }
        memset(new_files_ptr + stack->size, 0, 64 * sizeof(char*));
        stack->files = new_files_ptr;

        int* new_lines_ptr = realloc(stack->lines, sizeof(int) * new_capacity);
        if (!new_lines_ptr)
        {
            perror("Could not reallocate memory for lines array");
            exit(EXIT_FAILURE);
        }
        memset(new_lines_ptr + stack->size, 0, 64 * sizeof(int));
        stack->lines = new_lines_ptr;
        stack->capacity = new_capacity;
    }

    const uint index = stack->size++;
    stack->lines[index] = line;
    stack->callers[index] = caller;
    stack->files[index] = file;
    stack->callees[index] = callee;
}

void error_leave_function(struct ErrorCallStack* stack)
{
    const uint index = stack->size--;
    assert(index != 0);
}

struct ErrorCallStack error_thread_begin(const char* name, ...)
{
    struct ErrorCallStack result;
    va_list args, args_copy;
    va_start(args, name);
    va_copy(args_copy, args);
    const uint n_bytes = vsnprintf(NULL, 0, name, args) + 1;
    va_end(args);
    char* name_buffer = calloc(n_bytes, sizeof*name_buffer);
    if (!name_buffer)
    {
        perror("Could not allocate memory for thread name");
        exit(EXIT_FAILURE);
    }
    const uint written = vsprintf(name_buffer, name, args_copy);
    assert(written + 1 == n_bytes);
    va_end(args_copy);
    memset(&result, 0, sizeof(result));
    result.name = name_buffer;
    return result;
}

void error_thread_end(struct ErrorCallStack* stack)
{
    free(stack->callees);
    free(stack->callers);
    free(stack->files);
    free(stack->lines);
    free(stack->name);
    memset(stack, 0, sizeof(*stack));
}

static void DEFAULT_ERROR_HOOK_FUNCTION(const char* error_message, void* ignored)
{
    fprintf(stderr, "%s", error_message);
#ifndef NDEBUG
    assert(0);
#endif
    exit(EXIT_FAILURE);
}

static void (*error_hook_fn)(const char* error_message, void* param) = DEFAULT_ERROR_HOOK_FUNCTION;
static void* error_hook_param;

void report_error_message(struct ErrorCallStack* stack, const char* fmt, ...)
{
    assert(error_hook_fn);
    uint bytes_needed = 0;
    va_list args;
    if (!stack->name)
    {
        bytes_needed += snprintf(NULL, 0, "Error occurred:\n");
    }
    else
    {
        bytes_needed += snprintf(NULL, 0, "Error occurred in thread \"%s\":\n", stack->name);
    }
    va_start(args, fmt);
    bytes_needed += vsnprintf(NULL, 0, fmt, args);
    va_end(args);
    bytes_needed += snprintf(NULL, 0, "Call stack information:\n");
    for (uint i = stack->size; i != 0; --i)
    {
        bytes_needed += snprintf(NULL, 0, "function call %s in file %s - %d\n", stack->callees[i - 1], stack->files[i - 1], stack->lines[i - 1]);
    }
    char* buffer = calloc(bytes_needed, sizeof*buffer);
    if (!buffer)
    {
        perror("Could not allocate memory for error message");
        exit(EXIT_FAILURE);
    }

    uint written = 0;
    if (!stack->name)
    {
        written += sprintf(buffer + written, "Error occurred:\n");
    }
    else
    {
        written += sprintf(buffer + written, "Error occurred in thread \"%s\":\n", stack->name);
    }
    va_start(args, fmt);
    written += vsprintf(buffer + written, fmt, args);
    va_end(args);
    written += sprintf(buffer + written, "Call stack information:\n");
    for (uint i = stack->size; i != 0; --i)
    {
        written += sprintf(buffer + written, "function call %s in file %s - %d\n", stack->callees[i - 1], stack->files[i - 1], stack->lines[i - 1]);
    }
    assert(written == bytes_needed);
    error_hook_fn(buffer, error_hook_param);
    free(buffer);
}

void set_error_hook(void (* new_hook)(const char*, void*), void* param)
{
    error_hook_fn = new_hook;
    error_hook_param = param;
}
