#ifndef _COMMON_API_H_
#define _COMMON_API_H_

#ifdef _MSC_VER
#    define COMMON_API __declspec(dllexport)
#    define _ITERATOR_DEBUG_LEVEL 0
#else
#  define COMMON_API
#endif

#endif // API_H

