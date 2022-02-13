#ifndef TL_ASSERT_H
#define TL_ASSERT_H

#include <iostream>

#ifndef NDEBUG
#define TL_ASSERT(Expr, Msg) __M_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#define TL_ASSERT(Expr, Msg) ;
#endif

template <typename T>
void __M_Assert(const char* expr_str, bool expr, const char* file, int line,
                const T& msg) {
    if (!expr) {
        std::cerr << "Assert failed:\t" << msg << "\n"
                  << "Expected:\t" << expr_str << "\n"
                  << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}

#endif  // TL_ASSERT_H
