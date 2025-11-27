#ifndef PAR_CUH
#define PAR_CUH

#include <stdbool.h>

#ifdef __cplusplus  // 如果是C++编译器，添加extern "C"
extern "C" {
#endif

// 函数声明
void parse_command_line(int argc, char **argv);
int getparint(const char *key, int *val);
int getparfloat(const char *key, float *val);
int getparstring(const char *key, char *val);
int getparbool(const char *key, bool *val);
void err(const char *msg);

#ifdef __cplusplus
}
#endif

#endif // PAR_CUH