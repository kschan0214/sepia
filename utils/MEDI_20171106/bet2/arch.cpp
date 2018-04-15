#include <math.h>
#include "arch.h"
#if _MSC_VER < 1800
float roundf(float x)
{
  return x >= 0.0f ? floorf(x + 0.5f) : ceilf(x - 0.5f);
}
#endif
