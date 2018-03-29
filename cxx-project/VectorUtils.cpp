#include "VectorUtils.h"
#include "Vector.h"


double
VectorUtils::
sum(Vector const &vector)
{
  return vector.data[0] + vector.data[1] + vector.data[2];
}
