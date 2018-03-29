#include "Vector.h"
#include "VectorUtils.h"

#include <iostream>

int main()
{
  using namespace VectorUtils;
  std::cout << sum(Vector(1, 2, 3)) << std::endl;
}
