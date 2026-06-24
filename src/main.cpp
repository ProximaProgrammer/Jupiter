#include <iostream>

#include "core/Grid.hpp"

int main()
{
    Jupiter::Grid grid(
        64,
        64,
        128);

    std::cout
        << "Cells = "
        << static_cast<long long>(
            64LL * 64LL * 128LL)
        << '\n';

    return 0;
}