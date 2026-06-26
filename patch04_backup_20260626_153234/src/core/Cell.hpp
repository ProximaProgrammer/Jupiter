#pragma once

#include "CellGeometry.hpp"
#include "CellIndex.hpp"
#include "../state/ConservedState.hpp"
#include "../state/PrimitiveState.hpp"

namespace Jupiter
{

class Cell
{
public:
    CellIndex index;
    CellGeometry geometry;
    ConservedState conserved;
    PrimitiveState primitive;
};

} // namespace Jupiter
