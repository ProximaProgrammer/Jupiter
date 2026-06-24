#pragma once

#include "../state/ConservedState.hpp"
#include "../state/PrimitiveState.hpp"

#include "CellGeometry.hpp"
#include "CellIndex.hpp"

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

}