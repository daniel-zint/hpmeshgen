#pragma once

namespace Contour
{
	using ContourTypeT = int;
	namespace ContourType
	{
		constexpr ContourTypeT INSIDE = 0;

		constexpr ContourTypeT LAND = -100;
		constexpr ContourTypeT ISLAND = -101;
		constexpr ContourTypeT OPEN_SEA = -1100;
		constexpr ContourTypeT RIVER = -1101;
		constexpr ContourTypeT RADIATION = -1102;

		constexpr ContourTypeT NOTHING = -5000;

	}
}