#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "SimpleSparseMat.hpp"


BOOST_AUTO_TEST_SUITE(aa)

BOOST_AUTO_TEST_CASE(init_test1)
{
	constexpr int width = 4;
	constexpr int height = 3;
	const int matIn[height][width] = {
		{0, 1, 2, 1},
		{2, 3, 0, 5},
		{1, 0, 4, 0},
	};

	std::vector<ssmat::SparseEntry<int>> entries;
	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width; ++x)
		{
			if (matIn[y][x] != 0)
			{
				entries.emplace_back(x, y, matIn[y][x]);
			}
		}
	}

	const ssmat::SparseMat<int> sparseMat(entries);
	const auto& columnBeginIndices = sparseMat.getRowBeginIndices();
	const auto& xs = sparseMat.getXs();
	const auto& vs = sparseMat.getVs();

	BOOST_CHECK_EQUAL(vs.size(), xs.size());

	int matOut[height][width] = {};
	for (size_t y = 0; y + 1 < columnBeginIndices.size(); ++y)
	{
		const size_t rowBegin = columnBeginIndices[y];
		const size_t rowEnd = columnBeginIndices[y + 1];
		for (size_t i = rowBegin; i < rowEnd; ++i)
		{
			matOut[y][xs[i]] = vs[i];
		}
	}

	for (size_t y = 0; y < height; ++y)
	{
		for (size_t x = 0; x < width; ++x)
		{
			BOOST_CHECK_EQUAL(matIn[y][x], matOut[y][x]);
		}
	}
}


BOOST_AUTO_TEST_SUITE_END()
