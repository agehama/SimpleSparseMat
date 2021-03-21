#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "SimpleSparseMat.hpp"

#define FOR_MAT(W, H) for (size_t y = 0; y < H; ++y) for (size_t x = 0; x < W; ++x)

namespace
{
	template<typename T, size_t W, size_t H>
	ssmat::SparseMat<T> FromArray(const T(&vs)[H][W])
	{
		std::vector<ssmat::SparseEntry<T>> entries;
		entries.reserve(W * H);
		for (ssmat::IndexT y = 0; y < H; ++y)
		{
			for (ssmat::IndexT x = 0; x < W; ++x)
			{
				const T val = vs[y][x];
				if (val != static_cast<T>(0))
				{
					entries.emplace_back(x, y, val);
				}
			}
		}
		return ssmat::SparseMat<T>(entries);
	}

	template<typename T, size_t W, size_t H>
	void ToArray(const ssmat::SparseMat<T>& sparseMat, T(&matOut)[H][W])
	{
		const auto& columnBeginIndices = sparseMat.getRowBeginIndices();
		const auto& xs = sparseMat.getXs();
		const auto& vs = sparseMat.getVs();

		BOOST_CHECK_EQUAL(vs.size(), xs.size());

		for (size_t y = 0; y + 1 < columnBeginIndices.size(); ++y)
		{
			const size_t rowBegin = columnBeginIndices[y];
			const size_t rowEnd = columnBeginIndices[y + 1];
			for (size_t i = rowBegin; i < rowEnd; ++i)
			{
				matOut[y][xs[i]] = vs[i];
			}
		}
	}
}

BOOST_AUTO_TEST_SUITE(aa)

BOOST_AUTO_TEST_CASE(SparseMat_Init)
{
	constexpr int width = 4;
	constexpr int height = 3;
	const int matIn[height][width] = {
		{0, 1, 2, 1},
		{2, 3, 0, 5},
		{1, 0, 4, 0},
	};

	const auto sparseMat = FromArray(matIn);

	int matOut[height][width] = {};
	ToArray(sparseMat, matOut);

	FOR_MAT(width, height)
	{
		BOOST_CHECK_EQUAL(matIn[y][x], matOut[y][x]);
	}
}

BOOST_AUTO_TEST_CASE(SparseMat_Init2)
{
	auto data = ssmat::MakeEntries(
		{ 2,1,3,3,2,0,0,1 },
		{ 2,1,0,1,0,2,1,0 },
		{ 4,3,1,5,2,1,2,1 }
	);
	const auto entriesIn = SortEntries(data);

	const ssmat::SparseMat<int> sparseMat(entriesIn);
	const auto entriesOut = sparseMat.decompressEntries();

	BOOST_CHECK_EQUAL(entriesIn.size(), entriesOut.size());
	for (size_t i = 0; i < entriesIn.size(); ++i)
	{
		BOOST_CHECK_EQUAL(entriesIn[i], entriesOut[i]);
	}
}

BOOST_AUTO_TEST_CASE(SparseMat_Multiplication)
{
	constexpr int width = 4;
	constexpr int height = 3;
	const int matA[height][width] = {
		{0, 1, 2, 1},
		{2, 3, 0, 5},
		{1, 0, 4, 0},
	};
	const int matB[width][height] = {
		{0, 1, 2},
		{1, 4, 0},
		{0, 0, 1},
		{0, 3, 0},
	};
	int matC[height][height] = {};
	for (size_t i = 0; i < height; ++i)
	{
		for (size_t j = 0; j < height; ++j)
		{
			for (size_t k = 0; k < width; ++k)
			{
				matC[i][j] += matA[i][k] * matB[k][j];
			}
		}
	}

	const auto sparseMatA = FromArray(matA);
	const auto sparseMatB = FromArray(matB);
	const auto sparseMatC = sparseMatA * sparseMatB;

	int matC_[height][height] = {};
	ToArray(sparseMatC, matC_);

	FOR_MAT(height, height)
	{
		BOOST_CHECK_EQUAL(matC[y][x], matC_[y][x]);
	}

	BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()
