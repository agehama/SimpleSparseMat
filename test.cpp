#include <ostream>

#define NOMINMAX
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>

#include "SimpleSparseMat.hpp"

#define FOR_MAT(W, H) for (size_t y = 0; y < H; ++y) for (size_t x = 0; x < W; ++x)

namespace ssmat
{
	template<typename T>
	inline std::ostream& operator<<(std::ostream& s, const ssmat::SparseEntry<T>& entry)
	{
		return s << "(" << entry.x << ", " << entry.y << " | " << entry.v << ")";
	}
}

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
		for (size_t y = 0; y < sparseMat.rowCount(); ++y)
		{
			for (size_t i = sparseMat.rowBegin(y); i < sparseMat.rowEnd(y); ++i)
			{
				matOut[y][sparseMat.getX(i)] = sparseMat.getV(i);
			}
		}
	}

	template<typename T, size_t W_A, size_t H_A, size_t W_B>
	void MultipleMat(const T(&matA)[H_A][W_A], const T(&matB)[W_A][W_B], T(&matC)[H_A][W_B])
	{
		for (size_t i = 0; i < H_A; ++i)
		{
			for (size_t j = 0; j < W_B; ++j)
			{
				for (size_t k = 0; k < W_A; ++k)
				{
					matC[i][j] += matA[i][k] * matB[k][j];
				}
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
	MultipleMat(matA, matB, matC);

	const auto sparseMatA = FromArray(matA);
	const auto sparseMatB = FromArray(matB);
	const auto sparseMatC = sparseMatA * sparseMatB;

	int matC_[height][height] = {};
	ToArray(sparseMatC, matC_);

	FOR_MAT(height, height)
	{
		BOOST_CHECK_EQUAL(matC[y][x], matC_[y][x]);
	}
}

BOOST_AUTO_TEST_CASE(SparseMat_Summation)
{
	constexpr int width = 4;
	constexpr int height = 3;
	const int matA[height][width] = {
		{0, 1, 2, 1},
		{2, 3, 0, 5},
		{1, 0, 4, 0},
	};
	const int matB[height][width] = {
		{3, 1, 2, 1},
		{0, 0, 1, 0},
		{1, 0, 3, 4},
	};

	const auto sparseMatA = FromArray(matA);
	const auto sparseMatB = FromArray(matB);
	const auto sparseMatC = sparseMatA + sparseMatB;

	int matC_[height][width] = {};
	ToArray(sparseMatC, matC_);

	FOR_MAT(width, height)
	{
		BOOST_CHECK_EQUAL(matA[y][x] + matB[y][x], matC_[y][x]);
	}
}


BOOST_AUTO_TEST_CASE(SparseMat_Insert)
{
	constexpr int width = 4;
	constexpr int height = 3;
	int matIn[height][width] = {
		{0, 0, 2, 1},
		{2, 0, 0, 5},
		{0, 0, 4, 0},
	};

	auto sparseMatA = FromArray(matIn);

	matIn[1][2] = 3;
	matIn[2][3] = 1;
	sparseMatA.insert(2, 1, 3);
	sparseMatA.insert(3, 2, 1);

	int matOut[height][width] = {};
	ToArray(sparseMatA, matOut);

	FOR_MAT(height, height)
	{
		BOOST_CHECK_EQUAL(matIn[y][x], matOut[y][x]);
	}
}

BOOST_AUTO_TEST_CASE(SparseMat_Append)
{
	constexpr int width = 4;
	constexpr int height = 3;
	int matIn[height][width] = {
		{1, 0, 0, 0},
		{0, 0, 1, 0},
		{0, 1, 0, 0},
	};

	auto sparseMatA = FromArray(matIn);
	matIn[1][3] = 3;
	matIn[0][3] = 1;
	matIn[2][0] = 1;
	std::vector<ssmat::SparseEntry<int>> additionalEntries;
	additionalEntries.emplace_back(3, 1, 3);
	additionalEntries.emplace_back(3, 0, 1);
	additionalEntries.emplace_back(0, 2, 1);
	sparseMatA.append(additionalEntries);

	int matOut[height][width] = {};
	ToArray(sparseMatA, matOut);

	FOR_MAT(height, height)
	{
		BOOST_CHECK_EQUAL(matIn[y][x], matOut[y][x]);
	}
}

BOOST_AUTO_TEST_CASE(SparseMat_Fill0)
{
	ssmat::SparseMat<int> sparseMat;
	sparseMat.fill(3, 7, 0);

	const auto& rowBeginIndices = sparseMat.getRowBeginIndices();
	const auto& xs = sparseMat.getXs();
	const auto& vs = sparseMat.getVs();

	BOOST_CHECK_EQUAL(rowBeginIndices.size(), 0);
	BOOST_CHECK_EQUAL(xs.size(), 0);
	BOOST_CHECK_EQUAL(vs.size(), 0);
}

BOOST_AUTO_TEST_CASE(SparseMat_Fill1)
{
	constexpr int width = 3;
	constexpr int height = 7;

	ssmat::SparseMat<int> sparseMat;
	sparseMat.fill(width, height, 2);
	sparseMat.fill(height, width, 1);

	int matOut[height][width] = {};
	ToArray(sparseMat, matOut);

	FOR_MAT(width, height)
	{
		BOOST_CHECK_EQUAL(matOut[y][x], 1);
	}
}

BOOST_AUTO_TEST_SUITE_END()
