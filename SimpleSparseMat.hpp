#pragma once
#include <vector>

namespace ssmat
{
	using IndexT = std::int32_t;

	template<typename T>
	struct SparseEntry
	{
		SparseEntry(IndexT x, IndexT y, T v) :x(x), y(y), v(v) {}
		IndexT x, y;
		T v;
	};

	template<typename T>
	class SparseMat
	{
	public:
		SparseMat(const std::vector<SparseEntry<T>>& sortedEntries)
		{
			if (sortedEntries.empty())
			{
				return;
			}

			for (size_t i = 0; i < sortedEntries.size(); ++i)
			{
				const auto& entry = sortedEntries[i];
				while (rowBeginIndices.size() < entry.y + 1)
				{
					rowBeginIndices.emplace_back(static_cast<IndexT>(xs.size()));
				}

				xs.push_back(entry.x);
				vs.push_back(entry.v);
			}

			rowBeginIndices.emplace_back(static_cast<IndexT>(xs.size()));
		}

		const std::vector<IndexT>& getRowBeginIndices()const { return rowBeginIndices; }
		const std::vector<IndexT>& getXs()const { return xs; }
		const std::vector<T>& getVs()const { return vs; }

	private:
		std::vector<IndexT> rowBeginIndices;
		std::vector<IndexT> xs;
		std::vector<T> vs;
	};
}
