#pragma once
#include <vector>

namespace ssmat
{
	using IndexT = std::int32_t;

	template<typename T>
	struct SparseEntry
	{
		SparseEntry() = default;
		SparseEntry(IndexT x, IndexT y, T v) :x(x), y(y), v(v) {}
		IndexT x, y;
		T v;
		bool operator==(const SparseEntry& a)const { return x == a.x && y == a.y && v == a.v; }
	};

	template<typename T>
	inline std::ostream& operator<<(std::ostream& s, const SparseEntry<T>& entry)
	{
		return s << "(" << entry.x << ", " << entry.y << " | " << entry.v << ")";
	}

	template<typename T, size_t N>
	std::vector<SparseEntry<T>> MakeEntries(const IndexT(&xs)[N], const IndexT(&ys)[N], const T(&vs)[N])
	{
		std::vector<SparseEntry<T>> entries;
		entries.reserve(N);
		for (size_t i = 0; i < N; ++i)
		{
			entries.emplace_back(xs[i], ys[i], vs[i]);
		}
		return entries;
	}

	template<typename T>
	std::vector<SparseEntry<T>>& SortEntries(std::vector<SparseEntry<T>>& entries)
	{
		std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
			return a.y == b.y ? a.x < b.x : a.y < b.y;
			});
		return entries;
	}

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

		std::vector<SparseEntry<T>> decompressEntries()const
		{
			std::vector<SparseEntry<T>> entries(xs.size());
			for (size_t y = 0; y + 1 < rowBeginIndices.size(); ++y)
			{
				const size_t rowBegin = rowBeginIndices[y];
				const size_t rowEnd = rowBeginIndices[y + 1];
				for (size_t i = rowBegin; i < rowEnd; ++i)
				{
					auto& entry = entries[i];
					entry.y = y;
					entry.x = xs[i];
					entry.v = vs[i];
				}
			}

			return entries;
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
