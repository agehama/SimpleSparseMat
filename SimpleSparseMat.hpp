#pragma once
#include <vector>
#include <algorithm>

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

	enum class SparseFormat { CSR, CSC };

	template<typename T>
	std::vector<SparseEntry<T>>& SortEntries(std::vector<SparseEntry<T>>& entries, SparseFormat format)
	{
		switch (format)
		{
		case ssmat::SparseFormat::CSR:
			std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
				return a.y == b.y ? a.x < b.x : a.y < b.y;
				});
			break;
		case ssmat::SparseFormat::CSC:
			std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
				return a.x == b.x ? a.y < b.y : a.x < b.x;
				});
			break;
		default:
			break;
		}

		return entries;
	}

	template<typename T>
	class SparseMat
	{
	public:
		SparseMat() = default;

		SparseMat(const std::vector<SparseEntry<T>>& sortedEntries, SparseFormat format) : format(format)
		{
			if (sortedEntries.empty())
			{
				return;
			}

			switch (format)
			{
			case ssmat::SparseFormat::CSR:
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
				break;
			case ssmat::SparseFormat::CSC:
				for (size_t i = 0; i < sortedEntries.size(); ++i)
				{
					const auto& entry = sortedEntries[i];
					while (rowBeginIndices.size() < entry.x + 1)
					{
						rowBeginIndices.emplace_back(static_cast<IndexT>(xs.size()));
					}

					xs.push_back(entry.y);
					vs.push_back(entry.v);
				}
				break;
			default:
				break;
			}
		}

		void initFrom(const std::vector<SparseEntry<T>>& sortedEntries)
		{
			*this = SparseMat<T>(sortedEntries);
		}

		void fill(IndexT rowCount, IndexT columnCount, T value)
		{
			if (value == 0)
			{
				rowBeginIndices.clear();
				xs.clear();
				vs.clear();
			}
			else
			{
				const IndexT elemCount = rowCount * columnCount;
				rowBeginIndices.resize(rowCount);
				std::iota(rowBeginIndices.begin(), rowBeginIndices.end(), 0);
				std::transform(rowBeginIndices.begin(), rowBeginIndices.end(), rowBeginIndices.begin(), [columnCount](IndexT y) { return y * columnCount; });

				xs.resize(elemCount);
				std::iota(xs.begin(), xs.end(), 0);
				std::transform(xs.begin(), xs.end(), xs.begin(), [columnCount](IndexT i) { return i % columnCount; });

				vs.resize(elemCount);
				std::fill(vs.begin(), vs.end(), value);
			}
		}

		std::vector<SparseEntry<T>> decompressEntries()const
		{
			std::vector<SparseEntry<T>> entries(xs.size());

			switch (format)
			{
			case ssmat::SparseFormat::CSR:
				for (size_t y = 0; y < rowCount(); ++y)
				{
					for (size_t i = rowBegin(y); i < rowEnd(y); ++i)
					{
						auto& entry = entries[i];
						entry.y = y;
						entry.x = xs[i];
						entry.v = vs[i];
					}
				}
				break;
			case ssmat::SparseFormat::CSC:
				for (size_t x = 0; x < rowCount(); ++x)
				{
					for (size_t i = rowBegin(x); i < rowEnd(x); ++i)
					{
						auto& entry = entries[i];
						entry.y = xs[i];
						entry.x = x;
						entry.v = vs[i];
					}
				}
				break;
			default:
				break;
			}

			return entries;
		}

		SparseMat<T> operator*(const SparseMat<T>& B)const
		{
			std::vector<SparseEntry<T>> results;
			for (IndexT yA = 0; yA < rowCount(); ++yA)
			{
				std::map<IndexT, T> currentRow;
				for (IndexT yB = 0; yB < B.rowCount(); ++yB)
				{
					if (B.rowBegin(yB) == B.rowEnd(yB))
					{
						continue;
					}

					for (IndexT iA = rowBegin(yA); iA < rowEnd(yA); ++iA)
					{
						const IndexT xA = xs[iA];
						if (xA != yB)
						{
							continue;
						}

						for (IndexT iB = B.rowBegin(yB); iB < B.rowEnd(yB); ++iB)
						{
							const IndexT xB = B.xs[iB];
							currentRow[xB] += vs[iA] * B.vs[iB];
						}
					}
				}

				for (const auto& entry : currentRow)
				{
					results.emplace_back(entry.first, yA, entry.second);
				}
			}

			return SparseMat<T>(results, SparseFormat::CSR);
		}

		SparseMat<T> operator+(const SparseMat<T>& B)const
		{
			std::vector<SparseEntry<T>> results;

			const IndexT maxRowCount = std::max(rowCount(), B.rowCount());
			for (IndexT y = 0; y < maxRowCount; ++y)
			{
				std::map<IndexT, T> currentRow;

				if (y < rowBeginIndices.size())
				{
					for (IndexT i = rowBegin(y); i < rowEnd(y); ++i)
					{
						currentRow[xs[i]] = vs[i];
					}
				}

				if (y < B.rowBeginIndices.size())
				{
					for (IndexT i = B.rowBegin(y); i < B.rowEnd(y); ++i)
					{
						currentRow[B.xs[i]] += B.vs[i];
					}
				}

				for (const auto& entry : currentRow)
				{
					results.emplace_back(entry.first, y, entry.second);
				}
			}

			return SparseMat<T>(results, SparseFormat::CSR);
		}

		void transpose()
		{
			auto cooForm = decompressEntries();
			for (auto& entry : cooForm)
			{
				const IndexT temp = entry.x;
				entry.x = entry.y;
				entry.y = temp;
			}
			*this = SparseMat(SortEntries(cooForm, format), format);
		}

		void insert(IndexT x, IndexT y, T v)
		{
			auto cooForm = decompressEntries();
			cooForm.emplace_back(x, y, v);
			*this = SparseMat(SortEntries(cooForm, format), format);
		}

		void append(const std::vector<SparseEntry<T>>& entries)
		{
			auto cooForm = decompressEntries();
			cooForm.insert(cooForm.end(), entries.begin(), entries.end());
			*this = SparseMat(SortEntries(cooForm, format), format);
		}

		IndexT rowCount()const { return rowBeginIndices.size(); }
		IndexT rowBegin(IndexT row)const { return rowBeginIndices[row]; }
		IndexT rowEnd(IndexT row)const { return rowBeginIndices.size() <= row + 1 ? xs.size() : rowBeginIndices[row + 1]; }
		IndexT getX(IndexT i)const { return xs[i]; }
		T getV(IndexT i)const { return vs[i]; }

		const std::vector<IndexT>& getRowBeginIndices()const { return rowBeginIndices; }
		const std::vector<IndexT>& getXs()const { return xs; }
		const std::vector<T>& getVs()const { return vs; }
		std::vector<T>& getVs() { return vs; }
		SparseFormat getFormat()const { return format; }

	private:
		std::vector<IndexT> rowBeginIndices;
		std::vector<IndexT> xs;
		std::vector<T> vs;
		SparseFormat format = SparseFormat::CSR;
	};
}
