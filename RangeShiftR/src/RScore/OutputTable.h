/*----------------------------------------------------------------------------
 *
 *	Copyright (C) 2026 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Roslyn Henry, Théo Pannetier, Jette Wolff, Damaris Zurell
 *
 *	This file is part of RangeShifter.
 *
 *	RangeShifter is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	RangeShifter is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with RangeShifter. If not, see <https://www.gnu.org/licenses/>.
 *
 --------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------

OutputTable: one tabular output file (e.g. Pop, Inds, geneValues) that can be
written either as tab-separated text (the classic RangeShifter format) or as a
binary HDF5 file.

Usage:
	OutputTable t;
	t.addColumn("Rep", OutputTable::INT);
	t.addColumn("Epsilon", OutputTable::FLOAT);
	t.open(pathWithoutExtension, format);   // appends ".txt" or ".h5"
	t << rep << eps;  t.endRow();             // values in column order
	t.close();

Text format: exactly as before, i.e. the header line holds the column names
separated by tabs and every value is written with the default ostream
formatting of the C++ type passed to operator<< (so passing the same types as
the former `ofstream <<` code yields byte-identical files).

HDF5 format (only if compiled with RS_HDF5): one 1-D, chunked, compressed and
extendible dataset per column (STRING columns: variable-length strings) at the root of the file, named like the column.
The root attribute "columns" holds the column names in order, separated by
tabs, so readers can rebuild the table (see readRSOutput() in RangeShiftR).
Values are buffered in memory and written in chunks of `chunkRows` rows.

------------------------------------------------------------------------------*/

#ifndef OutputTableH
#define OutputTableH

#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <cstdint>

enum class OutputFormat : int { TEXT = 0, HDF5 = 1 };

// true if this build can write the given format
bool outputFormatAvailable(OutputFormat fmt);

class OutputTable {
public:
	enum ColType : std::uint8_t { INT, FLOAT, DOUBLE, STRING };

	OutputTable();
	~OutputTable();
	OutputTable(const OutputTable&) = delete;
	OutputTable& operator=(const OutputTable&) = delete;

	// --- schema (call before open) ---
	void clearColumns();
	void addColumn(const std::string& name, ColType type);
	std::size_t nColumns() const { return columns.size(); }

	// --- file handling ---
	// basePath: full path without extension; ".txt" or ".h5" is appended.
	// Returns false if the file could not be created (or HDF5 is unavailable).
	bool open(const std::string& basePath, OutputFormat fmt);
	bool isOpen() const { return opened; }
	void close();
	OutputFormat format() const { return fmt; }

	// --- writing values of the current row, in column order ---
	OutputTable& operator<<(int v);
	OutputTable& operator<<(short v) { return *this << static_cast<int>(v); }
	OutputTable& operator<<(unsigned int v);
	OutputTable& operator<<(long v);
	OutputTable& operator<<(float v);
	OutputTable& operator<<(double v);
	OutputTable& operator<<(const std::string& v); // for STRING columns
	void endRow();

private:
	struct Column {
		std::string name;
		ColType type;
	};
	std::vector<Column> columns;
	OutputFormat fmt = OutputFormat::TEXT;
	bool opened = false;
	std::size_t colIx = 0;     // index of the next value within the current row

	// text backend
	std::ofstream ofs;
	template <typename T> void writeText(T v);

	// binary backend (opaque, defined in OutputTable.cpp)
	struct H5Backend;
	std::unique_ptr<H5Backend> h5;
	template <typename T> void writeBinary(T v);
	void writeBinaryString(const std::string& v);
	void checkColumn();
};

#endif // OutputTableH
