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

#include "OutputTable.h"

#include <stdexcept>
#include <cstring>

#ifdef RS_HDF5
#include <hdf5.h>
#endif

bool outputFormatAvailable(OutputFormat fmt) {
	if (fmt == OutputFormat::TEXT) return true;
#ifdef RS_HDF5
	return fmt == OutputFormat::HDF5;
#else
	return false;
#endif
}

//---------------------------------------------------------------------------
// HDF5 backend
//---------------------------------------------------------------------------
#ifdef RS_HDF5
struct OutputTable::H5Backend {
	struct Col {
		hid_t dset = -1;
		hid_t memtype = -1;     // native type used in memory
		std::size_t esize = 0;  // bytes per value
		std::vector<unsigned char> buf;
		bool isString = false;
		std::vector<std::string> sbuf; // STRING columns
	};
	hid_t file = -1;
	std::vector<Col> cols;
	hsize_t nWritten = 0;       // rows already in the file
	std::size_t nBuffered = 0;  // complete rows in the buffers
	static constexpr std::size_t chunkRows = 32768;

	bool create(const std::string& path, const std::vector<OutputTable::Column>& columns) {
		H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); // no HDF5 error stack printing
		file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (file < 0) return false;

		hsize_t dims[1] = { 0 }, maxdims[1] = { H5S_UNLIMITED }, chunk[1] = { chunkRows };
		hid_t space = H5Screate_simple(1, dims, maxdims);
		hid_t dcpl = H5Pcreate(H5P_DATASET_CREATE);
		H5Pset_chunk(dcpl, 1, chunk);
		if (H5Zfilter_avail(H5Z_FILTER_SHUFFLE) > 0) H5Pset_shuffle(dcpl);
		if (H5Zfilter_avail(H5Z_FILTER_DEFLATE) > 0) H5Pset_deflate(dcpl, 4);

		bool ok = true;
		std::string names;
		for (const auto& c : columns) {
			Col col;
			hid_t filetype;
			switch (c.type) {
			case OutputTable::INT:    col.memtype = H5T_NATIVE_INT32;  filetype = H5T_STD_I32LE;  col.esize = 4; break;
			case OutputTable::FLOAT:  col.memtype = H5T_NATIVE_FLOAT;  filetype = H5T_IEEE_F32LE; col.esize = 4; break;
			case OutputTable::STRING:
				col.memtype = H5Tcopy(H5T_C_S1);
				H5Tset_size(col.memtype, H5T_VARIABLE);
				filetype = col.memtype;
				col.isString = true;
				break;
			default:                  col.memtype = H5T_NATIVE_DOUBLE; filetype = H5T_IEEE_F64LE; col.esize = 8; break;
			}
			col.dset = H5Dcreate2(file, c.name.c_str(), filetype, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
			if (col.dset < 0) ok = false;
			col.buf.reserve(chunkRows * col.esize);
			cols.push_back(std::move(col));
			if (!names.empty()) names += '\t';
			names += c.name;
		}
		H5Pclose(dcpl);
		H5Sclose(space);

		// column order as root attribute
		hid_t strtype = H5Tcopy(H5T_C_S1);
		H5Tset_size(strtype, names.size() > 0 ? names.size() : 1);
		hid_t aspace = H5Screate(H5S_SCALAR);
		hid_t attr = H5Acreate2(file, "columns", strtype, aspace, H5P_DEFAULT, H5P_DEFAULT);
		if (attr < 0 || H5Awrite(attr, strtype, names.c_str()) < 0) ok = false;
		if (attr >= 0) H5Aclose(attr);
		H5Sclose(aspace);
		H5Tclose(strtype);
		return ok;
	}

	template <typename T>
	void append(std::size_t ix, OutputTable::ColType type, T v) {
		unsigned char bytes[8];
		switch (type) {
		case OutputTable::INT:   { std::int32_t x = static_cast<std::int32_t>(v); std::memcpy(bytes, &x, 4); break; }
		case OutputTable::FLOAT: { float x = static_cast<float>(v); std::memcpy(bytes, &x, 4); break; }
		default:                 { double x = static_cast<double>(v); std::memcpy(bytes, &x, 8); break; }
		}
		auto& b = cols[ix].buf;
		b.insert(b.end(), bytes, bytes + cols[ix].esize);
	}

	void appendString(std::size_t ix, const std::string& v) {
		cols[ix].sbuf.push_back(v);
	}

	void rowCompleted() {
		if (++nBuffered >= chunkRows) flush();
	}

	void flush() {
		if (nBuffered == 0) return;
		hsize_t newdims[1] = { nWritten + nBuffered };
		hsize_t start[1] = { nWritten }, count[1] = { nBuffered };
		hid_t mspace = H5Screate_simple(1, count, nullptr);
		for (auto& c : cols) {
			H5Dset_extent(c.dset, newdims);
			hid_t fspace = H5Dget_space(c.dset);
			H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, nullptr, count, nullptr);
			herr_t status;
			if (c.isString) {
				std::vector<const char*> ptrs;
				ptrs.reserve(c.sbuf.size());
				for (const auto& s : c.sbuf) ptrs.push_back(s.c_str());
				status = H5Dwrite(c.dset, c.memtype, mspace, fspace, H5P_DEFAULT, ptrs.data());
				c.sbuf.clear();
			}
			else status = H5Dwrite(c.dset, c.memtype, mspace, fspace, H5P_DEFAULT, c.buf.data());
			if (status < 0) {
				H5Sclose(fspace);
				H5Sclose(mspace);
				throw std::runtime_error("Error writing HDF5 output file.");
			}
			H5Sclose(fspace);
			c.buf.clear();
		}
		H5Sclose(mspace);
		nWritten += nBuffered;
		nBuffered = 0;
	}

	void close() {
		if (file < 0) return;
		flush();
		for (auto& c : cols) {
			if (c.dset >= 0) H5Dclose(c.dset);
			if (c.isString) H5Tclose(c.memtype);
		}
		cols.clear();
		H5Fclose(file);
		file = -1;
	}

	~H5Backend() {
		try { close(); }
		catch (...) {}
	}
};
#else
struct OutputTable::H5Backend {}; // HDF5 support not compiled in
#endif // RS_HDF5

//---------------------------------------------------------------------------
// OutputTable
//---------------------------------------------------------------------------
OutputTable::OutputTable() {}

OutputTable::~OutputTable() { close(); }

void OutputTable::clearColumns() {
	if (opened) throw std::logic_error("OutputTable: cannot change columns of an open table.");
	columns.clear();
}

void OutputTable::addColumn(const std::string& name, ColType type) {
	if (opened) throw std::logic_error("OutputTable: cannot add a column to an open table.");
	columns.push_back(Column{ name, type });
}

bool OutputTable::open(const std::string& basePath, OutputFormat format) {
	close();
	fmt = format;
	colIx = 0;
	if (fmt == OutputFormat::TEXT) {
		ofs.open((basePath + ".txt").c_str());
		if (!ofs.is_open()) return false;
		for (std::size_t i = 0; i < columns.size(); i++) {
			if (i > 0) ofs << '\t';
			ofs << columns[i].name;
		}
		ofs << '\n';
		opened = true;
		return true;
	}
#ifdef RS_HDF5
	h5 = std::make_unique<H5Backend>();
	if (!h5->create(basePath + ".h5", columns)) {
		h5.reset();
		return false;
	}
	opened = true;
	return true;
#else
	return false; // binary output requested but not compiled in
#endif
}

void OutputTable::close() {
	if (!opened) return;
	if (fmt == OutputFormat::TEXT) {
		ofs.close();
		ofs.clear();
	}
#ifdef RS_HDF5
	else if (h5) {
		h5->close();
		h5.reset();
	}
#endif
	opened = false;
	colIx = 0;
}

void OutputTable::checkColumn() {
	if (colIx >= columns.size())
		throw std::logic_error("OutputTable: more values than columns in a row.");
}

template <typename T>
void OutputTable::writeText(T v) {
	if (colIx > 0) ofs << '\t';
	ofs << v;
	colIx++;
}

void OutputTable::writeBinaryString(const std::string& v) {
#ifdef RS_HDF5
	checkColumn();
	if (columns[colIx].type != STRING)
		throw std::logic_error("OutputTable: text value for numeric column " + columns[colIx].name);
	h5->appendString(colIx, v);
	colIx++;
#else
	(void)v;
#endif
}

template <typename T>
void OutputTable::writeBinary(T v) {
#ifdef RS_HDF5
	checkColumn();
	if (columns[colIx].type == STRING)
		throw std::logic_error("OutputTable: numeric value for text column " + columns[colIx].name);
	h5->append(colIx, columns[colIx].type, v);
	colIx++;
#else
	(void)v;
#endif
}

OutputTable& OutputTable::operator<<(int v) {
	if (fmt == OutputFormat::TEXT) writeText(v); else writeBinary(v);
	return *this;
}
OutputTable& OutputTable::operator<<(unsigned int v) {
	if (fmt == OutputFormat::TEXT) writeText(v); else writeBinary(v);
	return *this;
}
OutputTable& OutputTable::operator<<(long v) {
	if (fmt == OutputFormat::TEXT) writeText(v); else writeBinary(v);
	return *this;
}
OutputTable& OutputTable::operator<<(float v) {
	if (fmt == OutputFormat::TEXT) writeText(v); else writeBinary(v);
	return *this;
}
OutputTable& OutputTable::operator<<(double v) {
	if (fmt == OutputFormat::TEXT) writeText(v); else writeBinary(v);
	return *this;
}
OutputTable& OutputTable::operator<<(const std::string& v) {
	if (fmt == OutputFormat::TEXT) writeText(v); else writeBinaryString(v);
	return *this;
}

void OutputTable::endRow() {
	if (fmt == OutputFormat::TEXT) {
		ofs << '\n';
	}
#ifdef RS_HDF5
	else {
		if (colIx != columns.size())
			throw std::logic_error("OutputTable: row has fewer values than columns.");
		h5->rowCompleted();
	}
#endif
	colIx = 0;
}
