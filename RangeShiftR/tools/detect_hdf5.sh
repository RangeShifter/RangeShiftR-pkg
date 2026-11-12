#!/bin/sh
# Detects the optional HDF5 C library for RangeShiftR and writes src/<makevars>
# from src/<makevars>.in. Called by ./configure (Linux, macOS) and
# ./configure.win (Windows, Rtools).
#
# usage: sh tools/detect_hdf5.sh <makevars-name> <platform>
#        e.g.  sh tools/detect_hdf5.sh Makevars unix
#              sh tools/detect_hdf5.sh Makevars.win windows
#
# Without HDF5 the package is built exactly as before (text output only); the
# installation never fails because HDF5 is missing.
#
# Environment variables:
#   RANGESHIFTR_HDF5=no          do not use HDF5 even if it is available
#   HDF5_CPPFLAGS / HDF5_LIBS    use these flags instead of auto-detection

MAKEVARS="$1"
PLATFORM="$2"

: ${R_HOME=`R RHOME`}
if [ "$PLATFORM" = "windows" ]; then
  RBIN="${R_HOME}/bin${R_ARCH_BIN}/R.exe"
else
  RBIN="${R_HOME}/bin/R"
fi
CXX=`"$RBIN" CMD config CXX20 2>/dev/null`
[ -z "$CXX" ] && CXX=`"$RBIN" CMD config CXX`

# candidate flags, tried in this order until one compiles and links
CANDIDATES=""
add_candidate() { CANDIDATES="${CANDIDATES}${1}|${2}
"; }

if [ -n "${HDF5_LIBS}" ]; then
  add_candidate "${HDF5_CPPFLAGS}" "${HDF5_LIBS}"
fi
if pkg-config --exists hdf5 2>/dev/null; then
  if [ "$PLATFORM" = "windows" ]; then
    add_candidate "`pkg-config --cflags hdf5`" "`pkg-config --libs --static hdf5`"
  fi
  add_candidate "`pkg-config --cflags hdf5`" "`pkg-config --libs hdf5`"
fi
if command -v h5cc >/dev/null 2>&1; then
  H5_PREFIX=`h5cc -showconfig 2>/dev/null | sed -n 's/^ *Installation point: *//p'`
  [ -n "$H5_PREFIX" ] && add_candidate "-I${H5_PREFIX}/include" "-L${H5_PREFIX}/lib -lhdf5"
fi
if [ "$PLATFORM" = "windows" ]; then
  # the HDF5 library of Rtools (static) and its dependencies
  add_candidate "" "-lhdf5 -lz -lm -lsz"
  add_candidate "" "-lhdf5 -lsz -laec -lz -lm"
  add_candidate "" "-lhdf5 -lsz -laec -lz -lm -lshlwapi"
else
  # Debian/Ubuntu install the serial HDF5 library in a sub folder
  add_candidate "-I/usr/include/hdf5/serial" "-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5"
  add_candidate "" "-lhdf5"
fi

H5_CPPFLAGS=""
H5_LIBS=""
FOUND=no
if [ "${RANGESHIFTR_HDF5}" = "no" ]; then
  REASON=disabled
else
  tmpdir=`mktemp -d 2>/dev/null || echo ./rs_h5_test_$$`
  mkdir -p "$tmpdir"
  cat > "$tmpdir/h5test.cpp" <<'EOT'
#include <hdf5.h>
int main() { unsigned maj, min, rel; return H5get_libversion(&maj, &min, &rel) < 0; }
EOT
  OLDIFS="$IFS"
  IFS='
'
  for cand in $CANDIDATES; do
    IFS="$OLDIFS"
    c_flags=`echo "$cand" | cut -d'|' -f1`
    c_libs=`echo "$cand" | cut -d'|' -f2`
    if $CXX $c_flags "$tmpdir/h5test.cpp" -o "$tmpdir/h5test" $c_libs >/dev/null 2>&1; then
      H5_CPPFLAGS="-DRS_HDF5 $c_flags"
      H5_LIBS="$c_libs"
      FOUND=yes
      break
    fi
    IFS='
'
  done
  IFS="$OLDIFS"
  rm -rf "$tmpdir"
  REASON=notfound
fi

if [ "$FOUND" = "yes" ]; then
  echo "------------------------------------------------------------------"
  echo " RangeShiftR: HDF5 library found."
  echo " Binary output files are available: Simulation(OutputFormat = \"hdf5\")"
  echo "------------------------------------------------------------------"
elif [ "$REASON" = "disabled" ]; then
  echo "RangeShiftR: HDF5 support disabled (RANGESHIFTR_HDF5=no); text output only."
else
  echo "******************************************************************"
  echo " NOTE: The HDF5 library (libhdf5) was not found on this system."
  echo " RangeShiftR is installed WITHOUT HDF5 functionality:"
  echo " output files are written as text (OutputFormat = \"txt\") only."
  echo " Everything else works as usual."
  echo ""
  echo " To install RangeShiftR with HDF5 functionality, install libhdf5"
  echo " and then reinstall RangeShiftR:"
  case "$PLATFORM" in
    windows)
      echo "   Rtools normally includes HDF5; please update Rtools." ;;
    *)
      if [ "`uname -s 2>/dev/null`" = "Darwin" ]; then
        echo "   macOS:          brew install hdf5"
      else
        echo "   Debian/Ubuntu:  sudo apt install libhdf5-dev"
        echo "   Fedora/RHEL:    sudo dnf install hdf5-devel"
        echo "   conda:          conda install -c conda-forge hdf5"
      fi ;;
  esac
  echo "******************************************************************"
fi

sed -e "s|@HDF5_CPPFLAGS@|${H5_CPPFLAGS}|" -e "s|@HDF5_LIBS@|${H5_LIBS}|" \
  "src/${MAKEVARS}.in" > "src/${MAKEVARS}"
exit 0
