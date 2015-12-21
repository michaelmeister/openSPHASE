# Copyright 2011 Gregor Burger
# tested for linux only

find_library(TBB_TBB_LIBRARY NAMES tbb)
find_library(TBB_MALLOC_LIBRARY NAMES tbbmalloc)
set(TBB_LIBRARIES ${TBB_TBB_LIBRARY} ${TBB_MALLOC_LIBRARY})
find_path(TBB_INCLUDE_DIR NAMES tbb/tbb.h)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TBB DEFAULT_MSG TBB_LIBRARIES TBB_INCLUDE_DIR)
