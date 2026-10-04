if(CAIRO_INCLUDE_DIRS AND CAIRO_LIBRARIES)
  set(CAIRO_FOUND TRUE)
else()
  if(NOT WIN32)
    find_package(PkgConfig REQUIRED)
    if(Cairo_FIND_VERSION_COUNT GREATER 0)
      set(_cairo_version_cmp ">=${Cairo_FIND_VERSION}")
    endif()
    pkg_check_modules(_pc_cairo cairo${_cairo_version_cmp})
  else()
    set(_pc_cairo_FOUND TRUE)
  endif()

  if(_pc_cairo_FOUND)
    set(CAIRO_FOUND FALSE)
    find_library(CAIRO_LIBRARY cairo HINTS ${_pc_cairo_LIBRARY_DIRS})
    set(CAIRO_LIBRARIES "${CAIRO_LIBRARY}")
    find_path(CAIRO_INCLUDE_DIR cairo.h HINTS ${_pc_cairo_INCLUDE_DIRS} PATH_SUFFIXES cairo)
    set(CAIRO_INCLUDE_DIRS "${CAIRO_INCLUDE_DIR}")
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Cairo DEFAULT_MSG CAIRO_LIBRARIES CAIRO_INCLUDE_DIRS)
  endif()
endif()

mark_as_advanced(CAIRO_CFLAGS CAIRO_INCLUDE_DIRS CAIRO_LIBRARIES)

if(CAIRO_FOUND AND NOT TARGET Cairo::Cairo)
  add_library(Cairo::Cairo UNKNOWN IMPORTED)
  set_target_properties(Cairo::Cairo PROPERTIES
    IMPORTED_LOCATION "${CAIRO_LIBRARIES}"
    INTERFACE_INCLUDE_DIRECTORIES "${CAIRO_INCLUDE_DIRS}"
  )
endif()
