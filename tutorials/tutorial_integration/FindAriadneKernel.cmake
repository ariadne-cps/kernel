find_library(ARIADNE_KERNEL_LIBRARY NAMES ariadne-kernel)

find_package(Threads REQUIRED)
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(GMP QUIET gmp)
  pkg_check_modules(MPFR QUIET mpfr)
  pkg_check_modules(CAIRO QUIET cairo)
endif()

find_path(GMP_INCLUDE_DIR gmp.h HINTS ${GMP_INCLUDE_DIRS})
find_path(MPFR_INCLUDE_DIR mpfr.h HINTS ${MPFR_INCLUDE_DIRS})
find_path(CAIRO_INCLUDE_DIR cairo.h HINTS ${CAIRO_INCLUDE_DIRS} PATH_SUFFIXES cairo)
find_library(GMP_LIBRARY NAMES gmp libgmp HINTS ${GMP_LIBRARY_DIRS})
find_library(MPFR_LIBRARY NAMES mpfr libmpfr HINTS ${MPFR_LIBRARY_DIRS})
find_library(CAIRO_LIBRARY NAMES cairo libcairo HINTS ${CAIRO_LIBRARY_DIRS})

find_path(ARIADNE_KERNEL_INCLUDE_DIR ariadne-kernel.hpp PATH_SUFFIXES ariadne-kernel)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneKernel DEFAULT_MSG
  ARIADNE_KERNEL_LIBRARY
  ARIADNE_KERNEL_INCLUDE_DIR
  GMP_INCLUDE_DIR
  MPFR_INCLUDE_DIR
  CAIRO_INCLUDE_DIR
  GMP_LIBRARY
  MPFR_LIBRARY
  CAIRO_LIBRARY
)

if(AriadneKernel_FOUND)
  get_filename_component(ARIADNE_KERNEL_INCLUDE_PARENT_DIR ${ARIADNE_KERNEL_INCLUDE_DIR} DIRECTORY)
  set(ARIADNE_KERNEL_INCLUDE_DIRS
    ${ARIADNE_KERNEL_INCLUDE_PARENT_DIR}
    ${ARIADNE_KERNEL_INCLUDE_DIR}
    ${MPFR_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR}
    ${CAIRO_INCLUDE_DIR}
  )
  if(NOT TARGET AriadneKernel::ariadne-kernel)
    add_library(AriadneKernel::ariadne-kernel UNKNOWN IMPORTED)
    set_target_properties(AriadneKernel::ariadne-kernel PROPERTIES
      IMPORTED_LOCATION "${ARIADNE_KERNEL_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${ARIADNE_KERNEL_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${CAIRO_LIBRARY};${MPFR_LIBRARY};${GMP_LIBRARY};Threads::Threads"
    )
  endif()
endif()
