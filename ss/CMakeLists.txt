if(NOT TARGET bdd)
    ensure_submodule_initialized(bdd)
    add_subdirectory(bdd)
endif()

if(NOT TARGET algebra)
    ensure_submodule_initialized(algebra)
    add_subdirectory(algebra)
endif()
