# Copyright (c) 2026 Atsuto Seko
#
#
# Required preamble
#

cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_PREFIX_PATH $ENV{KIM_API_CMAKE_PREFIX_DIR})
find_package(KIM-API-ITEMS 2.2 REQUIRED CONFIG)

# End preamble

kim_api_items_setup_before_project(ITEM_TYPE "modelDriver")
project(Polymlp__MD_000000123456_000 LANGUAGES CXX)
kim_api_items_setup_after_project(ITEM_TYPE "modelDriver")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
find_package(OpenMP)

add_kim_api_model_driver_library(
  NAME                    ${PROJECT_NAME}
  CREATE_ROUTINE_NAME     "model_driver_create"
  CREATE_ROUTINE_LANGUAGE "cpp"
  )


target_sources(${PROJECT_NAME} PRIVATE
    polymlp_kim.cpp
    model_driver_polymlp.cpp
    polymlp/polymlp_api.cpp
    polymlp/polymlp_basis_function.cpp
    polymlp/polymlp_features.cpp
    polymlp/polymlp_features_mapping_utils.cpp
    polymlp/polymlp_features_polynomial.cpp
    polymlp/polymlp_features_utils.cpp
    polymlp/polymlp_functions_interface.cpp
    polymlp/polymlp_mapping.cpp
    polymlp/polymlp_mlipkk_spherical_harmonics.cpp
    polymlp/polymlp_model_params.cpp
    polymlp/polymlp_model_params_gtinv.cpp
    polymlp/polymlp_model_params_polynomial.cpp
    polymlp/polymlp_parse_potential.cpp
    polymlp/polymlp_parse_potential_legacy.cpp
    polymlp/polymlp_potential.cpp
    polymlp/polymlp_products.cpp
    polymlp/polymlp_read_gtinv.cpp
    polymlp/polymlp_gtinv_data.cpp
    polymlp/polymlp_gtinv_data_order3.cpp
    polymlp/polymlp_gtinv_data_order4.cpp
    polymlp/polymlp_gtinv_data_order5.cpp
    polymlp/polymlp_gtinv_data_order6.cpp
)

if(OpenMP_CXX_FOUND)
  target_link_libraries(
      ${PROJECT_NAME}
      PRIVATE OpenMP::OpenMP_CXX
  )
endif()
