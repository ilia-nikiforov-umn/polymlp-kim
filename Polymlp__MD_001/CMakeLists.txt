# Copyright (c) 2026 Atsuto Seko
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#
# Required preamble
#

cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_PREFIX_PATH $ENV{KIM_API_CMAKE_PREFIX_DIR})
find_package(KIM-API-ITEMS 2.2 REQUIRED CONFIG)

# End preamble

kim_api_items_setup_before_project(ITEM_TYPE "modelDriver")
project(Polymlp__MD_001 LANGUAGES CXX)
kim_api_items_setup_after_project(ITEM_TYPE "modelDriver")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT DEFINED CMAKE_CXX_FLAGS_RELEASE)
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -Wno-deprecated -Wsign-compare -DNDEBUG")
endif()

add_kim_api_model_driver_library(
  NAME                    ${PROJECT_NAME}
  CREATE_ROUTINE_NAME     "model_driver_create"
  CREATE_ROUTINE_LANGUAGE "cpp"
  )


target_sources(${PROJECT_NAME} PRIVATE
  polymlp/polymlp_api.cpp
  polymlp/polymlp_basis_function.cpp
  polymlp/polymlp_features.cpp
  polymlp/polymlp_features_polynomial.cpp
  polymlp/polymlp_features_mapping_utils.cpp
  polymlp/polymlp_features_utils.cpp
  polymlp/polymlp_functions_interface.cpp
  polymlp/polymlp_mapping.cpp
  polymlp/polymlp_mlipkk_spherical_harmonics.cpp
  polymlp/polymlp_model_params.cpp
  polymlp/polymlp_model_params_gtinv.cpp
  polymlp/polymlp_model_params_polynomial.cpp
  polymlp/polymlp_potential.cpp
  polymlp/polymlp_products.cpp
  polymlp/polymlp_read_gtinv.cpp
  polymlp/polymlp_gtinv_data.cpp
  polymlp/polymlp_gtinv_data_order1.cpp
  polymlp/polymlp_gtinv_data_order2.cpp
  polymlp/polymlp_gtinv_data_order3.cpp
  polymlp/polymlp_gtinv_data_order4.cpp
  polymlp/polymlp_gtinv_data_order5.cpp
  polymlp/polymlp_gtinv_data_order6.cpp
#  polymlp_kim.cpp
#  model_driver_polymlp.cpp
)
