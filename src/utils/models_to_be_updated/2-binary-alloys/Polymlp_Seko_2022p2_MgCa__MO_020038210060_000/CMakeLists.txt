#
# Author: Atsuto Seko
#

cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_PREFIX_PATH $ENV{KIM_API_CMAKE_PREFIX_DIR})
find_package(KIM-API-ITEMS 2.2 REQUIRED CONFIG)

kim_api_items_setup_before_project(ITEM_TYPE "portableModel")
project(Polymlp_Seko_2022p2_MgCa__MO_020038210060_000)
kim_api_items_setup_after_project(ITEM_TYPE "portableModel")

add_kim_api_model_library(
  NAME            ${PROJECT_NAME}
  DRIVER_NAME     "Polymlp__MD_000000123456_000"
  PARAMETER_FILES "polymlp.settings"
                  "polymlp.yaml"
  )
