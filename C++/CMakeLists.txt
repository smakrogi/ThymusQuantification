#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.4)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
#-----------------------------------------------------------------------------
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

#-----------------------------------------------------------------------------
# Module name
# - Usually identique to the directory name)
# - You shouldn't call 'project(...)'. See SlicerConfig.cmake for more details.
#-----------------------------------------------------------------------------
set(EXTENSION_NAME ThymusQuantification)
set(EXTENSION2_NAME MaskOverlapAnalysis)

#-----------------------------------------------------------------------------
# Prerequisites
#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)

  # Additional C/CXX flags - Should be defined before including Slicer_USE_FILE
  set(ADDITIONAL_C_FLAGS "" CACHE STRING "Additional ${EXTENSION_NAME} C Flags")
  set(ADDITIONAL_CXX_FLAGS "" CACHE STRING "Additional ${EXTENSION_NAME} CXX Flags")

  # set(ADDITIONAL_C_FLAGS "" CACHE STRING "Additional ${EXTENSION2_NAME} C Flags")
  # set(ADDITIONAL_CXX_FLAGS "" CACHE STRING "Additional ${EXTENSION2_NAME} CXX Flags")

  include(${Slicer_USE_FILE})

  set(EXTENSION_DEPENDS NA) # Specified as a space separated list or 'NA' if any
  set(EXTENSION_LICENSE_FILE ${Slicer_LICENSE_FILE})
  set(EXTENSION_README_FILE ${Slicer_README_FILE})

  include(SlicerEnableExtensionTesting)
endif()

# Extension meta information
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Slicer4:Developers:Projects:QtSlicer/Tutorials/ExtensionWriting")
set(EXTENSION_CATEGORY "NIH-NIA-3T")
set(EXTENSION_CONTRIBUTOR "Sokratis Makrogiannis")
set(EXTENSION_MAJOR_VERSION 1)
set(EXTENSION_MINOR_VERSION 0)
set(EXTENSION_PATCH_VERSION 0)
set(EXTENSION_STATUS "Beta")
set(EXTENSION_DESCRIPTION "This is an example of CLI module built as an extension")
set(EXTENSION_ACKNOWLEDGEMENTS "This work was supported by NIH/NIA/IRP")
set(EXTENSION_LICENSE_SHORT_DESCRIPTION ${Slicer_LICENSE_SHORT_DESCRIPTION})

# Extension2 meta information
set(EXTENSION2_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Slicer4:Developers:Projects:QtSlicer/Tutorials/ExtensionWriting")
set(EXTENSION2_CATEGORY "NIH-NIA-3T")
set(EXTENSION2_CONTRIBUTOR "Sokratis Makrogiannis")
set(EXTENSION2_MAJOR_VERSION 1)
set(EXTENSION2_MINOR_VERSION 0)
set(EXTENSION2_PATCH_VERSION 0)
set(EXTENSION2_STATUS "Beta")
set(EXTENSION2_DESCRIPTION "This is an example of CLI module built as an extension")
set(EXTENSION2_ACKNOWLEDGEMENTS "This work was supported by NIH/NIA/IRP")
set(EXTENSION2_LICENSE_SHORT_DESCRIPTION ${Slicer_LICENSE_SHORT_DESCRIPTION})

# Configure XML description
# configure_file(
#   ${EXTENSION_NAME}.xml.in
#   ${CMAKE_CURRENT_BINARY_DIR}/${EXTENSION_NAME}.xml
#   @ONLY
#   )

# configure_file(
#   ${EXTENSION2_NAME}.xml.in
#   ${CMAKE_CURRENT_BINARY_DIR}/${EXTENSION2_NAME}.xml
#   @ONLY
#   )

# Add module sub libraries
#add_subdirectory(Logic)

set(cli_module_export_directive "Q_SLICER_QTMODULES_LOADABLEEXTENSIONTEMPLATE_EXPORT")

# Additional includes
set(cli_module_include_directories
  )

# Additional source files
set(cli_module_SRCS
  ThymusSegmentationAlgorithms.cxx 
  OptionList.cxx 
  CommonFunctions.cxx
  )

set(cli_module_SRCS2
  ComputeOverlap.cxx
  )
  
# Additional Target libraries
set(cli_module_target_libraries
  ${ITK_LIBRARIES}
  )

#-----------------------------------------------------------------------------
# Build
#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${EXTENSION_NAME}
  #LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
  TARGET_LIBRARIES ${cli_module_target_libraries}
  INCLUDE_DIRECTORIES ${cli_module_include_directories}
  ADDITIONAL_SRCS ${cli_module_SRCS}
  #EXECUTABLE_ONLY
  )

SEMMacroBuildCLI(
  NAME ${EXTENSION2_NAME}
 #LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
  TARGET_LIBRARIES ${cli_module_target_libraries}
  INCLUDE_DIRECTORIES ${cli_module_include_directories}
  ADDITIONAL_SRCS ${cli_module_SRCS2}
  #EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if(BUILD_TESTING)
# add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
# Generate extension description file '<EXTENSION_NAME>.s4ext'
#-----------------------------------------------------------------------------
# if(NOT Slicer_SOURCE_DIR)
#   include(SlicerFunctionGenerateExtensionDescription)
#   slicerFunctionGenerateExtensionDescription(
#     EXTENSION_NAME ${EXTENSION_NAME}
#     EXTENSION_CATEGORY ${EXTENSION_CATEGORY}
#     EXTENSION_STATUS ${EXTENSION_STATUS}
#     EXTENSION_HOMEPAGE ${EXTENSION_HOMEPAGE}
#     EXTENSION_DESCRIPTION ${EXTENSION_DESCRIPTION}
#     EXTENSION_DEPENDS ${EXTENSION_DEPENDS}
#     DESTINATION_DIR ${CMAKE_CURRENT_BINARY_DIR}
#     SLICER_WC_REVISION ${Slicer_WC_REVISION}
#     SLICER_WC_ROOT ${Slicer_WC_ROOT}
#     )
# # endif()

# # if(NOT Slicer_SOURCE_DIR)
# #   include(SlicerFunctionGenerateExtensionDescription)
#   slicerFunctionGenerateExtensionDescription(
#     EXTENSION_NAME ${EXTENSION2_NAME}
#     EXTENSION_CATEGORY ${EXTENSION2_CATEGORY}
#     EXTENSION_STATUS ${EXTENSION2_STATUS}
#     EXTENSION_HOMEPAGE ${EXTENSION2_HOMEPAGE}
#     EXTENSION_DESCRIPTION ${EXTENSION2_DESCRIPTION}
#     EXTENSION_DEPENDS ${EXTENSION_DEPENDS}
#     DESTINATION_DIR ${CMAKE_CURRENT_BINARY_DIR}
#     SLICER_WC_REVISION ${Slicer_WC_REVISION}
#     SLICER_WC_ROOT ${Slicer_WC_ROOT}
#     )
# endif()

#-----------------------------------------------------------------------------
# Packaging
#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  include(${Slicer_EXTENSION_CPACK})
endif()
