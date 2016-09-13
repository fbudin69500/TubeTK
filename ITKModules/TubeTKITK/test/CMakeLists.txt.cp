# This file is empty in order to satisty ITKModuleExternal.cmake requirements
# when building against a version of ITK 4.9 that does not include the
# following change :
# InsightSoftwareConsortium/ITK@2c51c92

include_regular_expression( "^.*$" )

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( CompareImages_EXE
 ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )
 
 set( CompareTextFiles_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareTextFiles> )

######################  CompareTrainingMaskPython  ########################3
if( ITK_WRAP_PYTHON )
  set( MODULE_NAME ComputeTrainingMaskPython )
  # Test1 - VesselMask
  Midas3FunctionAddTestWithEnv( NAME ${MODULE_NAME}-Test1
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}Test.py
      MIDAS{ComputeTrainingMask-Test1.mha.md5}
      ${TEMP}/${MODULE_NAME}_VesselMaskTest1.mha
      ${TEMP}/${MODULE_NAME}_NotVesselMaskTest1.mha
      0.5
      2
    ENVIRONMENT ITK_BUILD_DIR=${ITK_DIR} TubeTK_BUILD_DIR=${PROJECT_BINARY_DIR}
      )

  # Test1 - Compare - VesselMask
  Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-VesselMask-Compare
    COMMAND ${CompareImages_EXE}
      -t ${TEMP}/${MODULE_NAME}_VesselMaskTest1.mha
      -b MIDAS{ComputeTrainingMask-Test1_VesselMask.mha.md5}
      -i 0.01 )

  set_property( TEST ${MODULE_NAME}-Test1-VesselMask-Compare
    APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

######################  ConvertImagesToCSV  ########################
  set( MODULE_NAME ConvertImagesToCSVPython )
  # Test1 
  Midas3FunctionAddTestWithEnv( NAME ${MODULE_NAME}-Test1
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}Test.py
      MIDAS{GDS0015_Large-TrainingMask.mha.md5}
      MIDAS{GDS0015_Large.mha.md5},MIDAS{ES0015_Large.mha.md5}
      ${TEMP}/${MODULE_NAME}Test1.csv
    ENVIRONMENT ITK_BUILD_DIR=${ITK_DIR} TubeTK_BUILD_DIR=${PROJECT_BINARY_DIR}
      )

  # Test1 - Compare 
  Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-Compare
    COMMAND ${CompareTextFiles_EXE}
      -t ${TEMP}/${MODULE_NAME}Test1.csv
	  -b MIDAS{ConvertImagesToCSVTest1.csv.md5} )
  set_property( TEST ${MODULE_NAME}-Test1-Compare
    APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )
endif()
