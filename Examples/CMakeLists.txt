##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

#
# Setup examples as tests
#
if( TubeTK_USE_PYTHON )
  TubeTKCheckPythonLibraries( REQUIRED LIBRARIES IPython tornado zmq jinja2
    tables matplotlib ERROR_MESSAGE
    "Set TubeTK_USE_EXAMPLES_AS_TESTS to OFF" )

  # Add notebook/example directories
  add_subdirectory( MergeAdjacentImages )

  if( TubeTK_USE_NUMPY_STACK )
    add_subdirectory( TubeNumPyArrayAndPropertyHistograms )
  endif()
endif()
