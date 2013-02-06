/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkSubSampleTubeSpatialObjectFilter_txx
#define __itkSubSampleTubeSpatialObjectFilter_txx

#include "itkSubSampleTubeSpatialObjectFilter.h"

namespace itk
{

template< typename TTubeSpatialObject >
SubSampleTubeSpatialObjectFilter< TTubeSpatialObject >
::SubSampleTubeSpatialObjectFilter():
  m_Sampling(1)
{
}


template< typename TTubeSpatialObject >
void
SubSampleTubeSpatialObjectFilter< TTubeSpatialObject >
::GenerateData()
{
  TubeSpatialObjectType * output = this->GetOutput();
  const TubeSpatialObjectType * input = this->GetInput();

  typedef typename TubeSpatialObjectType::PointListType PointListType;
  const PointListType & inputPoints = input->GetPoints();
  PointListType & outputPoints = output->GetPoints();
  const size_t numberOfInputPoints = inputPoints.size();
  size_t numberOfOutputPoints;
  if( this->m_Sampling == 1 )
    {
    numberOfOutputPoints = numberOfInputPoints / this->m_Sampling + 0;
    }
  else if( numberOfInputPoints % this->m_Sampling == 0 )
    {
    numberOfOutputPoints = numberOfInputPoints / this->m_Sampling + 1;
    }
  else
    {
    numberOfOutputPoints = numberOfInputPoints / this->m_Sampling + 2;
    }
  outputPoints.resize( numberOfOutputPoints );
  for( size_t inputIndex = 0, outputIndex = 0;
    outputIndex < numberOfOutputPoints - 1;
    ++outputIndex, inputIndex += this->m_Sampling )
    {
    outputPoints[outputIndex] = inputPoints[inputIndex];
    }
  outputPoints[numberOfOutputPoints - 1] = inputPoints[numberOfInputPoints - 1];
}

} // end namespace itk

#endif