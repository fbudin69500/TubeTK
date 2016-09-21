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

#ifndef __itktubeConvertImagesToCSVFilter_hxx
#define __itktubeConvertImagesToCSVFilter_hxx

#include "itktubeConvertImagesToCSVFilter.h"

namespace itk
{

namespace tube
{

/** Constructor */
template< class TInputImage >
ConvertImagesToCSVFilter< TInputImage >
::ConvertImagesToCSVFilter( void )
{
  m_InputImage = NULL;
  m_NumberRows = 0;
}

template< class TInputImage >
void
ConvertImagesToCSVFilter< TInputImage >
::Update( void )
{
  // number of pixels/stride or less
  const unsigned int ARows =
    m_InputImage->GetLargestPossibleRegion().GetNumberOfPixels() / m_Stride;
  const unsigned int ACols = m_ImageList.size() + 1;
  m_Output.set_size(ARows, ACols);
  std::vector< IteratorType * > iterList;
  for (unsigned int i = 0; i < m_NumImages; ++i)
    {
    iterList.push_back(new IteratorType(m_ImageList[i],
    m_ImageList[i]->GetLargestPossibleRegion()));
    }
  IteratorType maskIter(m_InputImage, m_InputImage->GetLargestPossibleRegion());
  unsigned int i = 0;
  while (!maskIter.IsAtEnd())
    {
    if (maskIter.Get() != 0)
      {
      for (i = 0; i<m_NumImages; ++i)
        {
        m_Output(m_NumberRows, i) = iterList[i]->Get();
        }
      m_Output(m_NumberRows, i) = maskIter.Get();
      m_NumberRows++;
      }
    for (int s = 0; s<m_Stride && !maskIter.IsAtEnd(); ++s)
      {
      for (unsigned int i = 0; i<m_NumImages; ++i)
        {
        ++(*iterList[i]);
        }
      ++maskIter;
      }
    }
  for (unsigned int i = 0; i<iterList.size(); ++i)
    {
    delete iterList[i];
    }
  iterList.clear();
}

template< class TInputImage >
void
ConvertImagesToCSVFilter< TInputImage >::
SetNthInput(InputImageType* image)
{
  this->m_ImageList.push_back (image);
}


/** PrintSelf */
template< class TInputImage>
void
ConvertImagesToCSVFilter< TInputImage >::
PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Stride = " << m_Stride << std::endl;
  os << indent << "NumImages = " << m_NumImages << std::endl;
}

} // End namespace tube

} // End namespace itk
#endif // End !defined( __itktubeConvertImagesToCSVFilter_hxx )
