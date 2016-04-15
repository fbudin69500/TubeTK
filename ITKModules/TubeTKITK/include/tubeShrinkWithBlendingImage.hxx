/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __tubeShrinkWithBlendingImage_hxx
#define __tubeShrinkWithBlendingImage_hxx

#include "tubeShrinkWithBlendingImage.h"

namespace tube {

/**
 *
 */
template< class TInputImage, class TOutputImage >
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::ShrinkWithBlendingImage( void )
{
  m_ShrinkWithBlendingFilter = ShrinkWithBlendingFilterType::New();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, itk::Indent  indent ) const
{
  m_ShrinkWithBlendingFilter->PrintSelf( os, indent);
}


/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetShrinkFactors(ShrinkFactorsType shrinkFactors)
{
  m_ShrinkWithBlendingFilter->SetShrinkFactors(shrinkFactors);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetShrinkFactor(unsigned int i, unsigned int factor)
{
  m_ShrinkWithBlendingFilter->SetShrinkFactor(i,factor);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetNewSize(InputSizeType newSize)
{
  m_ShrinkWithBlendingFilter->SetNewSize(newSize);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
typename ShrinkWithBlendingImage< TInputImage, TOutputImage >::InputSizeType
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetNewSize(void)
{
  return m_ShrinkWithBlendingFilter->GetNewSize();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetUseNewSize( bool useNewSize)
{
  m_ShrinkWithBlendingFilter->SetUseNewSize(useNewSize);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
bool
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetUseNewSize()
{
  return m_ShrinkWithBlendingFilter->GetUseNewSize();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
typename ShrinkWithBlendingImage< TInputImage, TOutputImage >::ShrinkFactorsType
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetShrinkFactors(void)
{
  return m_ShrinkWithBlendingFilter->GetShrinkFactors();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetOverlap(InputIndexType overlap )
{
  m_ShrinkWithBlendingFilter->SetOverlap(overlap);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
typename ShrinkWithBlendingImage< TInputImage, TOutputImage >::InputIndexType
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetOverlap(void)
{
  return m_ShrinkWithBlendingFilter->GetOverlap();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetBlendWithMean( bool blendWithMean )
{
  m_ShrinkWithBlendingFilter->SetBlendWithMean(blendWithMean);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
bool
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetBlendWithMean(void)
{
  return m_ShrinkWithBlendingFilter->GetBlendWithMean();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetBlendWithMax( bool blendWithMax)
{
  m_ShrinkWithBlendingFilter->SetBlendWithMax(blendWithMax);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
bool
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetBlendWithMax(void)
{
  return m_ShrinkWithBlendingFilter->GetBlendWithMax();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetBlendWithGaussianWeighting( bool blendWithGaussianWeighting)
{
  m_ShrinkWithBlendingFilter->SetBlendWithGaussianWeighting(blendWithGaussianWeighting);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
bool
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetBlendWithGaussianWeighting(void)
{
  return m_ShrinkWithBlendingFilter->GetBlendWithGaussianWeighting();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetUseLog( bool useLog)
{
  m_ShrinkWithBlendingFilter->SetUseLog(useLog);
}
/**
 *
 */
template< class TInputImage, class TOutputImage >
bool
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetUseLog(void)
{
  return m_ShrinkWithBlendingFilter->GetUseLog();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
typename ShrinkWithBlendingImage<TInputImage, TOutputImage>::PointImageType *
ShrinkWithBlendingImage<TInputImage, TOutputImage>
::GetPointImage(void)
{
  return m_ShrinkWithBlendingFilter->GetPointImage();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GenerateOutputInformation( void )
{
  m_ShrinkWithBlendingFilter->GenerateOutputInformation();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GenerateInputRequestedRegion( void )
{
  m_ShrinkWithBlendingFilter->GenerateInputRequestedRegion();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetInput( const TInputImage *inputImage )
{
  m_ShrinkWithBlendingFilter->SetInput(inputImage);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::Update()
{
  m_ShrinkWithBlendingFilter->Update();
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
typename TOutputImage::Pointer
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetOutput()
{
  return m_ShrinkWithBlendingFilter->GetOutput();
}

} // end namespace tube

#endif
