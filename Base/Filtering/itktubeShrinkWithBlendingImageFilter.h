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
/*=========================================================================
*
*  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
*
*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
*
*  For complete copyright, license and disclaimer of warranty information
*  please refer to the NOTICE file at the top of the ITK source tree.
*
*=========================================================================*/
#ifndef __itktubeShrinkWithBlendingImageFilter_h
#define __itktubeShrinkWithBlendingImageFilter_h

#include "itkShrinkImageFilter.h"

// Forward declare itkTubeTK class to allow friendship
namespace tube
{
template< typename tube_TInputImage, typename tube_TOutputImage >
class ShrinkWithBlendingImage;
}

namespace itk {

namespace tube {

/** \class ShrinkWithBlendingImageFilter
* \brief Reduce the size of an image by an integer factor in each
* dimension.
*
* ShrinkWithBlendingImageFilter reduces the size of an image by an integer
* factor in each dimension. The algorithm implemented is a max over the
* subsample. The output image size in each dimension is given by:
*
* outputSize[j] = max( std::floor(inputSize[j]/shrinkFactor[j]), 1 );
*
* NOTE: The physical centers of the input and output will be the
* same. Because of this, the Origin of the output may not be the same
* as the Origin of the input.
* Since this filter produces an image which is a different
* resolution, origin and with different pixel spacing than its input
* image, it needs to override several of the methods defined
* in ProcessObject in order to properly manage the pipeline execution
* model. In particular, this filter overrides
* ProcessObject::GenerateInputRequestedRegion() and
* ProcessObject::GenerateOutputInformation().
*
* This filter is implemented as a multithreaded filter.  It provides a
* ThreadedGenerateData() method for its implementation.
*
* \ingroup GeometricTransform Streamed
* \ingroup ITKImageGrid
*
*/
template< class TInputImage, class TOutputImage >
class ITK_EXPORT ShrinkWithBlendingImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef ShrinkWithBlendingImageFilter                   Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ShrinkWithBlendingImageFilter, ShrinkImageFilter );

  /** Typedef to images */
  typedef TOutputImage                          OutputImageType;
  typedef TInputImage                           InputImageType;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;

  typedef typename TOutputImage::IndexType      OutputIndexType;
  typedef typename TInputImage::IndexType       InputIndexType;
  typedef typename TOutputImage::OffsetType     OutputOffsetType;

  typedef typename TOutputImage::RegionType     OutputImageRegionType;
  typedef typename TInputImage::SizeType        InputSizeType;

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TInputImage::ImageDimension );

  itkStaticConstMacro( OutputImageDimension, unsigned int,
                       TOutputImage::ImageDimension );

  typedef Vector< float, ImageDimension >       PointImagePixelType;
  typedef Image< PointImagePixelType, OutputImageDimension >
                                                PointImageType;

  typedef FixedArray< unsigned int, ImageDimension > ShrinkFactorsType;

  /** Set the shrink factors. Values are clamped to
   * a minimum value of 1. Default is 1 for all dimensions. */
  void SetShrinkFactors( ShrinkFactorsType factors);
  void SetShrinkFactor(unsigned int i, unsigned int factor);
  void SetNewSize( InputSizeType newSize );
  itkGetMacro(NewSize, InputSizeType);

  /** Get the shrink factors. */
  itkGetConstReferenceMacro(ShrinkFactors, ShrinkFactorsType);

  itkSetMacro( Overlap, InputIndexType );
  itkGetMacro( Overlap, InputIndexType );

  itkSetMacro( BlendWithMean, bool );
  itkGetMacro( BlendWithMean, bool );

  itkSetMacro( BlendWithMax, bool );
  itkGetMacro( BlendWithMax, bool );

  itkSetMacro( BlendWithGaussianWeighting, bool );
  itkGetMacro( BlendWithGaussianWeighting, bool );

  itkSetMacro( UseLog, bool );
  itkGetMacro( UseLog, bool );

  itkSetMacro( UseNewSize, bool );
  itkGetMacro( UseNewSize, bool );
  itkBooleanMacro(UseNewSize);

  itkGetObjectMacro( PointImage, PointImageType );

  void GenerateOutputInformation( void );

  void GenerateInputRequestedRegion( void );

protected:
  ShrinkWithBlendingImageFilter( void );
  ~ShrinkWithBlendingImageFilter( void ) {}
  void PrintSelf( std::ostream & os, Indent indent ) const;

  void ThreadedGenerateData( const OutputImageRegionType &
    outputRegionForThread, ThreadIdType threadId );
  void UpdateShrinkFactors();

private:
  ShrinkWithBlendingImageFilter( const Self & ); //purposely not implemented
  void operator=( const Self & );            //purposely not implemented

  typename PointImageType::Pointer  m_PointImage;

  typename TInputImage::IndexType   m_Overlap;

  bool m_UseLog;

  bool m_BlendWithMean;
  bool m_BlendWithMax;
  bool m_BlendWithGaussianWeighting;
  bool m_UseNewSize;

  ShrinkFactorsType                  m_ShrinkFactors;
  InputSizeType                      m_NewSize;

  /** friendship facilitating itkTukeTK integration */
  template< typename tube_TInputImage, typename tube_TOutputImage >
  friend class ::tube::ShrinkWithBlendingImage;
};

} // end namespace tube

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeShrinkWithBlendingImageFilter.hxx"
#endif

#endif
