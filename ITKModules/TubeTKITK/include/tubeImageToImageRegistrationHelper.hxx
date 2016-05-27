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
#ifndef __tubeImageToImageRegistrationHelper_hxx
#define __tubeImageToImageRegistrationHelper_hxx

#include "tubeImageToImageRegistrationHelper.h"

namespace tube {

/**
 *
 */
template< class TImage, class TTransformPrecisionType >
ImageToImageRegistrationHelper<TImage,TTransformPrecisionType>
::ImageToImageRegistrationHelper( void )
{
  m_RegistrationHelper = FilterType::New();
}

/**
 *
 */
template< class TImage, class TTransformPrecisionType >
void
ImageToImageRegistrationHelper<TImage,TTransformPrecisionType>
::SetRegionOfInterest( const PointType & point1, const PointType & point2 )
{
  m_RegistrationHelper->SetRegionOfInterest(point1,point2);
}

/**
 *
 */
template< class TImage, class TTransformPrecisionType >
void
ImageToImageRegistrationHelper<TImage,TTransformPrecisionType>
::SetRegionOfInterest( const std::vector<float> & points )
{
  m_RegistrationHelper->SetRegionOfInterest(points);
}

/**
 *
 */
template< class TImage, class TTransformPrecisionType >
typename ImageToImageRegistrationHelper<TImage,TTransformPrecisionType>
::ImageType::ConstPointer
ImageToImageRegistrationHelper<TImage,TTransformPrecisionType>
::ResampleImage(
    InterpolationMethodEnumType interp, const TImage * movingImage,
    const MatrixTransformType * matrixTransform,
    const BSplineTransformType * bsplineTransform,
    PixelType defaultPixelValue)
{
  m_RegistrationHelper->ResampleImage( interp, movingImage,
                       matrixTransform, bsplineTransform, defaultPixelValue );
}

/**
 *
 */
template< class TImage, class TTransformPrecisionType >
typename ImageToImageRegistrationHelper<TImage,TTransformPrecisionType>
::ImageType::ConstPointer
ImageToImageRegistrationHelper<TImage,TTransformPrecisionType>
::GetFinalMovingImage(InterpolationMethodEnumType interp)
{
  m_RegistrationHelper->GetFinalMovingImage(interp);
}

/**
 *
 */
template< class TImage, class TTransformPrecisionType >
void
ImageToImageRegistrationHelper<TImage,TTransformPrecisionType>
::SetFixedLandmarks( const LandmarkVectorType & fixedLandmarks )
{
  m_RegistrationHelper->SetFixedLandmarks(fixedLandmarks);
  this-> Modified();
}

/**
 *
 */
template< class TImage, class TTransformPrecisionType >
void
ImageToImageRegistrationHelper<TImage,TTransformPrecisionType>
::SetMovingLandmarks( const LandmarkVectorType & movingLandmarks )
{
  m_RegistrationHelper->SetMovingLandmarks(movingLandmarks);
  this-> Modified();
}

/**
 *
 */
template< class TImage, class TTransformPrecisionType >
void
ImageToImageRegistrationHelper<TImage, TTransformPrecisionType>
::SetLoadedMatrixTransform( const MatrixTransformType & tfm )
{
  m_RegistrationHelper->SetLoadedMatrixTransform(tfm);
  this-> Modified();
}

/**
 *
 */
template< class TImage, class TTransformPrecisionType >
void
ImageToImageRegistrationHelper<TImage, TTransformPrecisionType>
::SetLoadedBSplineTransform( const BSplineTransformType & tfm )
{
  m_RegistrationHelper->SetLoadedBSplineTransform(tfm);
  this-> Modified();
}

/**
 *
 */
template <class TImage, class TTransformPrecisionType>
void
ImageToImageRegistrationHelper<TImage, TTransformPrecisionType>
::PrintSelfHelper( std::ostream & os, itk::Indent indent,
                   const std::string & basename,
                   MetricMethodEnumType metric,
                   InterpolationMethodEnumType interpolation ) const
{
  switch( metric )
    {
    case OptimizedRegistrationMethodType::MATTES_MI_METRIC:
      os << indent << basename << " Metric Method = MATTES_MI_METRIC"
        << std::endl;
      break;
    case OptimizedRegistrationMethodType::NORMALIZED_CORRELATION_METRIC:
      os << indent << basename
        << " Metric Method = CROSS_CORRELATION_METRIC" << std::endl;
      break;
    case OptimizedRegistrationMethodType::MEAN_SQUARED_ERROR_METRIC:
      os << indent << basename
        << " Metric Method = MEAN_SQUARED_ERROR_METRIC" << std::endl;
      break;
    default:
      os << indent << basename << " Metric Method = UNKNOWN" << std::endl;
      break;
    }
  os << indent << std::endl;

  switch( interpolation )
    {
    case OptimizedRegistrationMethodType::NEAREST_NEIGHBOR_INTERPOLATION:
      os << indent << basename
        << " Interpolation Method = NEAREST_NEIGHBOR_INTERPOLATION"
        << std::endl;
      break;
    case OptimizedRegistrationMethodType::LINEAR_INTERPOLATION:
      os << indent << basename
        << " Interpolation Method = LINEAR_INTERPOLATION" << std::endl;
      break;
    case OptimizedRegistrationMethodType::BSPLINE_INTERPOLATION:
      os << indent << basename
        << " Interpolation Method = BSPLINE_INTERPOLATION" << std::endl;
      break;
    case OptimizedRegistrationMethodType::SINC_INTERPOLATION:
      os << indent << basename
        << " Interpolation Method = SINC_INTERPOLATION" << std::endl;
      break;
    default:
      os << indent << basename
        << " Interpolation Method = UNKNOWN" << std::endl;
      break;
    }
}

/**
 *
 */
template <class TImage, class TTransformPrecisionType>
void
ImageToImageRegistrationHelper<TImage, TTransformPrecisionType>
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_RegistrationHelper->GetFixedImage() )
    {
    os << indent << "Fixed Image = "
       << m_RegistrationHelper->GetFixedImage() << std::endl;
    }
  if( m_RegistrationHelper->GetMovingImage() )
    {
    os << indent << "Moving Image = "
       << m_RegistrationHelper->GetMovingImage() << std::endl;
    }
  os << indent << std::endl;
  os << indent << "Use region of interest = "
     << m_RegistrationHelper->GetUseRegionOfInterest() << std::endl;
  os << indent << "Region of interest point1 = "
     << m_RegistrationHelper->GetRegionOfInterestPoint1()
     << std::endl;
  os << indent << "Region of interest point2 = "
     << m_RegistrationHelper->GetRegionOfInterestPoint2()
     << std::endl;
  os << indent << std::endl;
  os << indent << "Use Fixed Image Mask Object = "
     << m_RegistrationHelper->GetUseFixedImageMaskObject()
     << std::endl;
  os << indent << std::endl;
  if( m_RegistrationHelper->GetFixedImageMaskObject() )
    {
    os << indent << "Fixed Image Mask Object = "
       << m_RegistrationHelper->GetFixedImageMaskObject() << std::endl;
    }
  os << indent << "Use Moving Image Mask Object = "
     << m_RegistrationHelper->GetUseMovingImageMaskObject()
     << std::endl;
  os << indent << std::endl;
  if( m_RegistrationHelper->GetMovingImageMaskObject() )
    {
    os << indent << "Moving Image Mask Object = "
      << m_RegistrationHelper->GetMovingImageMaskObject() << std::endl;
    }
  os << indent << std::endl;
  os << indent << "Random Number Seed = "
     << m_RegistrationHelper->GetRandomNumberSeed()
     << std::endl;
  os << indent << std::endl;
  os << indent << "Enable Loaded Registration = "
     << m_RegistrationHelper->GetEnableLoadedRegistration()
     << std::endl;
  os << indent << "Enable Initial Registration = "
     << m_RegistrationHelper->GetEnableInitialRegistration()
     << std::endl;
  os << indent << "Enable Rigid Registration = "
     << m_RegistrationHelper->GetEnableRigidRegistration()
     << std::endl;
  os << indent << "Enable Affine Registration = "
     << m_RegistrationHelper->GetEnableAffineRegistration()
     << std::endl;
  os << indent << "Enable BSpline Registration = "
     << m_RegistrationHelper->GetEnableBSplineRegistration()
     << std::endl;
  os << indent << std::endl;
  os << indent << "Expected Offset (in Pixels) Magnitude = "
     << m_RegistrationHelper->GetExpectedOffsetPixelMagnitude()
     << std::endl;
  os << indent << "Expected Rotation Magnitude = "
     << m_RegistrationHelper->GetExpectedRotationMagnitude()
     << std::endl;
  os << indent << "Expected Scale Magnitude = "
     << m_RegistrationHelper->GetExpectedScaleMagnitude()
     << std::endl;
  os << indent << "Expected Skew Magnitude = "
     << m_RegistrationHelper->GetExpectedSkewMagnitude()
     << std::endl;
  os << indent << std::endl;
  os << indent << "Completed Initialization = "
     << m_RegistrationHelper->GetCompletedInitialization()
     << std::endl;
  os << indent << "Completed Resampling = "
     << m_RegistrationHelper->GetCompletedResampling()
     << std::endl;
  os << indent << std::endl;
  os << indent << "Rigid Metric Value = "
     << m_RegistrationHelper->GetRigidMetricValue()
     << std::endl;
  os << indent << "Affine Metric Value = "
     << m_RegistrationHelper->GetAffineMetricValue()
     << std::endl;
  os << indent << "BSpline Metric Value = "
     << m_RegistrationHelper->GetBSplineMetricValue()
     << std::endl;
  os << indent << "Final Metric Value = "
     << m_RegistrationHelper->GetFinalMetricValue()
     << std::endl;
  os << indent << std::endl;
  os << indent << "Report Progress = "
     << m_RegistrationHelper->GetReportProgress() << std::endl;
  os << indent << std::endl;
  if( m_RegistrationHelper->GetCurrentMovingImage() )
    {
    os << indent << "Current Moving Image = "
       << m_RegistrationHelper->GetCurrentMovingImage()
       << std::endl;
    }
  else
    {
    os << indent << "Current Moving Image = NULL" << std::endl;
    }
  if( m_RegistrationHelper->GetCurrentMatrixTransform() )
    {
    os << indent << "Current Matrix Transform = "
       << m_RegistrationHelper->GetCurrentMatrixTransform()
       << std::endl;
    }
  else
    {
    os << indent << "Current Matrix Transform = NULL" << std::endl;
    }
  if( m_RegistrationHelper->GetCurrentBSplineTransform() )
    {
    os << indent << "Current BSpline Transform = "
       << m_RegistrationHelper->GetCurrentBSplineTransform()
       << std::endl;
    }
  else
    {
    os << indent << "Current BSpline Transform = NULL" << std::endl;
    }
  os << indent << std::endl;
  if( m_RegistrationHelper->GetLoadedTransformResampledImage() )
    {
    os << indent << "Loaded Transform Resampled Image = "
       << m_RegistrationHelper->GetLoadedTransformResampledImage()
       << std::endl;
    }
  else
    {
    os << indent << "Loaded Transform Resampled Image = NULL" << std::endl;
    }
  if( m_RegistrationHelper->GetMatrixTransformResampledImage() )
    {
    os << indent << "Matrix Transform Resampled Image = "
       << m_RegistrationHelper->GetMatrixTransformResampledImage()
       << std::endl;
    }
  else
    {
    os << indent << "Matrix Transform Resampled Image = NULL" << std::endl;
    }
  if( m_RegistrationHelper->GetBSplineTransformResampledImage() )
    {
    os << indent << "BSpline Transform Resampled Image = "
       << m_RegistrationHelper->GetBSplineTransformResampledImage()
       << std::endl;
    }
  else
    {
    os << indent << "BSpline Transform Resampled Image = NULL"
       << std::endl;
    }
  os << indent << std::endl;
  if( m_RegistrationHelper->GetLoadedMatrixTransform() )
    {
    os << indent << "Loaded Matrix Transform = "
       << m_RegistrationHelper->GetLoadedMatrixTransform()
       << std::endl;
    }
  else
    {
    os << indent << "Loaded Matrix Transform = NULL" << std::endl;
    }
  if( m_RegistrationHelper->GetLoadedBSplineTransform() )
    {
    os << indent << "Loaded BSpline Transform = "
       << m_RegistrationHelper->GetLoadedBSplineTransform()
       << std::endl;
    }
  else
    {
    os << indent << "Loaded BSpline Transform = NULL" << std::endl;
    }
  os << indent << std::endl;

  switch( m_RegistrationHelper->GetInitialMethodEnum() )
    {
    case FilterType::INIT_WITH_NONE:
      os << indent << "Initial Registration Enum = INIT_WITH_NONE"
         << std::endl;
      break;
    case FilterType::INIT_WITH_CURRENT_RESULTS:
      os << indent
         << "Initial Registration Enum = INIT_WITH_CURRENT_RESULTS"
         << std::endl;
      break;
    case FilterType::INIT_WITH_IMAGE_CENTERS:
      os << indent
         << "Initial Registration Enum = INIT_WITH_IMAGE_CENTERS"
         << std::endl;
      break;
    case FilterType::INIT_WITH_CENTERS_OF_MASS:
      os << indent
         << "Initial Registration Enum = INIT_WITH_CENTERS_OF_MASS"
         << std::endl;
      break;
    case FilterType::INIT_WITH_SECOND_MOMENTS:
      os << indent
         << "Initial Registration Enum = INIT_WITH_SECOND_MOMENTS"
         << std::endl;
      break;
    default:
      os << indent << "Initial Registration Enum = UNKNOWN" << std::endl;
      break;
    }
  if( m_RegistrationHelper->GetInitialTransform() )
    {
    os << indent << "Initial Transform = "
       << m_RegistrationHelper->GetInitialTransform()
       << std::endl;
    }
  else
    {
    os << indent << "Initial Transform = NULL" << std::endl;
    }
  os << indent << std::endl;
  os << indent << "Rigid Sampling Ratio = "
     << m_RegistrationHelper->GetRigidSamplingRatio()
     << std::endl;
  os << indent << "Rigid Target Error = "
     << m_RegistrationHelper->GetRigidTargetError()
     << std::endl;
  os << indent << "Rigid Max Iterations = "
     << m_RegistrationHelper->GetRigidMaxIterations()
     << std::endl;
  PrintSelfHelper( os, indent, "Rigid",
                   m_RegistrationHelper->GetRigidMetricMethodEnum(),
            m_RegistrationHelper->GetRigidInterpolationMethodEnum() );
  os << indent << std::endl;
  if( m_RegistrationHelper->GetRigidTransform() )
    {
    os << indent << "Rigid Transform = "
       << m_RegistrationHelper->GetRigidTransform()
       << std::endl;
    }
  else
    {
    os << indent << "Rigid Transform = NULL" << std::endl;
    }
  os << indent << std::endl;
  os << indent << "Affine Sampling Ratio = "
     << m_RegistrationHelper->GetAffineSamplingRatio()
     << std::endl;
  os << indent << "Affine Target Error = "
     << m_RegistrationHelper->GetAffineTargetError()
     << std::endl;
  os << indent << "Affine Max Iterations = "
     << m_RegistrationHelper->GetAffineMaxIterations()
     << std::endl;
  PrintSelfHelper( os, indent, "Affine",
                   m_RegistrationHelper->GetAffineMetricMethodEnum(),
            m_RegistrationHelper->GetAffineInterpolationMethodEnum() );
  os << indent << std::endl;
  if( m_RegistrationHelper->GetAffineTransform() )
    {
    os << indent << "Affine Transform = "
       << m_RegistrationHelper->GetAffineTransform()
       << std::endl;
    }
  else
    {
    os << indent << "Affine Transform = NULL" << std::endl;
    }
  os << indent << std::endl;
  os << indent << "BSpline Sampling Ratio = "
     << m_RegistrationHelper->GetBSplineSamplingRatio()
     << std::endl;
  os << indent << "BSpline Target Error = "
     << m_RegistrationHelper->GetBSplineTargetError()
     << std::endl;
  os << indent << "BSpline Max Iterations = "
     << m_RegistrationHelper->GetBSplineMaxIterations()
     << std::endl;
  os << indent << "BSpline Control Point Pixel Spacing = "
     << m_RegistrationHelper->GetBSplineControlPointPixelSpacing()
     << std::endl;
  PrintSelfHelper( os, indent, "BSpline",
                   m_RegistrationHelper->GetBSplineMetricMethodEnum(),
           m_RegistrationHelper->GetBSplineInterpolationMethodEnum() );
  os << indent << std::endl;
  if( m_RegistrationHelper->GetBSplineTransform() )
    {
    os << indent << "BSpline Transform = "
       << m_RegistrationHelper->GetBSplineTransform()
       << std::endl;
    }
  else
    {
    os << indent << "BSpline Transform = NULL" << std::endl;
    }
  os << indent << std::endl;

}


} // end namespace tube

#endif
