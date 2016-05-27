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
#ifndef __tubeImageToImageRegistrationHelper_h
#define __tubeImageToImageRegistrationHelper_h

#include "itktubeImageToImageRegistrationHelper.h"
#include "itkObject.h"
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class ImageToImageRegistrationHelper
 *
 *  \ingroup TubeTKITK
 */


template< class TImage, class TTransformPrecisionType = double >
class ImageToImageRegistrationHelper:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ImageToImageRegistrationHelper                  Self;
  typedef itk::SmartPointer< Self >                       Pointer;
  typedef itk::SmartPointer< const Self >                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ImageToImageRegistrationHelper, Object );


  /** Typedef to images */
  typedef TImage                           ImageType;
  typedef TTransformPrecisionType          TransformPrecisionType;

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TImage::ImageDimension );

  typedef itk::tube::ImageToImageRegistrationHelper< ImageType, TransformPrecisionType >
  FilterType;
  //
  // Typedefs for the parameters of the registration methods
  //
  typedef typename  FilterType::MaskObjectType             MaskObjectType;
  typedef typename FilterType::PointType                   PointType;
  typedef typename FilterType::InterpolationMethodEnumType InterpolationMethodEnumType;
  //
  // Available Registration Methods
  //
  typedef typename FilterType::OptimizedRegistrationMethodType
  OptimizedRegistrationMethodType;

  typedef typename FilterType::MatrixTransformType             MatrixTransformType;
  typedef typename FilterType::BSplineTransformType            BSplineTransformType;
  typedef typename FilterType::PixelType                       PixelType;
  typedef typename FilterType::TransformListType               TransformListType;
  typedef typename FilterType::DisplacementFieldType           DisplacementFieldType;

  typedef typename FilterType::LandmarkVectorType              LandmarkVectorType;
  typedef typename FilterType::RigidTransformType              RigidTransformType;
  typedef typename FilterType::AffineTransformType             AffineTransformType;

  typedef typename FilterType::MetricMethodEnumType            MetricMethodEnumType;
  typedef typename FilterType::InitialMethodEnumType           InitialMethodEnumType;


  // **************
  // **************
  //  Specify the fixed and moving images
  // **************
  // **************
  tubeWrapSetConstObjectMacro( FixedImage, ImageType, RegistrationHelper );
  tubeWrapGetConstObjectMacro( FixedImage, ImageType, RegistrationHelper );

  tubeWrapSetConstObjectMacro( MovingImage, ImageType, RegistrationHelper );
  tubeWrapGetConstObjectMacro( MovingImage, ImageType, RegistrationHelper );

  // **************
  // **************
  tubeWrapSetMacro( RandomNumberSeed, unsigned int, RegistrationHelper );
  tubeWrapGetMacro( RandomNumberSeed, unsigned int, RegistrationHelper );
  // **************
  // **************
  //  Specify how the fixed image should be sampled when computing the metric and
  //    what ROI of the moving image is valid
  // **************
  // **************
  tubeWrapSetMacro( UseFixedImageMaskObject, bool, RegistrationHelper );
  tubeWrapGetConstMacro( UseFixedImageMaskObject, bool, RegistrationHelper );
  tubeWrapBooleanMacro( UseFixedImageMaskObject, RegistrationHelper );

  tubeWrapSetConstObjectMacro( FixedImageMaskObject, MaskObjectType,
                               RegistrationHelper );

  tubeWrapSetMacro( UseRegionOfInterest, bool, RegistrationHelper );
  tubeWrapGetMacro( UseRegionOfInterest, bool, RegistrationHelper );
  tubeWrapSetMacro( RegionOfInterestPoint1, PointType, RegistrationHelper );
  tubeWrapGetMacro( RegionOfInterestPoint1, PointType, RegistrationHelper );
  tubeWrapSetMacro( RegionOfInterestPoint2, PointType, RegistrationHelper );
  tubeWrapGetMacro( RegionOfInterestPoint2, PointType, RegistrationHelper );

  void SetRegionOfInterest( const PointType & point1, const PointType & point2 );

  void SetRegionOfInterest( const std::vector<float> & points );
  // **************
  //  Initialize the moving image mask as the region of initial overlap
  //  between the fixed and moving images
  // **************
  tubeWrapSetMacro( SampleFromOverlap, bool, RegistrationHelper );
  tubeWrapGetMacro( SampleFromOverlap, bool, RegistrationHelper );
  tubeWrapBooleanMacro( SampleFromOverlap, RegistrationHelper );

  tubeWrapSetMacro( SampleIntensityPortion, double, RegistrationHelper );
  tubeWrapGetConstMacro( SampleIntensityPortion, double, RegistrationHelper );
  // **************
  // **************
  //  Update
  // **************
  // **************
  tubeWrapCallMacro( Initialize, RegistrationHelper );
  /** This class provides an Update() method to fit the appearance of a
   * ProcessObject API, but it is not a ProcessObject.  */
  tubeWrapCallMacro( Update, RegistrationHelper );

  typename ImageType::ConstPointer  ResampleImage(
    InterpolationMethodEnumType interp
      = OptimizedRegistrationMethodType
        ::LINEAR_INTERPOLATION, const ImageType * movingImage = NULL,
    const MatrixTransformType * matrixTransform = NULL,
    const BSplineTransformType * bsplineTransform = NULL,
    PixelType defaultPixelValue = 0);

  // Returns the moving image resampled into the space of the fixed image
  typename ImageType::ConstPointer GetFinalMovingImage(InterpolationMethodEnumType interp
                                                       = OptimizedRegistrationMethodType
                                                         ::LINEAR_INTERPOLATION);

  // **************
  // **************
  // Compute registration "accuracy" by comparing a resampled moving image
  // with a baseline image.
  // **************
  // **************

  // Specify the baseline image.
  tubeWrapSetConstObjectMacro( BaselineImage, ImageType,
                               RegistrationHelper );
  // Bound the required accuracy for the registration test to "pass"
  tubeWrapSetMacro( BaselineNumberOfFailedPixelsTolerance, unsigned int,
                    RegistrationHelper );
  tubeWrapSetMacro( BaselineIntensityTolerance, PixelType, RegistrationHelper );
  tubeWrapSetMacro( BaselineRadiusTolerance, unsigned int, RegistrationHelper );

  // Must be called after setting the BaselineImage in order to resample
  //   the moving image into the BaselineImage space, compute differences,
  //   and determine if it passed the test within the specified tolerances
  tubeWrapCallMacro( ComputeBaselineDifference, RegistrationHelper );

  tubeWrapGetConstObjectMacro( BaselineDifferenceImage, ImageType, RegistrationHelper );
  tubeWrapGetConstObjectMacro( BaselineResampledMovingImage, ImageType,
                               RegistrationHelper );
  tubeWrapGetMacro( BaselineNumberOfFailedPixels, unsigned int, RegistrationHelper );
  tubeWrapGetMacro( BaselineTestPassed, bool, RegistrationHelper );

  // **************
  // **************
  // Process Control
  // **************
  // **************

  // **************
  // Control which steps of the registration pipeline are applied
  // **************

  tubeWrapSetMacro( EnableLoadedRegistration, bool, RegistrationHelper );
  tubeWrapGetConstMacro( EnableLoadedRegistration, bool, RegistrationHelper );
  tubeWrapBooleanMacro( EnableLoadedRegistration, RegistrationHelper );

  tubeWrapSetMacro( EnableInitialRegistration, bool, RegistrationHelper );
  tubeWrapGetConstMacro( EnableInitialRegistration, bool, RegistrationHelper );
  tubeWrapBooleanMacro( EnableInitialRegistration, RegistrationHelper );

  tubeWrapSetMacro( EnableRigidRegistration, bool, RegistrationHelper );
  tubeWrapGetConstMacro( EnableRigidRegistration, bool, RegistrationHelper );
  tubeWrapBooleanMacro( EnableRigidRegistration, RegistrationHelper );

  tubeWrapSetMacro( EnableAffineRegistration, bool, RegistrationHelper );
  tubeWrapGetConstMacro( EnableAffineRegistration, bool, RegistrationHelper );
  tubeWrapBooleanMacro( EnableAffineRegistration, RegistrationHelper );

  tubeWrapSetMacro( EnableBSplineRegistration, bool, RegistrationHelper );
  tubeWrapGetConstMacro( EnableBSplineRegistration, bool, RegistrationHelper );
  tubeWrapBooleanMacro( EnableBSplineRegistration, RegistrationHelper );

  // **************
  // Specify the optimizer
  // **************
  tubeWrapSetMacro( UseEvolutionaryOptimization, bool, RegistrationHelper );
  tubeWrapGetMacro( UseEvolutionaryOptimization, bool, RegistrationHelper );
  // **************
  // Specify the expected magnitudes within the transform.  Used to
  //   guide the operating space of the optimizers
  // **************
  tubeWrapSetMacro( ExpectedOffsetPixelMagnitude, double, RegistrationHelper );
  tubeWrapGetConstMacro( ExpectedOffsetPixelMagnitude, double, RegistrationHelper );

  tubeWrapSetMacro( ExpectedRotationMagnitude, double, RegistrationHelper);
  tubeWrapGetConstMacro( ExpectedRotationMagnitude, double, RegistrationHelper );

  tubeWrapSetMacro( ExpectedScaleMagnitude, double, RegistrationHelper );
  tubeWrapGetConstMacro( ExpectedScaleMagnitude, double, RegistrationHelper );

  tubeWrapSetMacro( ExpectedSkewMagnitude, double, RegistrationHelper );
  tubeWrapGetConstMacro( ExpectedSkewMagnitude, double, RegistrationHelper );
  // **************
  //  Return the current product of the registration pipeline
  // **************
  tubeWrapGetConstObjectMacro( CurrentMatrixTransform, MatrixTransformType,
                               RegistrationHelper );
  tubeWrapGetConstObjectMacro( CurrentBSplineTransform, BSplineTransformType,
                               RegistrationHelper );

  // The image used for registration is updated at certain points in the
  //   registration pipeline for speed and transform composition.
  // Specifically, the image is resmpled using the loaded transforms prior
  //   to running the initial registration method and the image is resampled
  //   after the affine registration / prior to running bspline registration
  // The result of these resamplings is available as the CurrentMovingImage.
  tubeWrapGetConstObjectMacro( CurrentMovingImage, ImageType, RegistrationHelper );
  tubeWrapGetConstObjectMacro( LoadedTransformResampledImage, ImageType,
                               RegistrationHelper );
  tubeWrapGetConstObjectMacro( MatrixTransformResampledImage, ImageType,
                               RegistrationHelper );
  tubeWrapGetConstObjectMacro( BSplineTransformResampledImage, ImageType,
                               RegistrationHelper );
  // **************
  //  Final metric value after the pipeline has completed
  // **************
  tubeWrapGetMacro( FinalMetricValue, double, RegistrationHelper );
  // **************
  //  Determine if progress messages should be sent to cout
  // **************
  tubeWrapSetMacro( ReportProgress, bool, RegistrationHelper );
  tubeWrapGetMacro( ReportProgress, bool, RegistrationHelper );
  tubeWrapBooleanMacro( ReportProgress, RegistrationHelper );

  tubeWrapSetMacro( MinimizeMemory, bool, RegistrationHelper );
  tubeWrapGetMacro( MinimizeMemory, bool, RegistrationHelper );
  tubeWrapBooleanMacro( MinimizeMemory, RegistrationHelper );
  //
  // Loaded transforms parameters
  //
  void SetLoadedMatrixTransform( const MatrixTransformType & tfm );
  tubeWrapGetConstObjectMacro( LoadedMatrixTransform, MatrixTransformType,
                               RegistrationHelper );
  void SetLoadedBSplineTransform( const BSplineTransformType & tfm );
  tubeWrapGetConstObjectMacro( LoadedBSplineTransform, BSplineTransformType,
                               RegistrationHelper );

  tubeWrapSetMacro( TransformList, TransformListType, RegistrationHelper );
  tubeWrapGetMacro( TransformList, TransformListType, RegistrationHelper );
  tubeWrapGetObjectMacro( DisplacementField, DisplacementFieldType, RegistrationHelper );

  //
  // Initial Parameters
  //
  void SetFixedLandmarks( const LandmarkVectorType & fixedLandmarks );
  void SetMovingLandmarks( const LandmarkVectorType & movingLandmarks );

  tubeWrapSetMacro( InitialMethodEnum, InitialMethodEnumType, RegistrationHelper );
  tubeWrapGetConstMacro( InitialMethodEnum, InitialMethodEnumType, RegistrationHelper );

  //
  // Rigid Parameters
  //

  tubeWrapSetMacro( RigidSamplingRatio, double, RegistrationHelper );
  tubeWrapGetConstMacro( RigidSamplingRatio, double, RegistrationHelper );

  tubeWrapSetMacro( RigidTargetError, double, RegistrationHelper );
  tubeWrapGetConstMacro( RigidTargetError, double, RegistrationHelper );

  tubeWrapSetMacro( RigidMaxIterations, unsigned int, RegistrationHelper );
  tubeWrapGetConstMacro( RigidMaxIterations, unsigned int, RegistrationHelper );

  tubeWrapSetMacro( RigidMetricMethodEnum, MetricMethodEnumType, RegistrationHelper );
  tubeWrapGetConstMacro( RigidMetricMethodEnum, MetricMethodEnumType,
                         RegistrationHelper );

  tubeWrapSetMacro( RigidInterpolationMethodEnum, InterpolationMethodEnumType,
                    RegistrationHelper );
  tubeWrapGetConstMacro( RigidInterpolationMethodEnum, InterpolationMethodEnumType,
                         RegistrationHelper );

  tubeWrapGetConstObjectMacro( RigidTransform, RigidTransformType,
                               RegistrationHelper );
  tubeWrapGetMacro( RigidMetricValue, double, RegistrationHelper );

  //
  // Affine Parameters
  //

  tubeWrapSetMacro( AffineSamplingRatio, double, RegistrationHelper );
  tubeWrapGetConstMacro( AffineSamplingRatio, double, RegistrationHelper );

  tubeWrapSetMacro( AffineTargetError, double, RegistrationHelper );
  tubeWrapGetConstMacro( AffineTargetError, double, RegistrationHelper );

  tubeWrapSetMacro( AffineMaxIterations, unsigned int, RegistrationHelper );
  tubeWrapGetConstMacro( AffineMaxIterations, unsigned int, RegistrationHelper );

  tubeWrapSetMacro( AffineMetricMethodEnum, MetricMethodEnumType, RegistrationHelper );
  tubeWrapGetConstMacro( AffineMetricMethodEnum, MetricMethodEnumType,
                         RegistrationHelper );

  tubeWrapSetMacro( AffineInterpolationMethodEnum, InterpolationMethodEnumType,
                    RegistrationHelper );
  tubeWrapGetConstMacro( AffineInterpolationMethodEnum, InterpolationMethodEnumType,
                         RegistrationHelper );

  tubeWrapGetConstObjectMacro( AffineTransform, AffineTransformType,
                               RegistrationHelper );
  tubeWrapGetMacro( AffineMetricValue, double, RegistrationHelper );

  //
  // BSpline Parameters
  //
  tubeWrapSetMacro( BSplineSamplingRatio, double, RegistrationHelper );
  tubeWrapGetConstMacro( BSplineSamplingRatio, double, RegistrationHelper );

  tubeWrapSetMacro( BSplineTargetError, double, RegistrationHelper );
  tubeWrapGetConstMacro( BSplineTargetError, double, RegistrationHelper );

  tubeWrapSetMacro( BSplineMaxIterations, unsigned int, RegistrationHelper );
  tubeWrapGetConstMacro( BSplineMaxIterations, unsigned int, RegistrationHelper );

  tubeWrapSetMacro( BSplineControlPointPixelSpacing, double, RegistrationHelper );
  tubeWrapGetConstMacro( BSplineControlPointPixelSpacing, double, RegistrationHelper );

  tubeWrapSetMacro( BSplineMetricMethodEnum, MetricMethodEnumType, RegistrationHelper );
  tubeWrapGetConstMacro( BSplineMetricMethodEnum, MetricMethodEnumType,
                         RegistrationHelper );

  tubeWrapSetMacro( BSplineInterpolationMethodEnum, InterpolationMethodEnumType,
                    RegistrationHelper );
  tubeWrapGetConstMacro( BSplineInterpolationMethodEnum, InterpolationMethodEnumType,
                         RegistrationHelper );

  tubeWrapGetConstObjectMacro( BSplineTransform, BSplineTransformType,
                               RegistrationHelper );
  tubeWrapGetMacro( BSplineMetricValue, double, RegistrationHelper );

protected:
  ImageToImageRegistrationHelper( void );
  ~ImageToImageRegistrationHelper() {}

  void PrintSelfHelper( std::ostream & os, itk::Indent indent,
    const std::string & basename, MetricMethodEnumType metric,
    InterpolationMethodEnumType interpolation ) const;
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkImageToImageRegistrationHelperFilter parameters **/
  ImageToImageRegistrationHelper(const Self &);
  void operator=(const Self &);


  typename FilterType::Pointer m_RegistrationHelper;

};
} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeImageToImageRegistrationHelper.hxx"
#endif

#endif // End !defined( __tubeImageToImageRegistrationHelper_h )
