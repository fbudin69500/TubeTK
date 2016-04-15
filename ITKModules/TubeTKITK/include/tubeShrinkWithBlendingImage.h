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
#ifndef __tubeShrinkWithBlendingImage_h
#define __tubeShrinkWithBlendingImage_h

#include "itktubeShrinkWithBlendingImageFilter.h"
#include "itkObject.h"


namespace tube
{
/** \class ShrinkWithBlendingImage
 *
 *  \ingroup TubeTKITK
 */

template< typename TInputImage, typename TOutputImage >
class ShrinkWithBlendingImage:
  public itk::ProcessObject
{
public:
    /** Standard class typedefs. */
    typedef ShrinkWithBlendingImage                         Self;
    typedef itk::SmartPointer< Self >                            Pointer;
    typedef itk::SmartPointer< const Self >                      ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( ShrinkWithBlendingImage, Object );


    /** Typedef to images */
    typedef TOutputImage                          OutputImageType;
    typedef TInputImage                           InputImageType;
//    typedef typename OutputImageType::Pointer     OutputImagePointer;
//    typedef typename InputImageType::Pointer      InputImagePointer;
//    typedef typename InputImageType::ConstPointer InputImageConstPointer;

//    typedef typename TOutputImage::IndexType      OutputIndexType;
    typedef typename TInputImage::IndexType       InputIndexType;
//    typedef typename TOutputImage::OffsetType     OutputOffsetType;

    itkStaticConstMacro( ImageDimension, unsigned int,
                         TInputImage::ImageDimension );

    itkStaticConstMacro( OutputImageDimension, unsigned int,
                         TOutputImage::ImageDimension );

    typedef itk::Vector< float, ImageDimension >       PointImagePixelType;
    typedef itk::Image< PointImagePixelType, OutputImageDimension >
                                                  PointImageType;

    typedef itk::FixedArray< unsigned int, ImageDimension > ShrinkFactorsType;

    /** Set the shrink factors. Values are clamped to
     * a minimum value of 1. Default is 1 for all dimensions. */
    void SetShrinkFactors(ShrinkFactorsType ShrinkFactors);
    void SetShrinkFactor(unsigned int i, unsigned int factor);

    /** Get the shrink factors. */
    ShrinkFactorsType GetShrinkFactors();

    void SetOverlap(InputIndexType overlap );
    InputIndexType  GetOverlap();

    void SetBlendWithMean( bool blendWithMean );
    bool GetBlendWithMean();

    void SetBlendWithMax( bool blendWithMax);
    bool GetBlendWithMax();

    void SetBlendWithGaussianWeighting( bool blendWithGaussianWeighting);
    bool GetBlendWithGaussianWeighting();

    void SetUseLog( bool useLog);
    bool GetUseLog();

    PointImageType* GetPointImage();

    void GenerateOutputInformation( void );

    void GenerateInputRequestedRegion( void );

    void SetInput( const InputImageType *inputImage );
    void Update(void);
    typename OutputImageType::Pointer GetOutput(void);

protected:
  ShrinkWithBlendingImage( void );
  ~ShrinkWithBlendingImage() {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
  /** itkShrinkWithBlendingImageFilter parameters **/
  ShrinkWithBlendingImage(const Self &);
  void operator=(const Self &);

  typedef itk::tube::ShrinkWithBlendingImageFilter< InputImageType, OutputImageType > ShrinkWithBlendingFilterType;
  typename ShrinkWithBlendingFilterType::Pointer m_ShrinkWithBlendingFilter;

};
} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeShrinkWithBlendingImage.hxx"
#endif

#endif // End !defined( __tubeShrinkWithBlendingImage_h )
