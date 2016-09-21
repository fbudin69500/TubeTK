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
#ifndef __itktubeConvertImagesToCSVFilter_h
#define __itktubeConvertImagesToCSVFilter_h

#include <itkObject.h>
#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>

#include "tubeMessage.h"

namespace itk
{
namespace tube
{
/** \class ConvertImagesToCSV
 */

 template< class TInputImage >
class ConvertImagesToCSVFilter : public Object
{
public:
  /** Standard class typedefs. */
  typedef ConvertImagesToCSVFilter      Self;
  typedef Object                        Superclass;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;


  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      ImagePointer;
  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename InputImageType::IndexType    IndexType;
  typedef ImageFileReader< InputImageType >     ReaderType;
  typedef ImageRegionIterator< InputImageType > IteratorType;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ConvertImagesToCSVFilter, Object);

  /** ImageDimension constants */
  itkStaticConstMacro( ImageDimension, unsigned int,
                      TInputImage::ImageDimension );

  itkSetObjectMacro(InputImage, InputImageType);
  itkGetObjectMacro(InputImage, InputImageType);

  itkGetMacro(Output, vnl_matrix <InputPixelType>);

  itkGetMacro(Stride, unsigned int);
  itkSetMacro(Stride, unsigned int);
  itkGetMacro(NumImages, unsigned int);
  itkSetMacro(NumImages, unsigned int);
  itkSetMacro(NumberRows, unsigned int);
  itkGetMacro(NumberRows, unsigned int);
  void SetNthInput(InputImageType*);

  void Update( void );

protected:

  ConvertImagesToCSVFilter ( void );
  ~ConvertImagesToCSVFilter(void) {};
  //void GenerateData( void );
  void PrintSelf( std::ostream& os, Indent indent ) const;


private:
  ConvertImagesToCSVFilter ( const Self& );
  void operator=( const Self& );

  typename InputImageType::Pointer                 m_InputImage;
  vnl_matrix <InputPixelType>                      m_Output;
  std::vector< typename InputImageType::Pointer >  m_ImageList;
  int                                              m_Stride;
  unsigned int                                     m_NumImages;
  unsigned int                                     m_NumberRows;
  
}; // End class ConvertImagesToCSVFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeConvertImagesToCSVFilter.hxx"
#endif

#endif // End !defined(_itktubeConvertImagesToCSVFilter_h)
