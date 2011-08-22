/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSpatialObjectReader.h"
#include "itkGroupSpatialObject.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include "itkTubeRidgeExtractor.h"

int itkTubeRidgeExtractorTest2( int argc, char * argv[] )
  {
  if( argc != 3 )
    {
    std::cout
      << "itkTubeRidgeExtractorTest <inputImage> <vessel.tre>"
      << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<float, 3>   ImageType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName( argv[1] );
  imReader->Update();

  ImageType::Pointer im = imReader->GetOutput();

  typedef itk::tube::RidgeExtractor<ImageType> RidgeOpType;
  RidgeOpType::Pointer ridgeOp = RidgeOpType::New();

  ridgeOp->SetDebug( true );

  ridgeOp->SetInputImage( im );
  ridgeOp->SetStepX( 0.75 );
  ridgeOp->SetScale( 2.0 );
  ridgeOp->SetExtent( 2.5 );
  ridgeOp->SetDynamicScale( true );

  typedef itk::SpatialObjectReader<>                   ReaderType;
  typedef itk::SpatialObject<>::ChildrenListType       ObjectListType;
  typedef itk::GroupSpatialObject<>                    GroupType;
  typedef itk::VesselTubeSpatialObject<>               TubeType;
  typedef TubeType::PointListType                      PointListType;
  typedef TubeType::PointType                          PointType;
  typedef TubeType::TubePointType                      TubePointType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();
  GroupType::Pointer group = reader->GetGroup();

  std::cout << "Number of children = " << group->GetNumberOfChildren()
    << std::endl;

  char tubeName[17];
  strcpy( tubeName, "Tube" );
  ObjectListType * tubeList = group->GetChildren( -1, tubeName );

  unsigned int numTubes = tubeList->size();
  std::cout << "Number of tubes = " << numTubes << std::endl;

  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator
    RandGenType;
  RandGenType::Pointer rndGen = RandGenType::New();
  rndGen->Initialize(); // set seed here

  RidgeOpType::IndexType imMinX = ridgeOp->GetExtractBoundMin();
  RidgeOpType::IndexType imMaxX = ridgeOp->GetExtractBoundMax();

  int failures = 0;
  for( unsigned int mcRun=0; mcRun<2; mcRun++ )
    {
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "***** Beginning tube test ***** " << std::endl;
    std::cout << "***** Beginning tube test ***** " << std::endl;
    unsigned int rndTubeNum = rndGen->GetUniformVariate( 0, 1 ) * numTubes;
    if( rndTubeNum > numTubes-1 )
      {
      rndTubeNum = numTubes-1;
      }
    ObjectListType::iterator tubeIter = tubeList->begin();
    for( unsigned int i=0; i<rndTubeNum; i++ )
      {
      ++tubeIter;
      }
    TubeType::Pointer tube = static_cast< TubeType * >(
      tubeIter->GetPointer() );
    std::cout << "Test tube = " << rndTubeNum << std::endl;

    PointListType tubePointList = tube->GetPoints();
    unsigned int numPoints = tubePointList.size();
    unsigned int rndPointNum = rndGen->GetUniformVariate( 0, 1 )
      * numPoints;
    if( rndPointNum > numPoints-1 )
      {
      rndPointNum = numPoints-1;
      }
    PointListType::iterator pntIter = tubePointList.begin();
    for( unsigned int i=0; i<rndPointNum; i++ )
      {
      ++pntIter;
      }
    TubePointType * pnt = static_cast< TubePointType * >(&(*pntIter));
    std::cout << "Test point = " << rndPointNum << std::endl;

    RidgeOpType::ContinuousIndexType x0;
    for( unsigned int i=0; i<ImageType::ImageDimension; i++)
      {
      x0[i] = pnt->GetPosition()[i];
      }
    std::cout << "Test index = " << x0 << std::endl;

    std::cout << "Setting min and max z." << std::endl;
    RidgeOpType::IndexType minX = imMinX;
    RidgeOpType::IndexType maxX = imMaxX;
    if( x0[2] - 5 > minX[2] )
      {
      minX[2] = x0[2] - 5;
      }
    if( x0[2] + 5 < maxX[2] )
      {
      maxX[2] = x0[2] + 5;
      }
    ridgeOp->SetExtractBoundMin( minX );
    ridgeOp->SetExtractBoundMax( maxX );

    RidgeOpType::ContinuousIndexType x1 = x0;
    if( !ridgeOp->LocalRidge( x1 ) )
      {
      std::cout << "Local ridge test failed.  No ridge found." << std::endl;
      std::cout << "   Source = " << x0 << std::endl;
      std::cout << "   Result = " << x1 << std::endl;
      ++failures;
      continue;
      }

    double diff = 0;
    for( unsigned int i=0; i<ImageType::ImageDimension; i++)
      {
      double tf = x0[i]-x1[i];
      diff += tf * tf;
      }
    diff = vcl_sqrt( diff );
    if( diff > 2 )
      {
      std::cout << "Local ridge test failed.  Local ridge too far."
        << std::endl;
      std::cout << "   Source = " << x0 << std::endl;
      std::cout << "   Result = " << x1 << std::endl;
      ++failures;
      continue;
      }

    std::cout << "Local ridge discovery a success!" << std::endl;
    std::cout << std::endl;
    std::cout << "***** Beginning tube extraction ***** " << std::endl;
    TubeType::Pointer xTube = ridgeOp->ExtractRidge( x1, mcRun );
    std::cout << "***** Ending tube extraction ***** " << std::endl;

    if( xTube.IsNull() )
      {
      std::cout << "Ridge extraction failed" << std::endl;
      ++failures;
      continue;
      }

    TubeType::Pointer xTube2;
    xTube2 = ridgeOp->ExtractRidge( x1, 101 );
    if( xTube2.IsNotNull() )
      {
      std::cout << "Ridge extracted twice - test failed" << std::endl;
      ++failures;
      continue;
      }

    if( !ridgeOp->DeleteTube( xTube ) )
      {
      std::cout << "Delete tube failed" << std::endl;
      ++failures;
      continue;
      }

    xTube = ridgeOp->ExtractRidge( x1, 101 );
    if( xTube.IsNull() )
      {
      std::cout << "Ridge extraction after delete failed." << std::endl;
      ++failures;
      continue;
      }

    if( !ridgeOp->DeleteTube( xTube ) )
      {
      std::cout << "Second delete tube failed" << std::endl;
      ++failures;
      continue;
      }

    if( xTube.IsNull() )
      {
      std::cout << "Ridge extraction failed" << std::endl;
      ++failures;
      continue;
      }

    ridgeOp->SmoothTube( xTube, 5 );

    if( !ridgeOp->AddTube( xTube ) )
      {
      std::cout << "Add tube failed" << std::endl;
      ++failures;
      continue;
      }

    if( !ridgeOp->DeleteTube( xTube ) )
      {
      std::cout << "Third delete tube failed" << std::endl;
      ++failures;
      continue;
      }
    }

  RidgeOpType::TubeMaskImageType::Pointer mask =
    ridgeOp->GetTubeMaskImage();
  itk::ImageRegionIterator< RidgeOpType::TubeMaskImageType > maskIt( mask,
    mask->GetLargestPossibleRegion() );
  while( !maskIt.IsAtEnd() )
    {
    if( maskIt.Get() != 0 )
      {
      std::cout << "Final mask not blank." << std::endl;
      std::cout << "Number of failures = " << failures << std::endl;
      return EXIT_FAILURE;
      }
    ++maskIt;
    }

  std::cout << "Number of failures = " << failures << std::endl;
  if( failures > 0 )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
  }
