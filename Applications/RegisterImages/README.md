TubeTK Shrink Image Application
=============================

---
*This file is part of [TubeTK](http://www.tubetk.org). TubeTK is developed by [Kitware, Inc.](http://www.kitware.com) and licensed under the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0).*

#### Overview:

Provides rigid, affine, and BSpline registration methods via a simple UI

Author(s): Stephen R Aylward (Kitware), Casey B Goodlett (Kitware)

Acknowledgements: This work is part of the National Alliance for Medical
   Image Computing (NAMIC), funded by the National Institutes of Health
   through the NIH Roadmap for Medical Research, Grant U54 EB005149.

#### Commandline Usage:

   RegisterImages   [--returnparameterfile
                    <std::string>]
                    [--processinformationaddress
                    <std::string>] [--xml] [--echo]
                    [--controlPointSpacing <int>]
                    [--bsplineSamplingRatio <float>]
                    [--bsplineMaxIterations <int>]
                    [--affineSamplingRatio <float>]
                    [--affineMaxIterations <int>]
                    [--rigidSamplingRatio <float>]
                    [--rigidMaxIterations <int>]
                    [--movingLandmarks
                    <std::vector<std::vector<float> >>]
                    ...  [--fixedLandmarks
                    <std::vector<std::vector<float> >>]
                    ...  [--interpolation
                    <NearestNeighbor|Linear|BSpline>]
                    [--minimizeMemory]
                    [--numberOfThreads <int>]
                    [--randomNumberSeed <int>]
                    [--fixedImageMask <std::string>]
                    [--sampleFromOverlap]
                    [--verbosityLevel <Silent|Standard
                    |Verbose>] [--expectedSkew <float>]
                    [--expectedScale <float>]
                    [--expectedRotation <float>]
                    [--expectedOffset <float>]
                    [--metric <MattesMI|NormCorr
                    |MeanSqrd>] [--registration <None
                    |Initial|Rigid|Affine|BSpline
                    |PipelineRigid|PipelineAffine
                    |PipelineBSpline>]
                    [--initialization <None|Landmarks
                    |ImageCenters|CentersOfMass
                    |SecondMoments>]
                    [--skipInitialRandomSearch]
                    [--saveDisplacementField
                    <std::string>] [--saveTransform
                    <std::string>] [--loadTransform
                    <std::string>] [--resampledImage
                    <std::string>] [--] [--version]
                    [-h] <std::string> <std::string>


Where:

   --returnparameterfile <std::string>
     Filename in which to write simple return parameters (int, float,
     int-vector, etc.) as opposed to bulk return parameters (image,
     geometry, transform, measurement, table).

   --processinformationaddress <std::string>
     Address of a structure to store process information (progress, abort,
     etc.). (default: 0)

   --xml
     Produce xml description of command line arguments (default: 0)

   --echo
     Echo the command line arguments (default: 0)

   --controlPointSpacing <int>
     Number of pixels between control points (default: 40)

   --bsplineSamplingRatio <float>
     Portion of the image to use in computing the metric during BSpline
     registration (default: 0.1)

   --bsplineMaxIterations <int>
     Maximum number of bspline optimization iterations (default: 20)

   --affineSamplingRatio <float>
     Portion of the image to use in computing the metric during affine
     registration (default: 0.02)

   --affineMaxIterations <int>
     Maximum number of affine optimization iterations (default: 50)

   --rigidSamplingRatio <float>
     Portion of the image to use in computing the metric during rigid
     registration (default: 0.01)

   --rigidMaxIterations <int>
     Maximum number of rigid optimization iterations (default: 100)

   --movingLandmarks <std::vector<std::vector<float> >>  (accepted multiple
      times)
     Ordered list of landmarks in the moving image

   --fixedLandmarks <std::vector<std::vector<float> >>  (accepted multiple
      times)
     Ordered list of landmarks in the fixed image

   --interpolation <NearestNeighbor|Linear|BSpline>
     Method for interpolation within the optimization process (default:
     Linear)

   --minimizeMemory
     Reduce the amount of memory required at the cost of increased
     computation time (default: 0)

   --numberOfThreads <int>
     Number of CPU threads to use (default: 0)

   --randomNumberSeed <int>
     Seed to generate a consistent random number sequence (default: 0)

   --fixedImageMask <std::string>
     Image which defines a mask for the fixed image

   --sampleFromOverlap
     Limit metric evaluation to the fixed image region overlapped by the
     moving image (default: 0)

   --verbosityLevel <Silent|Standard|Verbose>
     Level of detail of reporting progress (default: Standard)

   --expectedSkew <float>
     Expected misalignment after initialization (default: 0.01)

   --expectedScale <float>
     Expected misalignment after initialization (default: 0.05)

   --expectedRotation <float>
     Expected misalignment after initialization (default: 0.1)

   --expectedOffset <float>
     Expected misalignment after initialization (default: 10)

   --metric <MattesMI|NormCorr|MeanSqrd>
     Method to quantify image match (default: MattesMI)

   --registration <None|Initial|Rigid|Affine|BSpline|PipelineRigid
      |PipelineAffine|PipelineBSpline>
     Method for the registration process (default: PipelineAffine)

   --initialization <None|Landmarks|ImageCenters|CentersOfMass
      |SecondMoments>
     Method to prime the registration process (default: CentersOfMass)

   --skipInitialRandomSearch
     Skips initial random search (skips the evolutionary optimizer) during
     the registration process and uses only the gradient optimizer

   --saveDisplacementField <std::string>
     Save displacement field result from registration

   --saveTransform <std::string>
     Save the transform that results from registration

   --loadTransform <std::string>
     Load a transform that is immediately applied to the moving image

   --resampledImage <std::string>
     Registration results

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.

   <std::string>
     (required)  Image which defines the space into which the moving image
     is registered

   <std::string>
     (required)  The transform goes from the fixed image's space into the
     moving image's space
