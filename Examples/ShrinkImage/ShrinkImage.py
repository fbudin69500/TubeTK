#!/usr/bin/env python

# Copyright 2014 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

""" Shrink an image

Command line arguments (See command line help: -h):
---------------------------------------------------
    Required:
        --inputImage (string): input image.
        --outputImage (string): output image.
        --factor (integer): factor by which the image is shrunk
"""


import itk
import argparse
import sys

def ShrinkImage(inputFileName,outputFileName,factor):
  """Loads input image.
     Shrinks image.
     Writes output in given file name.
  """
  reader=itk.ImageFileReader.New(FileName=inputFileName)
  filter=itk.TubeTKITK.ShrinkWithBlendingImage.New(reader)
  filter.SetShrinkFactors([factor,factor,factor])
  writer=itk.ImageFileWriter.New(filter,FileName=outputFileName)
  writer.Update()


def main(argv=None):
    """Parsing command line arguments and reading input files."""
    if argv is None:
        argv = sys.argv
    parser = argparse.ArgumentParser(
            prog=argv[0],
            description=__doc__
    )
    parser.add_argument('-i', "--inputFileName", type=str, required=True, help="Input image")
    parser.add_argument('-o', "--outputFileName", type=str, required=True, help="Output image")
    parser.add_argument('-f', "--factor", type=int, required=True, help="Shrink factor")
    args = parser.parse_args(argv[1:])
    ShrinkImage(args.inputFileName,args.outputFileName,args.factor)


if __name__ == "__main__":
    sys.exit(main())
