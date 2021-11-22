#ifndef __icarus_signal_processing_EDGEDETECTION_CXX__
#define __icarus_signal_processing_EDGEDETECTION_CXX__

#include "EdgeDetection.h"
#include "MorphologicalFunctions2D.h"
#include "icarus_signal_processing/Filters/BilateralFilters.h"
#include "icarus_signal_processing/Filters/MiscUtils.h"


long icarus_signal_processing::EdgeDetection::CantorEnum(const int &x, const int &y) const
{
  int n = ((x + y) * (x + y + 1) / 2) + y;
  return n;
}

void icarus_signal_processing::EdgeDetection::Convolve2D(const Array2D<float> &input2D,
                                                         Array2D<float> &output2D,
                                                         const Array2D<float> &kernel) const
{
  // Input kernel must be normalized.

  int numChannels = input2D.size();
  int numTicks = input2D[0].size();

  int kernelX = kernel.size() / 2;
  int kernelY = kernel[0].size() / 2;

  int kernelWidth = kernel.size();
  int kernelHeight = kernel[0].size();

  for (int i = 0; i < numChannels; ++i)
  {
    for (int j = 0; j < numTicks; ++j)
    {

      for (int m = 0; m < kernelWidth; ++m)
      {
        for (int n = 0; n < kernelHeight; ++n)
        {

          int offsetX = m - kernelX;
          int offsetY = n - kernelY;

          int ix = i + offsetX;
          int iy = j + offsetY;

          if (0 <= ix && ix < numChannels && 0 <= iy && iy < numTicks)
          {
            output2D[i][j] += kernel[m][n] * input2D[ix][iy];
          }
        }
      }
    }
  }
  return;
}

void icarus_signal_processing::EdgeDetection::SobelX(const Array2D<float> &input2D,
                                                     Array2D<float> &gradient) const
{
  Array2D<float> kernel(3);
  for (auto &v : kernel)
    v.resize(3);

  // Define SobelX
  kernel[0][0] = 1.0;
  kernel[0][1] = 0.0;
  kernel[0][2] = -1.0;
  kernel[1][0] = 2.0;
  kernel[1][1] = 0.0;
  kernel[1][2] = -2.0;
  kernel[2][0] = 1.0;
  kernel[2][1] = 0.0;
  kernel[2][2] = -1.0;

  Convolve2D(input2D, gradient, kernel);

  return;
}

void icarus_signal_processing::EdgeDetection::SobelY(const Array2D<float> &input2D,
                                                     Array2D<float> &gradient) const
{
  Array2D<float> kernel(3);

  for (auto &v : kernel)
    v.resize(3);

  // Define SobelY
  kernel[0][0] = 1.0;
  kernel[0][1] = 2.0;
  kernel[0][2] = 1.0;
  kernel[1][0] = 0.0;
  kernel[1][1] = 0.0;
  kernel[1][2] = 0.0;
  kernel[2][0] = -1.0;
  kernel[2][1] = -2.0;
  kernel[2][2] = -1.0;

  Convolve2D(input2D, gradient, kernel);
  return;
}

void icarus_signal_processing::EdgeDetection::Sobel(const Array2D<float> &input2D,
                                                    Array2D<float> &sobelX,
                                                    Array2D<float> &sobelY,
                                                    Array2D<float> &gradient,
                                                    Array2D<float> &direction) const
{
  int numChannels = input2D.size();
  int numTicks = input2D[0].size();

  SobelX(input2D, sobelX);
  SobelY(input2D, sobelY);

  for (int i = 0; i < numChannels; ++i)
  {
    for (int j = 0; j < numTicks; ++j)
    {
      float g = sqrt(sobelX[i][j] * sobelX[i][j] + sobelY[i][j] * sobelY[i][j]);
      gradient[i][j] = g;
      float gradDir = atan2(sobelY[i][j], sobelX[i][j]);
      direction[i][j] = gradDir * 180.0 / M_PI;
    }
  }

  return;
}

void icarus_signal_processing::EdgeDetection::SobelRads(const Array2D<float> &input2D,
                                                    Array2D<float> &sobelX,
                                                    Array2D<float> &sobelY,
                                                    Array2D<float> &gradient,
                                                    Array2D<float> &direction) const
{
  int numChannels = input2D.size();
  int numTicks = input2D[0].size();

  SobelX(input2D, sobelX);
  SobelY(input2D, sobelY);

  for (int i = 0; i < numChannels; ++i)
  {
    for (int j = 0; j < numTicks; ++j)
    {
      float g = sqrt(sobelX[i][j] * sobelX[i][j] + sobelY[i][j] * sobelY[i][j]);
      gradient[i][j] = g;
      float gradDir = atan2(sobelY[i][j], sobelX[i][j]);
      direction[i][j] = gradDir;
    }
  }

  return;
}

void icarus_signal_processing::EdgeDetection::SepSobel(const Array2D<float> &input2D,
                                                    Array2D<float> &sobelX,
                                                    Array2D<float> &sobelY,
                                                    Array2D<float> &gradient,
                                                    Array2D<float> &direction) const
{
  int numChannels = input2D.size();
  int numTicks = input2D[0].size();

  SepSobelX(input2D, sobelX);
  SepSobelY(input2D, sobelY);

  for (int i = 0; i < numChannels; ++i)
  {
    for (int j = 0; j < numTicks; ++j)
    {
      float g = sqrt(sobelX[i][j] * sobelX[i][j] + sobelY[i][j] * sobelY[i][j]);
      gradient[i][j] = g;
      float gradDir = atan2(sobelY[i][j], sobelX[i][j]);
      direction[i][j] = gradDir * 180.0 / M_PI;
    }
  }

  return;
}

void icarus_signal_processing::EdgeDetection::SepSobelRads(const Array2D<float> &input2D,
                                                           Array2D<float> &sobelX,
                                                           Array2D<float> &sobelY,
                                                           Array2D<float> &gradient,
                                                           Array2D<float> &direction) const
{
  int numChannels = input2D.size();
  int numTicks = input2D[0].size();

  SepSobelX(input2D, sobelX);
  SepSobelY(input2D, sobelY);

  for (int i = 0; i < numChannels; ++i)
  {
    for (int j = 0; j < numTicks; ++j)
    {
      float g = sqrt(sobelX[i][j] * sobelX[i][j] + sobelY[i][j] * sobelY[i][j]);
      gradient[i][j] = g;
      float gradDir = atan2(sobelY[i][j], sobelX[i][j]);
      direction[i][j] = gradDir;
    }
  }

  return;
}


void icarus_signal_processing::EdgeDetection::SepSobelX(const Array2D<float> &input2D,
                                                        Array2D<float> &gradient) const
{
  const int numChannels = input2D.size();
  const int numTicks = input2D[0].size();

  for (int i = 0; i < numChannels; ++i)
  {
    SepSobelRow(input2D[i], gradient[i]);
  }

  for (int j = 0; j < numTicks; ++j) {
    VectorFloat columnIn(numChannels);
    VectorFloat columnOut(numChannels);
    for (int i = 0; i < numChannels; ++i) {
      columnIn[i] = gradient[i][j];
    }
    SepSobelCol(columnIn, columnOut);
    for (int i = 0; i < numChannels; ++i) {
      gradient[i][j] = columnOut[i];
    }
  }

  return;
}

void icarus_signal_processing::EdgeDetection::SepSobelY(const Array2D<float> &input2D,
                                                        Array2D<float> &gradient) const
{
  const int numChannels = input2D.size();
  const int numTicks = input2D[0].size();

  for (int i = 0; i < numChannels; ++i)
  {
    SepSobelCol(input2D[i], gradient[i]);
  }

  for (int j = 0; j < numTicks; ++j) {
    VectorFloat columnIn(numChannels);
    VectorFloat columnOut(numChannels);
    for (int i = 0; i < numChannels; ++i) {
      columnIn[i] = gradient[i][j];
    }
    SepSobelRow(columnIn, columnOut);
    for (int i = 0; i < numChannels; ++i) {
      gradient[i][j] = columnOut[i];
    }
  }

  return;
}

void icarus_signal_processing::EdgeDetection::SepSobelRow(const VectorFloat &inputRow,
                                                          VectorFloat &outputRow) const
{
  const int N = inputRow.size();

  // Boundary cases
  outputRow[0] = -inputRow[1];
  outputRow[N-1] = inputRow[N-2];

  for (int i=1; i<N-1; ++i) {
    outputRow[i] = inputRow[i-1] - inputRow[i+1];
  }

  return;
}

void icarus_signal_processing::EdgeDetection::SepSobelCol(const VectorFloat &inputRow,
                                                          VectorFloat &outputRow) const
{
  const int N = inputRow.size();

  // Boundary cases
  outputRow[0] = 2.0 * inputRow[0] + inputRow[1];
  outputRow[N-1] = 2.0 * inputRow[N-1] + inputRow[N-2];

  for (int i=1; i<N-1; ++i) {
    outputRow[i] = inputRow[i-1] + inputRow[i+1] + 2.0 * inputRow[i];
  }

  return;
}

void icarus_signal_processing::EdgeDetection::LSDGradX(const Array2D<float> &input2D,
                                                       Array2D<float> &output2D) const
{
  // Input kernel must be normalized.

  int numChannels = input2D.size();
  int numTicks = input2D[0].size();

  for (int i = 0; i < numChannels - 1; ++i)
  {
    for (int j = 0; j < numTicks - 1; ++j)
      output2D[i][j] = (input2D[i + 1][j] + input2D[i + 1][j + 1] -
                        input2D[i][j] - input2D[i][j + 1]) /
                       2.0;
  }

  // Use Periodic Boundary
  for (int j = 0; j < numTicks - 1; ++j)
    output2D[numChannels - 1][j] = (input2D[0][j] + input2D[0][j + 1] -
                                    input2D[numChannels - 1][j] - input2D[numChannels - 1][j + 1]) /
                                   2.0;

  for (int i = 0; i < numChannels - 1; ++i)
    output2D[i][numTicks - 1] = (input2D[i + 1][numTicks - 1] + input2D[i + 1][0] -
                                 input2D[i][numTicks - 1] - input2D[i][0]) /
                                2.0;

  output2D[numChannels - 1][numTicks - 1] = (input2D[0][numTicks - 1] + input2D[0][0] -
                                             input2D[numChannels - 1][numTicks - 1] - input2D[numChannels - 1][0]) /
                                            2.0;
  return;
}

void icarus_signal_processing::EdgeDetection::LSDGradY(const Array2D<float> &input2D,
                                                       Array2D<float> &output2D) const
{
  // Input kernel must be normalized.

  int numChannels = input2D.size();
  int numTicks = input2D[0].size();

  for (int i = 0; i < numChannels - 1; ++i)
  {
    for (int j = 0; j < numTicks - 1; ++j)
      output2D[i][j] = (input2D[i][j + 1] + input2D[i + 1][j + 1] -
                        input2D[i][j] - input2D[i + 1][j]) /
                       2.0;
  }

  // Use Periodic Boundary
  for (int j = 0; j < numTicks - 1; ++j)
    output2D[numChannels - 1][j] = (input2D[numChannels - 1][j + 1] + input2D[0][j + 1] -
                                    input2D[numChannels - 1][j] - input2D[0][j]) /
                                   2.0;

  for (int i = 0; i < numChannels - 1; ++i)
    output2D[i][numTicks - 1] = (input2D[i][0] + input2D[i + 1][0] -
                                 input2D[i][numTicks - 1] - input2D[i + 1][numTicks - 1]) /
                                2.0;

  output2D[numChannels - 1][numTicks - 1] = (input2D[numChannels - 1][0] + input2D[0][0] -
                                             input2D[numChannels - 1][numTicks - 1] - input2D[0][numTicks - 1]) /
                                            2.0;
  return;
}

void icarus_signal_processing::EdgeDetection::LSDGrad(const Array2D<float> &input2D,
                                                      Array2D<float> &gradX,
                                                      Array2D<float> &gradY,
                                                      Array2D<float> &gradient,
                                                      Array2D<float> &direction) const
{
  int numChannels = input2D.size();
  int numTicks = input2D[0].size();

  LSDGradX(input2D, gradX);
  LSDGradY(input2D, gradY);

  for (int i = 0; i < numChannels; ++i)
  {
    for (int j = 0; j < numTicks; ++j)
    {
      float g = sqrt(gradX[i][j] * gradX[i][j] + gradY[i][j] * gradY[i][j]);
      float gradDir = atan2(gradX[i][j], -gradY[i][j]);

      gradient[i][j] = g;
      direction[i][j] = gradDir * 180.0 / M_PI;
    }
  }
  return;
}

void icarus_signal_processing::EdgeDetection::EdgeNMS(const Array2D<float> &gradient2D,
                                                      const Array2D<float> &degrees2D,
                                                      Array2D<float> &output2D) const
{
  int numChannels = gradient2D.size();
  int numTicks = gradient2D[0].size();

  // Compensate for boundary
  for (int i = 1; i < numChannels - 1; ++i)
  {
    for (int j = 1; j < numTicks - 1; ++j)
    {
      if ((112.5 <= degrees2D[i][j] && degrees2D[i][j] < 157.5) || (-67.5 <= degrees2D[i][j] && degrees2D[i][j] < -22.5))
      {
        // Upper Diag Left - Lower Diag Right
        float v1 = gradient2D[i - 1][j + 1];
        float v2 = gradient2D[i + 1][j - 1];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2)
          output2D[i][j] = gradient2D[i][j];
        else
          output2D[i][j] = 0.0;
      }
      else if ((22.5 <= degrees2D[i][j] && degrees2D[i][j] < 67.5) || (-157.5 <= degrees2D[i][j] && degrees2D[i][j] < -112.5))
      {
        // Upper Diag Right - Lower Diag Left
        float v1 = gradient2D[i + 1][j + 1];
        float v2 = gradient2D[i - 1][j - 1];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2)
          output2D[i][j] = gradient2D[i][j];
        else
          output2D[i][j] = 0.0;
      }
      else if ((67.5 <= degrees2D[i][j] && degrees2D[i][j] < 112.5) || (-112.5 <= degrees2D[i][j] && degrees2D[i][j] < -67.5))
      {
        // Up-Down
        float v1 = gradient2D[i][j - 1];
        float v2 = gradient2D[i][j + 1];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2)
          output2D[i][j] = gradient2D[i][j];
        else
          output2D[i][j] = 0.0;
      }
      else
      {
        // Left-Right
        float v1 = gradient2D[i + 1][j];
        float v2 = gradient2D[i - 1][j];
        if (gradient2D[i][j] >= v1 && gradient2D[i][j] >= v2)
          output2D[i][j] = gradient2D[i][j];
        else
          output2D[i][j] = 0.0;
      }
    }
  }

  return;
}

void icarus_signal_processing::EdgeDetection::EdgeNMSInterpolation(const Array2D<float> &gradient2D,
                                                                   const Array2D<float> &gradX,
                                                                   const Array2D<float> &gradY,
                                                                   const Array2D<float> &degrees2D,
                                                                   Array2D<float> &output2D) const
{
    int numChannels = gradient2D.size();
    int numTicks = gradient2D[0].size();

    // Implementation was partly adapted from:
    // https://github.com/JustinLiang/ComputerVisionProjects/blob/master/CannyEdgeDetector/CannyEdgeDetector.m

    // Compensate for boundary
    for (int i = 1; i < numChannels - 1; ++i)
    {
        for (int j = 1; j < numTicks - 1; ++j)
        {
            if ((0 <= degrees2D[i][j] && degrees2D[i][j] < 45) || (-180 <= degrees2D[i][j] && -135 > degrees2D[i][j]))
            {
                float v1 = gradient2D[i + 1][j + 1];
                float v2 = gradient2D[i + 1][j];
                float v3 = gradient2D[i - 1][j];
                float v4 = gradient2D[i - 1][j - 1];
                float G = abs(gradY[i][j] / gradX[i][j]);
                if ((gradient2D[i][j] >= v1 + (v2 - v1) * G) &&
                    (gradient2D[i][j] >= v4 + G * (v3 - v4)))
                {
                    output2D[i][j] = gradient2D[i][j];
                }
                else
                {
                    output2D[i][j] = 0;
                }
            }
            else if ((-45 <= degrees2D[i][j] && degrees2D[i][j] < 0) || (135 <= degrees2D[i][j] && 180 >= degrees2D[i][j]))
            {
                float v1 = gradient2D[i + 1][j - 1];
                float v2 = gradient2D[i + 1][j];
                float v3 = gradient2D[i - 1][j];
                float v4 = gradient2D[i - 1][j + 1];
                float G = abs(gradY[i][j] / gradX[i][j]);
                if ((gradient2D[i][j] >= v1 + (v2 - v1) * G) &&
                    (gradient2D[i][j] >= v4 + G * (v3 - v4)))
                {
                    output2D[i][j] = gradient2D[i][j];
                }
                else
                {
                    output2D[i][j] = 0;
                }
            }
            else if ((45 <= degrees2D[i][j] && degrees2D[i][j] < 90) || (-135 <= degrees2D[i][j] && -90 > degrees2D[i][j]))
            {
                float v1 = gradient2D[i + 1][j + 1];
                float v2 = gradient2D[i][j + 1];
                float v3 = gradient2D[i][j - 1];
                float v4 = gradient2D[i - 1][j - 1];
                float G = abs(gradX[i][j] / gradY[i][j]);
                if ((gradient2D[i][j] >= v1 + (v2 - v1) * G) &&
                    (gradient2D[i][j] >= v4 + G * (v3 - v4)))
                {
                    output2D[i][j] = gradient2D[i][j];
                }
                else
                {
                    output2D[i][j] = 0;
                }
            }
            else
            {
                float v1 = gradient2D[i - 1][j + 1];
                float v2 = gradient2D[i][j + 1];
                float v3 = gradient2D[i][j - 1];
                float v4 = gradient2D[i + 1][j - 1];
                float G = abs(gradX[i][j] / gradY[i][j]);
                if ((gradient2D[i][j] >= v1 + (v2 - v1) * G) &&
                    (gradient2D[i][j] >= v4 + G * (v3 - v4)))
                {
                    output2D[i][j] = gradient2D[i][j];
                }
                else
                {
                    output2D[i][j] = 0;
                }
            }
        }
    }

    return;
}

void icarus_signal_processing::EdgeDetection::DoubleThresholding(const Array2D<float> &doneNMS2D,
                                                                 Array2D<bool>        &binary2D,
                                                                 std::vector<int>     &strongEdgeRows,
                                                                 std::vector<int>     &strongEdgeCols,
                                                                 std::vector<int>     &weakEdgeRows,
                                                                 std::vector<int>     &weakEdgeCols,
                                                                 float                 lowThreshold,
                                                                 float                 highThreshold) const
{
    int numChannels = doneNMS2D.size();
    int numTicks = doneNMS2D[0].size();

    // Implementation was partly adapted from:
    // https://github.com/JustinLiang/ComputerVisionProjects/blob/master/CannyEdgeDetector/CannyEdgeDetector.m

    strongEdgeRows.reserve(numChannels * numTicks);
    strongEdgeCols.reserve(numChannels * numTicks);
    weakEdgeRows.reserve(numChannels * numTicks);
    weakEdgeCols.reserve(numChannels * numTicks);
    // Enumerator for strong / weak edges

    // Compensate for boundary
    for (int i = 0; i < numChannels; ++i)
    {
        for (int j = 0; j < numTicks; ++j)
        {
            if (doneNMS2D[i][j] >= highThreshold)
            {
                binary2D[i][j] = true;
                strongEdgeRows.push_back(i);
                strongEdgeCols.push_back(j);
            }
            else if ((doneNMS2D[i][j] < highThreshold) && (doneNMS2D[i][j] >= lowThreshold))
            {
                weakEdgeRows.push_back(i);
                weakEdgeCols.push_back(j);
            }
            else
            {
                binary2D[i][j] = false;
            }
        }
    }
    return;
}

void icarus_signal_processing::EdgeDetection::HysteresisThresholding(const Array2D<bool>    &binary2D,
                                                                     const std::vector<int> &strongEdgeRows,
                                                                     const std::vector<int> &strongEdgeCols,
                                                                     const std::vector<int> &weakEdgeRows,
                                                                     const std::vector<int> &weakEdgeCols,
                                                                     Array2D<bool>          &output2D) const
{
    int numChannels = binary2D.size();
    int numTicks = binary2D[0].size();

    bool converged = false;

    // Initialize fixed point iteration buffer
    Array2D<bool> tempBuffer(numChannels);
    for (auto &v : tempBuffer) v.resize(numTicks);

    if (strongEdgeRows.size() != strongEdgeCols.size()) throw std::runtime_error("HysteresisThresholding: strong edge rows not same as columnns");
    if (weakEdgeRows.size()   != weakEdgeCols.size())   throw std::runtime_error("HysteresisThresholding: weak edge rows not same as columnns");

    for (size_t i = 0; i < strongEdgeRows.size(); ++i)
    {
        tempBuffer[strongEdgeRows[i]][strongEdgeCols[i]] = true;
        output2D[strongEdgeRows[i]][strongEdgeCols[i]]   = true;
    }

    // Construct strong edge / weak edge 2D binary arrays from indices
    // REMARK: binary2D already contains strong edges.
    // Array2D<bool> strongEdges2D(numChannels);
    // for (auto& v : strongEdge2D) {
    //   v.resize(numTicks);
    // }
    Array2D<bool> weakEdges2D(numChannels);
    for (auto &v : weakEdges2D)
    {
        v.resize(numTicks);
    }

    for (size_t i = 0; i < weakEdgeRows.size(); ++i)
    {
        weakEdges2D[weakEdgeRows[i]][weakEdgeCols[i]] = true;
    }

    // This fixed point iteration always converges.
    while (!converged)
    {
        // Dilation + compute overlap
        getDilation2D(output2D, 3, 3, tempBuffer);

        for (int i = 0; i < numChannels; ++i)
        {
            for (int j = 0; j < numTicks; ++j)
            {
                // Pixels have to be true in both dilated binary image and weak edges
                tempBuffer[i][j] = output2D[i][j] || (tempBuffer[i][j] && weakEdges2D[i][j]);
            }
        }
        converged = true;
        // Check array equality
        for (int i = 0; i < numChannels; ++i)
        {
            for (int j = 0; j < numTicks; ++j)
            {
                if (tempBuffer[i][j] != output2D[i][j])
                {
                    converged = false;
                    break;
                }
            }

            if (!converged) break;
        }

        // Update buffer
        for (int i = 0; i < numChannels; ++i)
        {
            for (int j = 0; j < numTicks; ++j)
            {
                output2D[i][j] = tempBuffer[i][j];
            }
        }
    }
    return;
}

void icarus_signal_processing::EdgeDetection::HysteresisThresholdingFast(const Array2D<float> &doneNMS2D,
                                                                         const float lowThreshold,
                                                                         const float highThreshold,
                                                                         Array2D<bool> &outputROI) const
{
    /*
  Hysteresis Thresholding using Disjoint Set Forest Data Structure and Union-Find

  Reference:
  Artur Nowakowski, Wladyslaw Skarbek, "Fast computation of thresholding hysteresis for edge detection," 
  Proc. SPIE 6159, Photonics Applications in Astronomy, Communications, Industry, and 
  High-Energy Physics Experiments IV, 615948 (26 April 2006); https://doi.org/10.1117/12.674959

  This Hysteresis Thresholding includes double thresholding, and performs both in one pass. 
  */
    const int numChannels = doneNMS2D.size();
    const int numTicks = doneNMS2D[0].size();

    const int forestSize = numChannels * numTicks;

    DisjointSetForest forest(forestSize, forestSize);
    forest.MakeSet();

    // 1. Initialize Strong Edges

    for (int i = 0; i < numChannels; ++i)
    {

        for (int j = 0; j < numTicks; ++j)
        {

            int flatIndex = i * numTicks + j;

            if (doneNMS2D[i][j] >= highThreshold)
                forest.parent[flatIndex] = forestSize;
        }
    }

    for (int i = 0; i < numChannels; ++i)
    {

        for (int j = 0; j < numTicks; ++j)
        {

            int flatIndex = i * numTicks + j;
            int lowerBoundx = std::max(i - 1, 0);
            int upperBoundx = std::min(i + 2, (int)numChannels);
            int lowerBoundy = std::max(j - 1, 0);
            int upperBoundy = std::min(j + 2, (int)numTicks);

            // Process strong edges and its neighbors
            if (doneNMS2D[i][j] >= highThreshold)
            {
                // Assign every strong edge to single root node (no need for unions)
                // Root node has index forestSize (last element of array with size forestSize + 1)
                if (forest.Find(flatIndex) == flatIndex)
                {

                    forest.parent[flatIndex] = forestSize;
                }
                // Handle neighboring weak edges

                for (int k = lowerBoundx; k < upperBoundx; ++k)
                {

                    for (int l = lowerBoundy; l < upperBoundy; ++l)
                    {

                        int flatIndexNeigh = k * numTicks + l;
                        const float &grad = doneNMS2D[k][l];

                        if (grad >= lowThreshold)
                        {

                            forest.Union(flatIndexNeigh, flatIndex);
                        }
                    }
                }
            }
            // Process weak edges
            else if ((doneNMS2D[i][j] < highThreshold) && (doneNMS2D[i][j] >= lowThreshold))
            {

                for (int k = lowerBoundx; k < upperBoundx; ++k)
                {

                    for (int l = lowerBoundy; l < upperBoundy; ++l)
                    {

                        int flatIndexNeigh = k * numTicks + l;
                        const float &grad = doneNMS2D[k][l];

                        if (grad >= lowThreshold)
                        {

                            forest.Union(flatIndexNeigh, flatIndex);
                        }
                    }
                }
            }
            else
                continue;
        }
    }

    for (int flatIdx = 0; flatIdx < forestSize; ++flatIdx)
    {

        int rep = forest.Find(flatIdx);
        int row = (flatIdx / numTicks);
        int col = (flatIdx % numTicks);

        if (rep == forestSize)
            outputROI[row][col] = true;

        else
            outputROI[row][col] = false;
    }
    return;
}

void icarus_signal_processing::EdgeDetection::SparseHysteresisThresholding(const Array2D<float> &doneNMS2D,
                                                                           const float lowThreshold,
                                                                           const float highThreshold,
                                                                           Array2D<bool> &outputROI) const
{
    /*
  Hysteresis Thresholding using Disjoint Set Forest Data Structure and Union-Find

  Reference:
  Artur Nowakowski, Wladyslaw Skarbek, "Fast computation of thresholding hysteresis for edge detection," 
  Proc. SPIE 6159, Photonics Applications in Astronomy, Communications, Industry, and 
  High-Energy Physics Experiments IV, 615948 (26 April 2006); https://doi.org/10.1117/12.674959

  This Hysteresis Thresholding includes double thresholding, and performs both in one pass. 
  */
    const int numChannels = doneNMS2D.size();
    const int numTicks = doneNMS2D.at(0).size();

    std::unordered_map<long, EdgeCandidate> edges;

    // 1. Initialize Strong Edges

    int count_id = 1;
    // First Pass
    for (int i=0; i<numChannels; ++i) 
    {
        for (int j=0; j<numTicks; ++j)
        {
            if (doneNMS2D[i][j] >= highThreshold) 
            {
              // if edge is a strong edge, assign to root node 
              // ID of 0 is reserved for root note reference
              EdgeCandidate strongEdge(i, j, 0, true);
              long key = CantorEnum(i, j);
              edges.emplace(std::make_pair(key, strongEdge));
              count_id++;
            }
            else if ( (doneNMS2D[i][j] < highThreshold) && 
                      (doneNMS2D[i][j] >= lowThreshold)) 
            {
              EdgeCandidate weakEdge(i, j, count_id, false);
              long key = CantorEnum(i, j);
              edges.emplace(std::make_pair(key, weakEdge));
              count_id++;
            }
            else continue;
        }
    }

    const int forestSize = edges.size();
    DisjointSetForest forest(forestSize);
    forest.MakeSet();

    // Assign all strong edge to root node.
    for (auto& node : edges) 
    {
        EdgeCandidate &edge = node.second;
        int i = edge.row;
        int j = edge.col;

        int lowerBoundx = std::max(i-1, 0);
        int upperBoundx = std::min(i+2, (int) numChannels);
        int lowerBoundy = std::max(j-1, 0);
        int upperBoundy = std::min(j+2, (int) numTicks);

        if (edge.edgeType)
        {
            if (forest.Find(edge.id) == edge.id) 
            {
              // Assign strong edge to root node "0"
              forest.parent[edge.id] = 0;
            }
             // Handle neighbors
            for (int k=lowerBoundx; k<upperBoundx; ++k) 
            {
                for (int l=lowerBoundy; l<upperBoundy; ++l) 
                {
                    const float &grad = doneNMS2D[k][l];
                    if (grad >= lowThreshold) 
                    {
                            long key = CantorEnum(k, l);
                            forest.Union(edges[key].id, edge.id);
                    }
                }
            }
        }
        else // Process Weak Edges 
        {
            for (int k=lowerBoundx; k<upperBoundx; ++k) 
            {
                for (int l=lowerBoundy; l<upperBoundy; ++l) 
                {
                    const float &grad = doneNMS2D[k][l];
                    if (grad >= lowThreshold) 
                    {
                      long key = CantorEnum(k, l);
                      forest.Union(edges[key].id, edge.id);
                    }
                }
            }
        }
    }

    for (auto& node: edges) 
    {
        EdgeCandidate &edge = node.second;
        int rep = forest.Find(edge.id);
        const int &row = edge.row;
        const int &col = edge.col;
        if (rep == 0) outputROI[row][col] = true;
    }
    return;
}

void icarus_signal_processing::EdgeDetection::Canny(const Array2D<float> &waveLessCoherent,
                                                    Array2D<bool>        &output2D,
                                                    const unsigned int    sx,
                                                    const unsigned int    sy,
                                                    const float           sigma_x,
                                                    const float           sigma_y,
                                                    const float           sigma_r,
                                                    const float           lowThreshold,
                                                    const float           highThreshold,
                                                    const char            mode) const
{
    /*
  This implementation of canny edge detection replaces gaussian smoothing with
  edge-preserving bilateral filtering.
  */
    int numChannels = waveLessCoherent.size();
    int numTicks = waveLessCoherent[0].size();

    icarus_signal_processing::BilateralFilters filter;

    Array2D<float> temp(numChannels);
    for (auto &v : temp)
    {
        v.resize(numTicks);
    }

    Array2D<bool> boolTemp(numChannels);
    for (auto &v : boolTemp)
    {
        v.resize(numTicks);
    }

    Array2D<float> gradient(numChannels);
    for (auto &v : gradient)
    {
        v.resize(numTicks);
    }

    Array2D<float> direction(numChannels);
    for (auto &v : direction)
    {
        v.resize(numTicks);
    }

    Array2D<float> sobelX(numChannels);
    for (auto &v : sobelX)
    {
        v.resize(numTicks);
    }

    Array2D<float> sobelY(numChannels);
    for (auto &v : sobelY)
    {
        v.resize(numTicks);
    }

    Array2D<float> morphed2D(numChannels);
    for (auto &v : morphed2D)
    {
        v.resize(numTicks);
    }

    SepSobelRads(waveLessCoherent, sobelX, sobelY, gradient, direction);

    // 1. Perform Edge-Preserving Smoothing
    // Here sx = 7, sy = 7, sigma_x = 10, sigma_y = 10, sigma_r = 30
    filter.directional(waveLessCoherent, direction, temp, sx, sy, sigma_x, sigma_y, sigma_r, 360);

    // 2. Run Morphological FIlter
    if (mode == 'e')
    {
        icarus_signal_processing::Erosion2D erosion2D(sx,sy);
        erosion2D(temp.begin(), temp.size(), morphed2D.begin());
    }
    else if (mode == 'd')
    {
        icarus_signal_processing::Dilation2D dilation2D(sx,sy);
        dilation2D(temp.begin(), temp.size(), morphed2D.begin());
    }
    else
    {
        icarus_signal_processing::Gradient2D gradient2D(sx,sy);
        gradient2D(temp.begin(), temp.size(), morphed2D.begin());
    }

    // 2. Run Sobel Edge Filtering
    SepSobel(morphed2D, sobelX, sobelY, gradient, direction);

    // 3. NMS Edge with Interpolation
    EdgeNMSInterpolation(gradient, sobelX, sobelY, direction, temp);

    // 4. Run double thresholding

    // std::vector<int> strongEdgeRows;
    // std::vector<int> strongEdgeCols;
    // std::vector<int> weakEdgeRows;
    // std::vector<int> weakEdgeCols;

    HysteresisThresholdingFast(temp, lowThreshold, highThreshold, output2D);

    // DoubleThresholding(temp, boolTemp,
    //                    strongEdgeRows, strongEdgeCols,
    //                    weakEdgeRows, weakEdgeCols, lowThreshold, highThreshold);

    // HysteresisThresholding(boolTemp, strongEdgeRows, strongEdgeCols, weakEdgeRows, weakEdgeCols, output2D);

    return;
}

void icarus_signal_processing::EdgeDetection::gradientRegionGrow(const Array2D<float> &direction,
                                                                 const int             anchorX,
                                                                 const int             anchorY,
                                                                 const int             regionID,
                                                                 const float           tolerance,
                                                                 const unsigned int    windowX,
                                                                 const unsigned int    windowY,
                                                                 Array2D<int>         &partitions) const
{
    assert(tolerance > 0);
    assert(regionID > 0);
    int numChannels = direction.size();
    int numTicks = direction[0].size();

    float Sx = 0.0;
    float Sy = 0.0;
    float meanAngle = direction[anchorX][anchorY];
    std::queue<int> anchorsX;
    std::queue<int> anchorsY;

    anchorsX.push(anchorX);
    anchorsY.push(anchorY);

    unsigned int count = 1;

    while (!anchorsX.empty())
    {
        // std::cout << anchorsX.size() << std::endl;
        // std::cout << anchorsY.size() << std::endl;
        // Neighbors of current pixel
        int ix = anchorsX.front();
        int iy = anchorsY.front();
        anchorsX.pop();
        anchorsY.pop();

        int lbx = ix - (int)windowX;
        int ubx = ix + (int)windowX;
        int lby = iy - (int)windowY;
        int uby = iy + (int)windowY;
        int lowerBoundx = std::max(lbx, 0);
        int upperBoundx = std::min(ubx, (int)numChannels);
        int lowerBoundy = std::max(lby, 0);
        int upperBoundy = std::min(uby, (int)numTicks);

        for (int n = lowerBoundx; n < upperBoundx; ++n)
        {
            for (int m = lowerBoundy; m < upperBoundy; ++m)
            {
                float diffAngle = (float)abs(meanAngle - direction[n][m]);
                // std::cout << "Mean Angle = " << meanAngle << std::endl;
                // std::cout << "Diff Angle = " << diffAngle << std::endl;
                if (partitions[n][m] == 0 && diffAngle < tolerance)
                {
                    anchorsX.push(n);
                    anchorsY.push(m);
                    count++;
                    partitions[n][m] = regionID;
                    Sx = Sx + cos(direction[n][m] * M_PI / 180.0);
                    Sy = Sy + sin(direction[n][m] * M_PI / 180.0);
                    meanAngle = atan2(Sy, Sx) * 180.0 / M_PI;
                }
            }
        }
    }
    return;
}

void icarus_signal_processing::EdgeDetection::regionGrow2D(const Array2D<float>   &direction,
                                                           const std::vector<int> &anchorsX,
                                                           const std::vector<int> &anchorsY,
                                                           const float             tolerance,
                                                           const unsigned int      windowX,
                                                           const unsigned int      windowY,
                                                           Array2D<int>           &partitions) const
{
    assert(anchorsX.size() == anchorsY.size());

    int groupID = 1;
    for (size_t i = 0; i < anchorsX.size(); ++i)
    {
        // std::cout << "Label = " << groupID << std::endl;
        int ix = anchorsX[i];
        int iy = anchorsY[i];
        gradientRegionGrow(direction, ix, iy, groupID, tolerance, windowX, windowY, partitions);
        groupID++;
    }
    return;
}

void icarus_signal_processing::EdgeDetection::getDilation1D(const std::vector<bool>& inputWaveform,
                                                            const unsigned int       structuringElement,
                                                            std::vector<bool>&       dilationVec) const
{
    const size_t N = inputWaveform.size();
    const size_t k = (size_t)structuringElement;

    assert(dilationVec.size() == N);

    if (N <= k)
    {
        std::cout << "Input array size " << N << " must be greater than structuring element size " << k << std::endl;
        return;
    }
    const size_t windowSize = k/2;
    const size_t paddingSize = (k - (N % k)) % k;
    const size_t bufferSize = N + 2 * windowSize + paddingSize;
    std::vector<bool> suffixArr(bufferSize);
    std::vector<bool> prefixArr(bufferSize);

    // Padding Operations on Buffers
    for (size_t i = 0; i < windowSize; ++i)
    {
        suffixArr[i] = false;
        prefixArr[i] = false;
    }

    for (size_t i = N + windowSize; i < bufferSize; ++i)
    {
        suffixArr[i] = false;
        prefixArr[i] = false;
    }

    // Compute Prefix and Suffix Buffers
    for (size_t i = 0; i < N + paddingSize; ++i)
    {
        if (i % k == 0)
        {
            prefixArr[i + windowSize] = inputWaveform[i];
        }
        else if ((i % k == 0) && (i < N))
        {
            prefixArr[i + windowSize] = (prefixArr[i + windowSize - 1] || inputWaveform[i]);
        }
        else continue;
    }

    for (size_t i = N + paddingSize; i != 0; --i)
    {
        if (i > N)
        {
            // Compensate for divisibility padding (must be -inf)
            continue;
        }
        else if (i % k == 0)
        {
            suffixArr[i + windowSize - 1] = inputWaveform[i - 1];
        }
        else
        {
            suffixArr[i + windowSize - 1] = (suffixArr[i + windowSize] || inputWaveform[i - 1]);
        }
    }

    int prefixIndex = 0;
    int suffixIndex = 0;

    for (size_t i = windowSize; i < N + windowSize; ++i)
    {
        prefixIndex = i + windowSize;
        suffixIndex = i - windowSize;
        dilationVec[i - windowSize] = (prefixArr[prefixIndex] || suffixArr[suffixIndex]);
    }
    return;
}


void icarus_signal_processing::EdgeDetection::getDilation2D(const std::vector<std::vector<bool> >& waveform2D,
                                                            const unsigned int                     structuringElementx,
                                                            const unsigned int                     structuringElementy,
                                                            std::vector<std::vector<bool> >&       dilation2D) const
{
    size_t numChannels = waveform2D.size();
    size_t nTicks = waveform2D[0].size();

    assert(dilation2D.size() == numChannels);
    assert(dilation2D[0].size() == nTicks);

    for (size_t i = 0; i < numChannels; ++i)
    {
        getDilation1D(waveform2D[i], structuringElementy, dilation2D[i]);
    }
    for (size_t j = 0; j < nTicks; ++j)
    {
        std::vector<bool> column(numChannels);
        std::vector<bool> columnOut(numChannels);
        for (size_t i = 0; i < numChannels; ++i)
        {
            column[i] = dilation2D[i][j];
        }

        getDilation1D(column, structuringElementx, columnOut);

        for (size_t i = 0; i < numChannels; ++i)
        {
            dilation2D[i][j] = columnOut[i];
        }
    }
    return;
}


#endif
