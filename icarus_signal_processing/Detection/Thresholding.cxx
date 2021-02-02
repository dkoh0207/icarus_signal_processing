#ifndef __icuars_signal_processing_LOCALTHRESHOLDING_CXX__
#define __icuars_signal_processing_LOCALTHRESHOLDING_CXX__

#include "Thresholding.h"

void icuars_signal_processing::Thresholding::computeIntegralImage(const Array2D<float> &waveform2D,
                                                                  Array2D<float>       &integral2D) const
{
    int numChannels = waveform2D.size();
    int numTicks = waveform2D.at(0).size();

    for (int i = 0; i < numChannels; ++i)
    {
        float sumPixel = 0.0;
        for (int j = 0; j < numTicks; ++j)
        {
            sumPixel = sumPixel + waveform2D[i][j];
            if (i == 0)
            {
                integral2D[i][j] = sumPixel;
            }
            else
            {
                integral2D[i][j] = integral2D[i - 1][j] + sumPixel;
            }
        }
    }
    return;
}

void icuars_signal_processing::Thresholding::Niblack(const Array2D<float> &waveform2D,
                                                     Array2D<bool>        &outputROI,
                                                     const float           k,
                                                     const int             sx,
                                                     const int             sy) const
{
    int numChannels = waveform2D.size();
    int numTicks = waveform2D.at(0).size();

    // Initialize integral image buffer for O(wh) computation of sum/std.
    Array2D<float> sumInt2D;
    sumInt2D.resize(numChannels);
    for (auto &v : sumInt2D)
    {
        v.resize(numTicks);
    }

    Array2D<float> squareInt2D;
    squareInt2D.resize(numChannels);
    for (auto &v : squareInt2D)
    {
        v.resize(numTicks);
    }

    Array2D<float> squareWf2D;
    squareWf2D.resize(numChannels);
    for (auto &v : squareWf2D)
    {
        v.resize(numTicks);
    }

    for (int i = 0; i < numChannels; ++i)
    {
        for (int j = 0; j < numChannels; ++j)
        {
            squareWf2D[i][j] = waveform2D[i][j] * waveform2D[i][j];
        }
    }

    computeIntegralImage(waveform2D, sumInt2D);
    computeIntegralImage(squareWf2D, squareInt2D);

    for (int i = 0; i < numChannels; ++i)
    {
        for (int j = 0; j < numTicks; ++j)
        {
            int lbx = i - (int)sx;
            int ubx = i + (int)sx;
            int lby = j - (int)sy;
            int uby = j + (int)sy;
            int lowerBoundx = std::max(lbx, 0);
            int upperBoundx = std::min(ubx, (int)numChannels);
            int lowerBoundy = std::max(lby, 0);
            int upperBoundy = std::min(uby, (int)numTicks);

            // Compute local mean and standard dev in constant time
            float localMean = (sumInt2D[upperBoundx][upperBoundy] +
                               sumInt2D[lowerBoundx][lowerBoundy] -
                               sumInt2D[lowerBoundx][upperBoundy] -
                               sumInt2D[upperBoundy][lowerBoundy]) /
                              ((float)sx * sy);

            float localMeanSq = (squareInt2D[upperBoundx][upperBoundy] +
                                 squareInt2D[lowerBoundx][lowerBoundy] -
                                 squareInt2D[lowerBoundx][upperBoundy] -
                                 squareInt2D[upperBoundy][lowerBoundy]) /
                                ((float)sx * sy);

            float localStd = localMeanSq - (localMean * localMean);
            float localThreshold = localMean + k * localStd;
            if (waveform2D[i][j] > localThreshold)
            {
                outputROI[i][j] = true;
            }
        }
    }
    return;
}

void icuars_signal_processing::Thresholding::Sauvola(const Array2D<float> &waveform2D,
                                                     Array2D<bool>        &outputROI,
                                                     const float           k,
                                                     const float           R,
                                                     const int             sx,
                                                     const int             sy) const
{
    int numChannels = waveform2D.size();
    int numTicks = waveform2D.at(0).size();

    // Initialize integral image buffer for O(wh) computation of sum/std.
    Array2D<float> sumInt2D;
    sumInt2D.resize(numChannels);
    for (auto &v : sumInt2D)
    {
        v.resize(numTicks);
    }

    Array2D<float> squareInt2D;
    squareInt2D.resize(numChannels);
    for (auto &v : squareInt2D)
    {
        v.resize(numTicks);
    }

    Array2D<float> squareWf2D;
    squareWf2D.resize(numChannels);
    for (auto &v : squareWf2D)
    {
        v.resize(numTicks);
    }

    for (int i = 0; i < numChannels; ++i)
    {
        for (int j = 0; j < numChannels; ++j)
        {
            squareWf2D[i][j] = waveform2D[i][j] * waveform2D[i][j];
        }
    }

    computeIntegralImage(waveform2D, sumInt2D);
    computeIntegralImage(squareWf2D, squareInt2D);

    for (int i = 0; i < numChannels; ++i)
    {
        for (int j = 0; j < numTicks; ++j)
        {
            int lbx = i - (int)sx;
            int ubx = i + (int)sx;
            int lby = j - (int)sy;
            int uby = j + (int)sy;
            int lowerBoundx = std::max(lbx, 0);
            int upperBoundx = std::min(ubx, (int)numChannels);
            int lowerBoundy = std::max(lby, 0);
            int upperBoundy = std::min(uby, (int)numTicks);

            // Compute local mean and standard dev in constant time
            float localMean = (sumInt2D[upperBoundx][upperBoundy] +
                               sumInt2D[lowerBoundx][lowerBoundy] -
                               sumInt2D[lowerBoundx][upperBoundy] -
                               sumInt2D[upperBoundy][lowerBoundy]) /
                              ((float)sx * sy);

            float localMeanSq = (squareInt2D[upperBoundx][upperBoundy] +
                                 squareInt2D[lowerBoundx][lowerBoundy] -
                                 squareInt2D[lowerBoundx][upperBoundy] -
                                 squareInt2D[upperBoundy][lowerBoundy]) /
                                ((float)sx * sy);

            float localStd = localMeanSq - (localMean * localMean);
            float localThreshold = localMean * (1.0 + k * (localStd / R - 1.0));
            if (waveform2D[i][j] > localThreshold)
            {
                outputROI[i][j] = true;
            }
        }
    }
    return;
}

// void icuars_signal_processing::Thresholding::Otsu(
//   const Array2D<float>& waveform2D,
//   Array2D<bool>& outputBinary) const
// {
//   // Run Histogram Normalization
//   std::vector<int> histogram(256);

//   int numChannels = waveform2D.size();
//   int numTicks = waveform2D.at(0).size();

//   float minPixVal = std::numeric_limits<float>::max();
//   float maxPixVal = -std::numeric_limits<float>::max();

//   // Find maximum and minimum pixel intensity
//   for (int i=0; i<numChannels; ++i) {
//     for (int j=0; j<numTicks; ++j) {
//       if (waveform2D[i][j] < minPixVal) {
//         minPixVal = waveform2D[i][j];
//       }
//       if (waveform2D[i][j] > maxPixVal) {
//         maxPixVal = waveform2D[i][j];
//       }
//     }
//   }
//   // Slope for histogram normalization to 0-255
//   float m = 255.0 / ((float) maxPixVal - minPixVal);
//   // Construct Histogram
//   for (int i=0; i<numChannels; ++i) {
//     for (int j=0; j<numTicks; ++j) {
//       float normPix = std::round(m * (waveform2D[i][j] - minPixVal));
//       int normPixInt = (int) normPix;
//       histogram[normPixInt]++;
//     }
//   }
//   return;
// }

void icuars_signal_processing::Thresholding::globalMean(const Array2D<float> &waveform2D,
                                                        Array2D<bool>        &outputBinary,
                                                        float                 k) const
{
    float mean = 0;
    int numChannels = waveform2D.size();
    int numTicks = waveform2D.at(0).size();
    int count = 0;

    // Find maximum and minimum pixel intensity
    for (int i = 0; i < numChannels; ++i)
    {
        for (int j = 0; j < numTicks; ++j)
        {
            mean += waveform2D[i][j];
            count++;
        }
    }
    mean = mean / count;
    for (int i = 0; i < numChannels; ++i)
    {
        for (int j = 0; j < numTicks; ++j)
        {
            if (waveform2D[i][j] > mean * k)
            {
                outputBinary[i][j] = true;
            }
        }
    }
    return;
}

#endif
