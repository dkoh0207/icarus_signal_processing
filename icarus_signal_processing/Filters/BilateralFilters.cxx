#ifndef __icarus_signal_processing_BILATERALFILTERS_CXX__
#define __icarus_signal_processing_BILATERALFILTERS_CXX__

#include "BilateralFilters.h"
#include "MiscUtils.h"


void icarus_signal_processing::BilateralFilters::bilateral(
  const Array2D<float>& input2D,
  Array2D<float>& output2D,
  const unsigned int sx,
  const unsigned int sy,
  float sigma_x,
  float sigma_y,
  float sigma_r) const
{
  /*
  Brute force implementation of Bilateral Filtering
  WARNING: O(N^2) complexity, algorithm may be slow (~17s for 1050 x 4098)
  */
  int numChannels = input2D.size();
  int numTicks = input2D[0].size();
  int xHalfWindowSize = sx / 2;
  int yHalfWindowSize = sy / 2;

  // Precompute weights
  std::vector<float> weightX(xHalfWindowSize+1);
  std::vector<float> weightY(yHalfWindowSize+1);

  weightX[0] = 1.0;
  weightY[0] = 1.0;

  for (int n=1; n<=xHalfWindowSize; ++n) {
    weightX[n] = exp(-n*n / (2.0 * sigma_x * sigma_x));
  }

  for (int n=1; n<=yHalfWindowSize; ++n) {
    weightY[n] = exp(-n*n / (2.0 * sigma_y * sigma_y));
  }

  output2D.resize(numChannels);
  for (auto& v : output2D) {
    v.resize(numTicks);
  }

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {

      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      int lowerBoundx = std::max(lbx, 0);
      int upperBoundx = std::min(ubx, (int) numChannels);
      int lowerBoundy = std::max(lby, 0);
      int upperBoundy = std::min(uby, (int) numTicks);

      float newPixel = 0.0;
      float normFactor = 0.0;

      for (int ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (int iy=lowerBoundy; iy<upperBoundy; ++iy) {
          float w = 1.0;
          w = weightX[abs(ix - i)] * weightY[abs(iy - j)];
          w *= exp(-(input2D[i][j] - input2D[ix][iy])*(input2D[i][j] - input2D[ix][iy]) / (2.0 * sigma_r * sigma_r));
          newPixel += w * input2D[ix][iy];
          normFactor += w;
        }
      }

      newPixel = newPixel / normFactor;
      output2D[i][j] = newPixel;
    }
  }
  return;
}


void icarus_signal_processing::BilateralFilters::directional(
  const Array2D<float>& input2D,
  const Array2D<float>& direction,
  Array2D<float>& output2D,
  const unsigned int sx,
  const unsigned int sy,
  float sigma_x,
  float sigma_y,
  float sigma_theta,
  const int angleMesh = 360) const
{
  /*
  Brute force implementation of Bilateral Filtering
  WARNING: O(sx * sy) complexity, algorithm may be slow (~17s for 1050 x 4098)
  */
  int numChannels = input2D.size();
  int numTicks = input2D[0].size();
  int xHalfWindowSize = sx / 2;
  int yHalfWindowSize = sy / 2;

  // Precompute weights
  std::vector<float> weightX(xHalfWindowSize+1);
  std::vector<float> weightY(yHalfWindowSize+1);
  std::vector<float> weightAngle(angleMesh+2);

  weightX[0] = 1.0;
  weightY[0] = 1.0;
  weightAngle[0] = 1.0;

  for (int n=1; n<=xHalfWindowSize; ++n) {
    weightX[n] = exp(-n*n / (2.0 * sigma_x * sigma_x));
  }

  for (int n=1; n<=yHalfWindowSize; ++n) {
    weightY[n] = exp(-n*n / (2.0 * sigma_y * sigma_y));
  }

  for (int n=1; n<=angleMesh+1; ++n) {
    weightAngle[n] = exp(-n*n / (2.0 * sigma_theta * sigma_theta));
  }

  output2D.resize(numChannels);
  for (auto& v : output2D) {
    v.resize(numTicks);
  }

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {

      int lbx = i - (int) xHalfWindowSize;
      int ubx = i + (int) xHalfWindowSize;
      int lby = j - (int) yHalfWindowSize;
      int uby = j + (int) yHalfWindowSize;
      int lowerBoundx = std::max(lbx, 0);
      int upperBoundx = std::min(ubx, (int) numChannels);
      int lowerBoundy = std::max(lby, 0);
      int upperBoundy = std::min(uby, (int) numTicks);

      float newPixel = 0.0;
      float normFactor = 0.0;

      for (int ix=lowerBoundx; ix<upperBoundx; ++ix) {
        for (int iy=lowerBoundy; iy<upperBoundy; ++iy) {
          float w = 1.0;
          w = weightX[abs(ix - i)] * weightY[abs(iy - j)];
          // w *= exp(-(direction[i][j] - direction[ix][iy])*(direction[i][j] - direction[ix][iy]) 
          //    / (2.0 * sigma_theta * sigma_theta));
          int dtheta = (int) abs(direction[i][j] - direction[ix][iy]);
          w *= weightAngle[dtheta];
          newPixel += w * input2D[ix][iy];
          normFactor += w;
        }
      }

      newPixel = newPixel / normFactor;
      output2D[i][j] = newPixel;
    }
  }
  return;
}


void icarus_signal_processing::BilateralFilters::Gaussian(
  const Array2D<float>& input2D,
  const float sigma_x,
  const float sigma_y,
  Array2D<float>& output2D) const
{
  Gaussian<float>(input2D, sigma_x, sigma_y, output2D);
  return;
}

void icarus_signal_processing::BilateralFilters::Gaussian(
  const Array2D<std::complex<float>>& input2D,
  const float sigma_x,
  const float sigma_y,
  Array2D<std::complex<float>>& output2D) const
{
  Gaussian<std::complex<float>>(input2D, sigma_x, sigma_y, output2D);
  return;
}

template <typename T>
void icarus_signal_processing::BilateralFilters::Gaussian(
  const Array2D<T>& input2D,
  const float sigma_x,
  const float sigma_y,
  Array2D<T>& output2D) const
{
  /*
  The constant time O(1) recursive gaussian filter
  See Young and Vliet, Recursive Implementation of the Gaussian Filter
  */
  if (sigma_x < 0.5 || sigma_y < 0.5) throw std::runtime_error("BilateralFilters: sigma_x and/or sigma_y < 0.5");

  size_t numChannels = input2D.size();
  size_t numTicks = input2D[0].size();

  if (output2D.size() != numChannels) throw std::runtime_error("CSUBilateral: output array does not match numChannels");
  if (output2D[0].size() != numTicks) throw std::runtime_error("CSUBilateral: output length does not match numTicks");

  // Initialize Recursive Gaussian Parameters
  std::vector<float> filterCoeffsX = computeGausCoeffs(sigma_x);
  std::vector<float> filterCoeffsY = computeGausCoeffs(sigma_y);
  
  // Gaussian filter is separable, so 1D gaussian applied row/column in seq.
  for (size_t i=0; i<numChannels; ++i) {
    Gaussian1D(input2D[i], filterCoeffsY, output2D[i]);
  }
  for (size_t j=0; j<numTicks; ++j) {
    std::vector<T> columnIn(numChannels);
    std::vector<T> columnOut(numChannels);
    for (size_t i=0; i<numChannels; ++i) {
      columnIn[i] = output2D[i][j];
    }
    Gaussian1D(columnIn, filterCoeffsX, columnOut);
    for (size_t i=0; i<numChannels; ++i) {
      output2D[i][j] = columnOut[i];
    }
  }

  return;
}


std::vector<float> icarus_signal_processing::BilateralFilters::computeGausCoeffs(
  const float sigma) const
{
  // Compute q-value
  float q = 0.0;
  if (sigma >= 2.5) {
    q = 0.98711 * sigma  - 0.96330;
  } 
  else {
    q = 3.97156-4.14554 * sqrt(1 - 0.26891 * sigma);
  }
  // Compute linear filter coefficients
  float b0 = 1.57825 + 2.44413*q + 1.4281*std::pow(q, 2.0) + 0.422205*std::pow(q, 3.0);
  float b1 = (2.44413*q) + (2.85619*std::pow(q, 2.0)) + (1.26661*std::pow(q, 3.0));
  float b2 = -(1.4281*std::pow(q, 2.0)) - (1.26661*std::pow(q, 3.0));
  float b3 = 0.422205*std::pow(q, 3.0);
  float B = 1 - (b1 + b2 + b3) / b0;
  std::vector<float> filterCoeffs = {b0, b1, b2, b3, B};
  return filterCoeffs;
}


void icarus_signal_processing::BilateralFilters::Gaussian1D(
  const std::vector<float>& input1D,
  const std::vector<float>& filterCoeffs,
  std::vector<float>& output1D) const
{
  Gaussian1D<float>(input1D, filterCoeffs, output1D);
  return;
}

void icarus_signal_processing::BilateralFilters::Gaussian1D(
  const std::vector<std::complex<float>>& input1D,
  const std::vector<float>& filterCoeffs,
  std::vector<std::complex<float>>& output1D) const
{
  Gaussian1D<std::complex<float>>(input1D, filterCoeffs, output1D);
  return;
}

template <typename T>
void icarus_signal_processing::BilateralFilters::Gaussian1D(
  const std::vector<T>& input1D,
  const std::vector<float>& filterCoeffs,
  std::vector<T>& output1D) const
{
  if (filterCoeffs.size() != 5)          throw std::runtime_error("Gaussian1D: filterCoeffs not of size 5");
  if (output1D.size() != input1D.size()) throw std::runtime_error("Gaussian1D: input not same size as output");

  size_t N = input1D.size();

  float b0 = filterCoeffs[0];
  float b1 = filterCoeffs[1];
  float b2 = filterCoeffs[2];
  float b3 = filterCoeffs[3];
  float B = filterCoeffs[4];

  std::vector<T> w(N+6); // Padding left 3 + right 3
  w[0] = input1D[0];
  w[1] = input1D[0];
  w[2] = input1D[0];
  w[N] = input1D[0];
  w[N+1] = input1D[0];
  w[N+2] = input1D[0];

  std::vector<T> c(N+6); // Padding left 3 + right 3
  c[0] = input1D[0];
  c[1] = input1D[0];
  c[2] = input1D[0];
  c[N] = input1D[0];
  c[N+1] = input1D[0];
  c[N+2] = input1D[0];

  // Do forward pass
  for (size_t i=3; i<N+3; ++i) {
    w[i] = B*input1D[i] + (b1*w[i-1] + b2*w[i-2] + b3*w[i-3]) / b0;
  }
  // Do backward pass
  for (size_t i=N+2; i>2; --i) {
    c[i] = B*w[i] + (b1*c[i+1] + b2*c[i+2] + b3*c[i+3]) / b0;
  }
  
  for (size_t i=0; i<N; ++i) {
    output1D[i] = c[i+3];
  }
  return;
}


float icarus_signal_processing::BilateralFilters::Discretize(
  const Array2D<float>& input2D,
  const unsigned int maxPixelValue,
  Array2D<short>& output2D) const
{
  float m = Discretize<short>(input2D, maxPixelValue, output2D);
  return m;
}

float icarus_signal_processing::BilateralFilters::Discretize(
  const Array2D<float>& input2D,
  const unsigned int maxPixelValue,
  Array2D<float>& output2D) const
{
  float m = Discretize<float>(input2D, maxPixelValue, output2D);
  return m;
}

template <typename T>
float icarus_signal_processing::BilateralFilters::Discretize(
  const Array2D<float>& input2D,
  const unsigned int maxPixelValue,
  Array2D<T>& output2D) const
{

  size_t numChannels = input2D.size();
  size_t numTicks = input2D[0].size();

  if (output2D.size() != numChannels) throw std::runtime_error("Discretize: output size does not match numChannels");
  if (output2D[0].size() != numTicks) throw std::runtime_error("Discretize: output length does not match numTicks");
  if (maxPixelValue == 0)             throw std::runtime_error("Discretize: maxPixelVale is zero");

  icarus_signal_processing::MiscUtils utils;
  float lb = utils.computeMinimum(input2D);
  float ub = utils.computeMaximum(input2D);

  float m = ((float) maxPixelValue) / (ub - lb);

  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<numTicks; ++j) {
      T newVal = (T) std::round(m * (input2D[i][j] - lb));
      output2D[i][j] = newVal;
    }
  }
  return m;
}


void icarus_signal_processing::BilateralFilters::CSUBilateral(
  const Array2D<float>& input2D,
  Array2D<float>& output2D,
  const unsigned int T,
  const float sigma_s,
  const float sigma_r) const
{

  size_t numChannels = input2D.size();
  size_t numTicks = input2D[0].size();
  std::complex<float> I(0, 1);

  if (output2D.size() != numChannels) throw std::runtime_error("CSUBilateral: output array does not match numChannels");
  if (output2D[0].size() != numTicks) throw std::runtime_error("CSUBilateral: output length does not match numTicks");

  float gamma = M_PI / (2.0 * ((float) T));
  float rho = gamma * sigma_r;
  size_t N = 3;
  if (sigma_r < (1.0 / gamma)) {
    N = (size_t) std::ceil( std::pow(gamma * sigma_r, -2.0));
    std::cout << "N = " << N << std::endl;
    std::cout << "Gamma = " << gamma << std::endl;
    std::cout << "sigma_r = " << sigma_r << std::endl;
    N = 5;
  }

  std::cout << "Initializing Auxiliary Images..." << std::endl;
  
  Array3D<std::complex<float>> H(N+1);
  Array3D<std::complex<float>> G(N+1);
  Array3D<std::complex<float>> D(N+1);

  Array3D<std::complex<float>> Hbar(N+1);
  Array3D<std::complex<float>> Gbar(N+1);

  for (size_t n=0; n<=N; ++n) {
    H[n].resize(numChannels);
    G[n].resize(numChannels);
    D[n].resize(numChannels);
    Hbar[n].resize(numChannels);
    Gbar[n].resize(numChannels);
    for (size_t i=0; i<numChannels; ++i) {
      H[n][i].resize(numTicks);
      G[n][i].resize(numTicks);
      D[n][i].resize(numTicks);
      Hbar[n][i].resize(numTicks);
      Gbar[n][i].resize(numTicks);
    }
  }

  icarus_signal_processing::MiscUtils utils;

  float w = gamma / (rho * sqrt((float) N));
  float c0 = std::pow(2.0, -(float)N);

  std::cout << "Computing Auxiliary Images..." << std::endl;

  // Set up auxiliary images H, G, D;
  for (size_t n=0; n<=N; ++n) {
    float binom = (float) utils.nChoosek(N, n);
    for (size_t i=0; i<numChannels; ++i) {
      for (size_t j=0; j<numTicks; ++j) {
        float q = (float) (2 * (int) n - (int) N);
        std::complex<float> z = q * I * w * input2D[i][j];
        H[n][i][j] = std::exp(z);
        G[n][i][j] = H[n][i][j] * input2D[i][j];
        D[n][i][j] = c0 * binom * std::exp(-z);
      }
    }
  }

  std::cout << "Running Gaussian Filtering..." << std::endl;

  // Run O(1) Gaussian Filtering
  for (size_t n=0; n<=N; ++n) {
    Gaussian(H[n], sigma_s, sigma_s, Hbar[n]);
    Gaussian(G[n], sigma_s, sigma_s, Gbar[n]);
  }

  std::cout << "Computing Output..." << std::endl;

  // Compute output
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<numTicks; ++j) {
      std::complex<float> summedGbarX;
      std::complex<float> summedHbarX;
      std::complex<float> fbarX;
      for (size_t n=0; n<=N; ++n) {
        summedGbarX += D[n][i][j] * Gbar[n][i][j];
        summedHbarX += D[n][i][j] * Hbar[n][i][j];
      }
      fbarX = summedGbarX / summedHbarX;
      output2D[i][j] = fbarX.real();
    }
  }
  return;

}


void icarus_signal_processing::BilateralFilters::GPABilateral(
  const Array2D<float>& input2D,
  Array2D<float>& output2D,
  const unsigned int T,
  const float sigma_x,
  const float sigma_y,
  const float sigma_r,
  const float epsilon) const
{
  size_t numChannels = input2D.size();
  size_t numTicks = input2D[0].size();

  if (output2D.size() != numChannels) throw std::runtime_error("GPABilateral: output array does not match numChannels");
  if (output2D[0].size() != numTicks) throw std::runtime_error("GPABilateral: output length does not match numTicks");

  unsigned int N0 = estimateGPAOrder(sigma_r, epsilon, (float) T);
  std::cout << "N0 = " << N0 << std::endl;

  Array2D<float> h(numChannels);
  Array2D<float> F(numChannels);
  Array2D<float> G(numChannels);
  Array2D<float> P(numChannels);
  Array2D<float> Q(numChannels);
  Array2D<float> H(numChannels);
  Array2D<float> Fbar(numChannels);

  for (auto& v : h) {
    v.resize(numTicks);
  }
  for (auto& v : F) {
    v.resize(numTicks);
  }
  for (auto& v : G) {
    v.resize(numTicks);
  }
  for (auto& v : P) {
    v.resize(numTicks);
  }
  for (auto& v : Q) {
    v.resize(numTicks);
  }
  for (auto& v : H) {
    v.resize(numTicks);
  }
  for (auto& v : Fbar) {
    v.resize(numTicks);
  }


  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<numTicks; ++j) {
      h[i][j] = input2D[i][j] - T;
      F[i][j] = exp(-h[i][j]*h[i][j] / (2.0*sigma_r*sigma_r));
      G[i][j] = 1.0;
      P[i][j] = 0.0;
      Q[i][j] = 0.0;
      H[i][j] = h[i][j] / sigma_r;
    }
  }
  Gaussian(F, sigma_x, sigma_y, Fbar);
  for (size_t n=1; n<=N0; ++n) {
    for (size_t i=0; i<numChannels; ++i) {
      for (size_t j=0; j<numTicks; ++j) {
        Q[i][j] += G[i][j] * Fbar[i][j];
        F[i][j] = H[i][j] * F[i][j];
      }
    }
    Gaussian(F, sigma_x, sigma_y, Fbar);
    for (size_t i=0; i<numChannels; ++i) {
      for (size_t j=0; j<numTicks; ++j) {
        P[i][j] += G[i][j] * Fbar[i][j];
        G[i][j] = H[i][j] * G[i][j] / ((float) n);
      }
    }
  }
  for (size_t i=0; i<numChannels; ++i) {
    for (size_t j=0; j<numTicks; ++j) {
      output2D[i][j] = sigma_r * (P[i][j] / Q[i][j]) + (float) T;
    }
  }
  return;
}


unsigned int icarus_signal_processing::BilateralFilters::estimateGPAOrder(
  float sigma_r,
  float epsilon,
  float T) const
{
  if (sigma_r >= 70) {
    unsigned int N0 = 10;
    return N0;
  }
  else {
    float lambda = (T / sigma_r) * (T / sigma_r);
    float p = 1 + log(lambda);
    float q = -lambda - log(epsilon);
    float t = q / (exp(1.0) * lambda);
    float W0 = t - t*t + 3*std::pow(t, 3.0) / 2.0 - 8.0*std::pow(t, 4.0) / 4.0;
    float N = q / W0;
    if (sigma_r < 30) {
      for (int k=1; k<4; ++k) {
        N = N - (N * log(N) - p* N - q) / (log(N) + 1 - p);
      }
    }
    unsigned int N0 = (unsigned int) std::round(N);
    return N0;
  }
}


#endif
