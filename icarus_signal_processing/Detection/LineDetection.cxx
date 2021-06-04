#ifndef __icarus_signal_processing_LINEDETECTION_CXX__
#define __icarus_signal_processing_LINEDETECTION_CXX__

#include "LineDetection.h"
#include "MorphologicalFunctions2D.h"


void icarus_signal_processing::LineDetection::HoughTransform(
  const Array2D<bool>& binary2D,
  Array2D<int>& accumulator2D,
  const unsigned int thetaSteps) const
{
  const double pi = 3.141592653589793238462643383279502884;

  int width = binary2D.size();
  int height = binary2D.at(0).size();

//  std::cout << "width (i) = " << width << std::endl;
//  std::cout << "height (j) = " << height << std::endl;

  float dtheta = ((float) pi) / ((float) thetaSteps);
  int diagLength = (int) std::round(std::sqrt(width*width + height*height));

  accumulator2D.resize(2 * diagLength);
  for (auto& v : accumulator2D) {
    v.resize(thetaSteps);
  }
  // Precompute trig table for optimization. 
  std::vector<float> trigtab(thetaSteps * 2);
  for (size_t i=0; i<thetaSteps; ++i) {
    trigtab[i * 2] = (float) std::cos(i * dtheta);
    trigtab[i * 2+1] = (float) std::sin(i * dtheta);
  }
  // // Collect non-zero points from binary image.
  std::vector<Point> edgeList;
  for (int i=0; i<width; ++i) {
    for (int j=0; j<height; ++j) {
      if (binary2D[i][j]) {
        Point pt = {i, j};
        edgeList.push_back(pt);
      }
    }
  }

  float r = 0;
  int ir = 0;

  // Update accumulator
  for (Point& pt : edgeList) {
    for (size_t itheta=0; itheta<thetaSteps; ++itheta) {
      // std::cout << pt.x << " " << pt.y << std::endl;
      r = (float) (pt.x * trigtab[itheta*2] + pt.y * trigtab[itheta*2 + 1]);
      ir = diagLength + (int) std::round(r);
      if (std::abs( (int) ir ) > 2 * diagLength) {
        std::cout << ir << std::endl;
      }
      accumulator2D[ir][itheta] += 1;
    }
  }

  return;
}


int icarus_signal_processing::LineDetection::CartesianHoughTransform(
  const Array2D<bool>& binary2D,
  Array2D<int>& accumulator2D,
  const float maxAngleDev,
  const int thetaSteps) const
{
  const double pi = 3.141592653589793238462643383279502884;

  // assert(maxAngleDev < 90);

  const int numChannels = binary2D.size();
  const int numTicks = binary2D.at(0).size();

  const float maxAngle = ((float) maxAngleDev) / ((float) 180) * (float) pi;
  const float dtheta = 2.0 * maxAngle / ((float) 2 * thetaSteps);

  std::vector<float> slopeTab(thetaSteps * 2);
  for (int i=0; i< thetaSteps; ++i) {
    // 0 ~ maxAngle
    slopeTab[i * 2] = (float) std::tan(i * dtheta);
    // 0 ~ -maxAngle
    slopeTab[i * 2+1] = (float) -std::tan(i * dtheta);
  }

  int padding = ((int) numChannels * std::tan(maxAngle) + 1);

  // Accumulator for slope and tick-intercept 
  accumulator2D.resize(2 * thetaSteps);
  for (auto& v : accumulator2D) {
    v.resize(numTicks + 2 * padding);
  }

  for (int i=0; i<numChannels; ++i) {
    for (int j=0; j<numTicks; ++j) {
      if (binary2D[i][j]) {
        // 0 ~ -maxAngle
        for (int m=1; m<thetaSteps; ++m) {
          const int intercept = (int) (slopeTab[2*(thetaSteps-m)+1] * ((float) i));
          accumulator2D[m][padding + j - intercept] += 1;
        }
        for (int m=0; m<thetaSteps; ++m) {
          // Process positive angles
          const int intercept = (int) (slopeTab[2*m] * ((float) i));
          accumulator2D[thetaSteps + m][padding + j - intercept] += 1;
        }
      }
    }
  }

  // Run NMS
  return padding;
}



void icarus_signal_processing::LineDetection::simpleFastNMS(
  const Array2D<int>& accumulator2D,
  std::vector<int>& rhoIndex,
  std::vector<int>& thetaIndex,
  const int threshold,
  const int sx,
  const int sy) const
{
  simpleFastNMS<int>(accumulator2D, rhoIndex, thetaIndex, threshold, sx, sy);
}

void icarus_signal_processing::LineDetection::simpleFastNMS(
  const Array2D<long>& accumulator2D,
  std::vector<int>& rhoIndex,
  std::vector<int>& thetaIndex,
  const int threshold,
  const int sx,
  const int sy) const
{
  simpleFastNMS<long>(accumulator2D, rhoIndex, thetaIndex, threshold, sx, sy);
}

template <typename T>
void icarus_signal_processing::LineDetection::simpleFastNMS(
  const Array2D<T>& accumulator2D,
  std::vector<int>& rhoIndex,
  std::vector<int>& thetaIndex,
  const T threshold,
  const int fStructuringElementX,
  const int fStructuringElementY) const
{
  const int numTheta = accumulator2D.size();
  const int numIntercept = accumulator2D.at(0).size();

  icarus_signal_processing::Morph2DFast morph2d;

  Array2D<T> tempBuffer(numTheta);
  for (auto& v : tempBuffer) {
    v.resize(numIntercept);
  }

  // The local maximum of the accumulator array is obtained by
  // selecting points whose dilation matches the original image.

  icarus_signal_processing::Dilation2D(fStructuringElementX, fStructuringElementY)(accumulator2D, numTheta, tempBuffer.begin());

  std::unordered_map<long, LabeledPoint> edges;

  int count_id = 1;

  for (int i=0; i<numTheta; ++i) {
    for (int j=0; j<numIntercept; ++j) {
      if ((std::abs(accumulator2D[i][j] - tempBuffer[i][j]) < 0.0001) && 
            (accumulator2D[i][j] > (T) threshold)) {
        long cantorEnum = ((i + j) * (i + j + 1) / 2) + j;
        LabeledPoint p = {i, j, count_id};
        rhoIndex.push_back(j);
        thetaIndex.push_back(i);
        edges.emplace(std::make_pair(cantorEnum, p));
        count_id++;
      }
    }
  }
  return;
}


void icarus_signal_processing::LineDetection::FindPeaksNMS(
  const Array2D<int>& accumulator2D,
  std::vector<int>& rhoIndex,
  std::vector<int>& thetaIndex,
  const unsigned int threshold,
  const unsigned int angleWindow,
  const unsigned int maxLines,
  const unsigned int windowSize) const
{

  // std::cout << maxLines << std::endl;

  rhoIndex.reserve(maxLines);
  thetaIndex.reserve(maxLines);

  int rhoSteps = accumulator2D.size();
  int thetaSteps = accumulator2D.at(0).size();
  for (int i=0; i < rhoSteps; ++i) {
    for (int j=0; j < thetaSteps; ++j) {
      int lbx = i - (int) windowSize;
      int ubx = i + (int) windowSize;
      int lby = j - (int) windowSize;
      int uby = j + (int) windowSize;
      int lowerBoundx = std::max(lbx, 0);
      int upperBoundx = std::min(ubx, (int) rhoSteps);
      int lowerBoundy = std::max(lby, 0);
      int upperBoundy = std::min(uby, (int) thetaSteps);
      int centerVal = accumulator2D[i][j];
      if (centerVal < (int) threshold) continue;
      int angleLowerBound = thetaSteps / 2 - angleWindow;
      int angleUpperBound = thetaSteps / 2 + angleWindow;
      if (angleLowerBound > j or angleUpperBound < j) {
        continue;
      }
      else {
        bool isNonMax = false;
        for (int ix=lowerBoundx; ix<upperBoundx; ++ix) {
          for (int iy=lowerBoundy; iy<upperBoundy; ++iy) {
            int currentVal = accumulator2D[ix][iy];
            if (currentVal > centerVal) {
              isNonMax = true;
            }
            if (isNonMax) break;
          }
          if (isNonMax) break;
        }
        if (isNonMax) {
          continue;
        } else {
          rhoIndex.push_back(i);
          thetaIndex.push_back(j);
        }
      }
    }
  }
  return;
}


void icarus_signal_processing::LineDetection::FindVerticalSegments(
  const Array2D<bool>& selectVals,
  Array2D<bool>& newMask2D,
  const size_t thetaSteps,
  const std::vector<int>& rhoIndex,
  const std::vector<int>& thetaIndex,
  const unsigned int maxGap,
  const unsigned int minLineLength,
  const unsigned int margin) const
{
  // Not yet completed, use ProbabilisticHT for now
  const double pi = 3.141592653589793238462643383279502884;

  int numChannels = selectVals.size();
  int numTicks = selectVals.at(0).size();

  newMask2D.resize(numChannels);
  for (auto& v : newMask2D) {
    v.resize(numTicks);
  }

  if (rhoIndex.size() == 0) {
    return;
  }

  float dtheta = (float) pi / (float) thetaSteps;

  for (size_t i=0; i<rhoIndex.size(); ++i) {
    float theta = dtheta * (float) thetaIndex[i] - (float) pi / 2.0; 
    int xIntercept = (int) std::round(rhoIndex[i] * 1.0 / std::cos(theta));
    int x = xIntercept;
    int y = 0;
    float error = 0.0;

//    std::cout << "theta = " << theta << std::endl;

    float slope = -1.0 * std::tan(theta);
    int sign = (slope > 0) - (slope < 0);
    // Since we are only interested in near vertical segments, we move along
    // the y-axis.
    float delta = std::abs(slope);

//    std::cout << "Slope = " << slope << std::endl;
//    std::cout << "Sign = " << sign << std::endl;

    std::vector<int> tempVecx;
    std::vector<int> tempVecy;

    tempVecx.reserve(maxGap);
    tempVecy.reserve(maxGap);

    unsigned int skipCounts = 0;
    bool vecEmpty = true;

    while (x < numChannels && y < numTicks) {
      // std::cout << "x = " << x << " | " << "y = " << y << std::endl;
      bool pixel = selectVals[x][y];
      if (pixel && vecEmpty) {
        // Currently at beginning of candidate line segment.
        newMask2D[x][y] = true;
        vecEmpty = false;
        skipCounts = 0;
      }
      else if (!pixel && vecEmpty) {
        y++;
        error += delta;
        if (error > 0.5) {
          x += sign;
          error = error - 1.0;
        }
        continue;
      }
      else if (pixel && !vecEmpty) {
        for (size_t v=0; v<tempVecx.size(); ++v) {
          newMask2D[tempVecx[v]][tempVecy[v]] = true;
        }
        tempVecx.clear();
        tempVecy.clear();
        vecEmpty = true;
        skipCounts = 0;
        newMask2D[x][y] = true;
      }
      else {
        if (skipCounts < maxGap) {
          tempVecx.push_back(x);
          tempVecy.push_back(y);
          skipCounts++;
        } 
        else {
          tempVecx.clear();
          tempVecy.clear();
          vecEmpty = true;
          skipCounts = 0;
        }
      }
      y++;
      error += delta;
      if (error > 0.5) {
        x += sign;
        error = error - 1.0;
      }
    }
  }
//  std::cout << margin << std::endl;
//  std::cout << minLineLength << std::endl;
  return;
}


icarus_signal_processing::Line icarus_signal_processing::LineDetection::getLine(
  const Array2D<bool>& binary2D,
  const float theta,
  const float r,
  const float eps,
  const unsigned int angleWindow) const
{

  int width = binary2D.size();
  int height = binary2D.at(0).size();

  const double pi = 3.141592653589793238462643383279502884;
  int theta_deg = (int) std::round(theta * 180.0 / (float) pi);
  // Ignore lines which are near parallel to tick axis (these don't cause problems)
  if (std::abs(theta_deg) < angleWindow || std::abs(180 - theta_deg) < angleWindow) {
    Line l = {-1, -1, 0};
    return l;
  }
  else {
    float slope = -1.0 / std::tan(theta);
    int xIntercept = (int) r / (std::cos(theta) + eps);
    int yIntercept = (int) r / (std::sin(theta) + eps);
    // Bounds
    Point pt1 = {(int) xIntercept, 0};
    Point pt2 = {0, (int) yIntercept};
    Point pt3 = { (width - 1), (int) (slope * float(width - 1) + (float) yIntercept) };
    Point pt4 = { (int) ( (float) (height - 1 - yIntercept) / slope) , height-1};

    std::vector<Point> pts = {pt1, pt2, pt3, pt4};

    int x0 = width;
    int y0 = height;

    if (std::abs(slope) < 1.0) {
      // Then we draw lines by incrementing x, so we find the smallest x
      // index inside the given image.
      for (auto& pt : pts) {
        if ( !((0 <= pt.x) && (pt.x < width) && (0 <= pt.y) && (pt.y < height)) ) {
          continue;
        }
        if (pt.x <= x0) {
          x0 = pt.x;
          y0 = pt.y;
        } 
      }
    }
    else {
      // Then we draw lines by incrementing y, so we find the smallest y
      // index inside the given image.
      for (auto& pt : pts) {
        if ( !((0 <= pt.x) && (pt.x < width) && (0 <= pt.y) && (pt.y < height)) ) {
          continue;
        }
        if (pt.y <= y0) {
          x0 = pt.x;
          y0 = pt.y;
        } 
      }
    }
    Line l = { (int) x0, (int) y0, slope};
    return l;
  }
}


void icarus_signal_processing::LineDetection::drawLine(
  Array2D<bool>& newSelectVals,
  const Line line,
  const unsigned int dilationX,
  const unsigned int dilationY) const
{
  int x = line.x0;
  int y = line.y0;
  float slope = line.slope;
  float error = 0.0;
  float derror = std::abs(slope);
  int sign =  ( (int) (slope > 0) - (int) (slope < 0) );
  int width = newSelectVals.size();
  int height = newSelectVals.at(0).size();

  // Bresenham's Line Drawing Algorithm

  if (derror < 1.0) {
    while ((0 <= x) && (x < width) && (0 <= y) && (y < height)) {
      newSelectVals[x][y] = true;
      int lbx = (int) x - (int) dilationX;
      int ubx = (int) x + (int) dilationX;
      int lby = (int) y - (int) dilationY;
      int uby = (int) y + (int) dilationY;
      int lowerBoundx = std::max(lbx, 0);
      int upperBoundx = std::min(ubx, (int) width);
      int lowerBoundy = std::max(lby, 0);
      int upperBoundy = std::min(uby, (int) height);
      for (int i=lowerBoundx; i<upperBoundx; ++i) {
        for (int j=lowerBoundy; j<upperBoundy; ++j) {
          newSelectVals[i][j] = true;
        }
      }
      error += derror;
      ++x;
      if (error >= 0.5) {
        y += sign;
        error = error - 1.0;
      }
    }
  }
  else {
    derror = 1.0 / derror;
    while ((0 <= x) && (x < width) && (0 <= y) && (y < height)) {
      newSelectVals[x][y] = true;
      int lbx = (int) x - (int) dilationX;
      int ubx = (int) x + (int) dilationX;
      int lby = (int) y - (int) dilationY;
      int uby = (int) y + (int) dilationY;
      int lowerBoundx = std::max(lbx, 0);
      int upperBoundx = std::min(ubx, (int) width);
      int lowerBoundy = std::max(lby, 0);
      int upperBoundy = std::min(uby, (int) height);
      for (int i=lowerBoundx; i<upperBoundx; ++i) {
        for (int j=lowerBoundy; j<upperBoundy; ++j) {
          newSelectVals[i][j] = true;
        }
      }
      error += derror;
      ++y;
      if (error >= 0.5) {
        x += sign;
        error = error - 1.0;
      }
    }
  }
  return;
}


void icarus_signal_processing::LineDetection::drawLine2(
  Array2D<bool>& newSelectVals,
  const int &interceptIndex,
  const float &theta, //radians
  const int &padding) const
{
  const int numChannels = newSelectVals.size();
  const int numTicks = newSelectVals.at(0).size();

  int channelOffset = 0;
  int tickAxisIntercept = 0;

  float eps = 0.001;
  float slope = std::tan(theta);
  int sign =  ( (int) (slope > 0) - (int) (slope < 0) );

  // Compute tick axis intercept
  if ( (padding < interceptIndex) && (interceptIndex <= (numTicks + padding) ))
  {  
    tickAxisIntercept = interceptIndex - padding;
  } 
  else if (padding > interceptIndex) 
  {
    channelOffset = ( (float) (padding - interceptIndex) ) / (std::abs(slope) + eps);
    tickAxisIntercept = 0;
  }
  else {
    tickAxisIntercept = numTicks;
    channelOffset = ( (float) (interceptIndex - numTicks - padding)) / (std::abs(slope) + eps);
  }


  // Bresenham's Line Drawing Algorithm
  // Since we are only interested in isochronous tracks, safe to assume
  // change in angle is greater than change in ticks. 

  int ix = channelOffset;
  int iy = tickAxisIntercept;

  float error = 0.0;
  float slopeAbs = std::abs(slope);

  while ( (ix < numChannels) && (iy < numTicks) && (iy >= 0) ) {
    newSelectVals[ix][iy] = true;
    error += slopeAbs;
    ix++; 
    if (error > 0.5) {
      iy += sign;
      error = error - 1.0;
    }
  }

  return;

}


void icarus_signal_processing::LineDetection::refineSelectVals(
  const Array2D<bool>& selectVals,
  Array2D<bool>& refinedSelectVals,
  const size_t thetaSteps,
  const unsigned int threshold,
  const unsigned int angleWindow,
  const unsigned int maxLines,
  const unsigned int windowSize,
  const unsigned int dilationX,
  const unsigned int dilationY,
  const float eps) const
{
  int numChannels = selectVals.size();
  int numTicks = selectVals.at(0).size();
  const double pi = 3.141592653589793238462643383279502884;

  Array2D<int> accumulator2D;
  HoughTransform(selectVals, accumulator2D, thetaSteps);

  std::vector<int> rhoIndex;
  std::vector<int> thetaIndex;

  FindPeaksNMS(accumulator2D, rhoIndex, thetaIndex, 
               threshold, angleWindow, maxLines, windowSize);

  int diagLength = (int) std::round(
    std::sqrt(numChannels * numChannels + numTicks * numTicks));

  float rho = 0.0;
  float angle = 0.0;

//  std::cout << rhoIndex.size() << std::endl;
//  std::cout << thetaIndex.size() << std::endl;

  for (size_t i=0; i<rhoIndex.size(); ++i) {
    rho = rhoIndex[i] - diagLength;
    angle = thetaIndex[i] * pi / ( (float) thetaSteps );
    Line l = getLine(selectVals, angle, rho, eps, angleWindow);
//    std::cout << "x0 = " << l.x0 << ", y0 = " << l.y0 << ", slope = " << l.slope << std::endl;
    drawLine(refinedSelectVals, l, dilationX, dilationY);
  }
  return;
}

#endif
