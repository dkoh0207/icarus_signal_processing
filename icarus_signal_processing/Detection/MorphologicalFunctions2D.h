/**
 * \file MorphologicalFunctions2D.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for objects to perform various morphological filter operations
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __SIGPROC_TOOLS_MorphologicalFunctions2D_H__
#define __SIGPROC_TOOLS_MorphologicalFunctions2D_H__

//#include "ICARUSSigProcDefs.h"
#include <vector>
#include <memory>

#include "icarus_signal_processing/ICARUSSigProcDefs.h"

namespace icarus_signal_processing 
{
/**
 * \class IMorphologicalFunctions2D
 * 
 * \ingroup icarus_signal_processing
 * 
 * \brief Interface class for objects
 */

template <typename T> using WaveArray = std::vector<std::vector<T>>;

class IMorphologicalFunctions2D
{
public:
   /**
   *  @brief  Virtual Destructor
   */
   virtual ~IMorphologicalFunctions2D() noexcept = default;
 
   /**
   *  @brief Interface functions which provided templated access
   *
   *  @param fragment            The artdaq fragment to process
   */
   virtual void operator()(ArrayFloat::const_iterator, const unsigned int, ArrayFloat::iterator)  const = 0;
   virtual void operator()(ArrayBool::const_iterator,  const unsigned int, ArrayBool::iterator )  const = 0;
};

using Filter2DFunctionVec = std::vector<std::unique_ptr<IMorphologicalFunctions2D>>; 

/**
* \class Dilation2D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the dilation of a waveform
*/
class Dilation2D : virtual public IMorphologicalFunctions2D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Dilation2D(const unsigned int structuringElementX, const unsigned int structuringElementY) : 
                        fStructuringElementX(structuringElementX), fStructuringElementY(structuringElementY) 
             {}

    /**
     *  @brief  Destructor
     */
    ~Dilation2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(ArrayFloat::const_iterator,  const unsigned int, ArrayFloat::iterator)  const override;
    void operator()(ArrayBool::const_iterator,   const unsigned int, ArrayBool::iterator )  const override;

private:
    template <typename T> void getDilation2D(typename WaveArray<T>::const_iterator, const unsigned int, typename WaveArray<T>::iterator) const;

    unsigned int fStructuringElementX;
    unsigned int fStructuringElementY;
};

/**
* \class Erosion2D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the erosion of a waveform
*/
class Erosion2D : virtual public IMorphologicalFunctions2D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Erosion2D(const unsigned int structuringElementX, const unsigned int structuringElementY) : 
                        fStructuringElementX(structuringElementX), fStructuringElementY(structuringElementY) 
             {}

    /**
     *  @brief  Destructor
     */
    ~Erosion2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(ArrayFloat::const_iterator,  const unsigned int, ArrayFloat::iterator)  const override;
    void operator()(ArrayBool::const_iterator,   const unsigned int, ArrayBool::iterator )  const override;

private:
    template <typename T> void getErosion2D(typename WaveArray<T>::const_iterator, const unsigned int, typename WaveArray<T>::iterator) const;

    unsigned int fStructuringElementX;
    unsigned int fStructuringElementY;
};

/**
* \class Gradient2D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the gradient of a wabeform
*/
class Gradient2D : virtual public IMorphologicalFunctions2D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Gradient2D(const unsigned int structuringElementX, const unsigned int structuringElementY) : 
                        fStructuringElementX(structuringElementX), fStructuringElementY(structuringElementY) 
             {}

    /**
     *  @brief  Destructor
     */
    ~Gradient2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(ArrayFloat::const_iterator,  const unsigned int, ArrayFloat::iterator)  const override;
    void operator()(ArrayBool::const_iterator,   const unsigned int, ArrayBool::iterator )  const override;

private:
    template <typename T> void getGradient2D(typename WaveArray<T>::const_iterator, const unsigned int, typename WaveArray<T>::iterator) const;

    unsigned int fStructuringElementX;
    unsigned int fStructuringElementY;
};

/**
* \class Average2D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the average of a waveform
*/
class Average2D : virtual public IMorphologicalFunctions2D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Average2D(const unsigned int, const unsigned int)
//    explicit Average2D(const unsigned int structuringElementX, const unsigned int structuringElementY)  :
//                        fStructuringElementX(structuringElementX), fStructuringElementY(structuringElementY) 
             {}

    /**
     *  @brief  Destructor
     */
    ~Average2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(ArrayFloat::const_iterator,  const unsigned int, ArrayFloat::iterator)  const override;
    void operator()(ArrayBool::const_iterator,   const unsigned int, ArrayBool::iterator )  const override;

private:
    template <typename T> void getAverage2D(typename WaveArray<T>::const_iterator, const unsigned int, typename WaveArray<T>::iterator) const;

//    unsigned int fStructuringElementX;
//    unsigned int fStructuringElementY;
};

/**
* \class Median2D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the median of a waveform
*/
class Median2D : virtual public IMorphologicalFunctions2D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Median2D(const unsigned int structuringElementX, const unsigned int structuringElementY, const unsigned int pixelRange) :
                fStructuringElementX(structuringElementX), fStructuringElementY(structuringElementY), fPixelRange(pixelRange) 
                {}

    /**
     *  @brief  Destructor
     */
    ~Median2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(ArrayFloat::const_iterator,  const unsigned int, ArrayFloat::iterator)  const override;
    void operator()(ArrayBool::const_iterator,   const unsigned int, ArrayBool::iterator )  const override;

private:
    template <typename T> void getMedian2D(typename WaveArray<T>::const_iterator, const unsigned int, typename WaveArray<T>::iterator) const;

    unsigned int fStructuringElementX;
    unsigned int fStructuringElementY;
    unsigned int fPixelRange;
};

/**
* \class Opening2D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the opening of a waveform
*/
class Opening2D : virtual public IMorphologicalFunctions2D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Opening2D(const unsigned int structuringElementX, const unsigned int structuringElementY) : 
                        fStructuringElementX(structuringElementX), fStructuringElementY(structuringElementY) 
             {}

    /**
     *  @brief  Destructor
     */
    ~Opening2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(ArrayFloat::const_iterator waveformIn,  const unsigned int nChannels, ArrayFloat::iterator opening)  const override
    {
        ArrayFloat temp(nChannels);

        Erosion2D (fStructuringElementX,fStructuringElementY)(waveformIn, nChannels, temp.begin());
        Dilation2D(fStructuringElementX,fStructuringElementY)(temp.begin(),       nChannels, opening);
    }

    void operator()(ArrayBool::const_iterator waveformIn, const unsigned int nChannels, ArrayBool::iterator opening) const override
    {
        ArrayBool temp(nChannels);

        Erosion2D (fStructuringElementX,fStructuringElementY)(waveformIn, nChannels, temp.begin());
        Dilation2D(fStructuringElementX,fStructuringElementY)(temp.begin(),       nChannels, opening);
    }

private:
    unsigned int fStructuringElementX;
    unsigned int fStructuringElementY;

};

/**
* \class Closing2D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the closing of a waveform
*/
class Closing2D : virtual public IMorphologicalFunctions2D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Closing2D(const unsigned int structuringElementX, const unsigned int structuringElementY) : 
                        fStructuringElementX(structuringElementX), fStructuringElementY(structuringElementY) 
             {}

    /**
     *  @brief  Destructor
     */
    ~Closing2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(ArrayFloat::const_iterator waveformIn, const unsigned int nChannels, ArrayFloat::iterator closing)  const override
    {
        ArrayFloat temp(nChannels);

        Dilation2D(fStructuringElementX,fStructuringElementY)(waveformIn, nChannels, temp.begin());
        Erosion2D( fStructuringElementX,fStructuringElementY)(temp.begin(),       nChannels, closing);
    }

    void operator()(ArrayBool::const_iterator waveformIn, const unsigned int nChannels, ArrayBool::iterator closing)  const override
    {
        ArrayBool temp(nChannels);

        Dilation2D(fStructuringElementX,fStructuringElementY)(waveformIn, nChannels, temp.begin());
        Erosion2D( fStructuringElementX,fStructuringElementY)(temp.begin(),       nChannels, closing);
    }

private:
    unsigned int fStructuringElementX;
    unsigned int fStructuringElementY;

};

}

#endif
/** @} */ // end of doxygen group 