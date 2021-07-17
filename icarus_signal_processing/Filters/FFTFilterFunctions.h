/**
 * \file FFTFilterFunctions.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for objects to perform FFT filtering (high, low pass filters)
 *
 * @author usher@slac.stanford.edu
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __SIGPROC_TOOLS_FFTFILTERFUNCTIONS_H__
#define __SIGPROC_TOOLS_FFTFILTERFUNCTIONS_H__

#include <vector>
#include <memory>
#include <complex>

namespace icarus_signal_processing 
{
/**
 * \class IFFTFilterFunctions
 * 
 * \ingroup icarus_signal_processing
 * 
 * \brief Interface class for objects
 */

template <class T> class ICARUSFFT;
template <class T> using Waveform = std::vector<T>;

class IFFTFilterFunction
{
public:
   /**
   *  @brief  Virtual Destructor
   */
   virtual ~IFFTFilterFunction() noexcept = default;
 
   /**
   *  @brief Interface functions which provided templated access
   *
   *  @param Waveform  The waveform to process
   */
   virtual void operator()(Waveform<float>&)  const = 0;
   virtual void operator()(Waveform<double>&) const = 0;
};

using FFTFilterFunctionVec = std::vector<std::unique_ptr<IFFTFilterFunction>>; 

/**
* \class High Pass Filter 
* 
* \ingroup icarus_signal_processing
* 
* \brief Performs a high pass filter on incoming waveform
*/
class HighPassFFTFilter : virtual public IFFTFilterFunction
{
public:
    /**
     *  @brief  Constructor
     */
    explicit HighPassFFTFilter(const double&, const double&);

    /**
     *  @brief  Destructor
     */
    ~HighPassFFTFilter();
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(Waveform<float>&)  const override;
    void operator()(Waveform<double>&) const override;

private:
    void highPassFilter(Waveform<float>&) const;

    using KernelVec = std::vector<std::complex<float>>;

    KernelVec fHighPassFilterKernel;

    // Keep track of the FFT 
    std::unique_ptr<icarus_signal_processing::ICARUSFFT<float>> fFFT; ///< Object to handle thread safe FFT
};

/**
* \class Low Pass Filter
* 
* \ingroup icarus_signal_processing
* 
* \brief Performs a low pass filter on incoming waveform
*/
class LowPassFFTFilter : virtual public IFFTFilterFunction
{
public:
    /**
     *  @brief  Constructor
     */
    explicit LowPassFFTFilter()  {}

    /**
     *  @brief  Destructor
     */
    ~LowPassFFTFilter() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(Waveform<float>&)  const override;
    void operator()(Waveform<double>&) const override;

private:
    void lowPassFilter(Waveform<float>&) const;

    using KernelVec = std::vector<Waveform<float>>;

    KernelVec fFilterKernelVec;
};

/**
* \class Window Filter
* 
* \ingroup icarus_signal_processing
* 
* \brief Performs a windowing filter on incoming waveform
*/
class WindowFFTFilter : virtual public IFFTFilterFunction
{
public:
    /**
     *  @brief  Constructor
     */
    explicit WindowFFTFilter(const std::pair<double,double>&, const std::pair<double,double>&);

    /**
     *  @brief  Destructor
     */
    ~WindowFFTFilter();
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(Waveform<float>&)  const override;
    void operator()(Waveform<double>&) const override;

private:
    void windowFilter(Waveform<float>&) const;

    using KernelVec = std::vector<std::complex<float>>;

    KernelVec fWindowFilterKernel;

    // Keep track of the FFT 
    std::unique_ptr<icarus_signal_processing::ICARUSFFT<float>> fFFT; ///< Object to handle thread safe FFT
};

/**
* \class No Filter 
*  
* \ingroup icarus_signal_processing
* 
* \brief This does not filtering (but allows to be included in a list of filter functions)
*/
class NoFFTFilter : virtual public IFFTFilterFunction
{
public:
    /**
     *  @brief  Constructor
     */
    explicit NoFFTFilter() {}

    /**
     *  @brief  Destructor
     */
    ~NoFFTFilter() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(Waveform<float>&)  const override {return;}
    void operator()(Waveform<double>&) const override {return;}

private:
};

/**
* \class low pass Butterworth filter
*  
* \ingroup icarus_signal_processing
* 
* \brief Implements a low pass Butterworth filter
*/
class LowPassButterworthFilter : virtual public IFFTFilterFunction
{
public:
    /**
     *  @brief  Constructor
     */
    explicit LowPassButterworthFilter(unsigned int, unsigned int, unsigned int);

    /**
     *  @brief  Destructor
     */
    ~LowPassButterworthFilter();
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(Waveform<float>&)  const override;
    void operator()(Waveform<double>&) const override;

private:
    void lowPassButterworthFilter(Waveform<float>&) const;

    unsigned int                     fThreshold;
    unsigned int                     fOrder;
    unsigned int                     fFrequencySize;
    std::vector<std::complex<float>> fFilterVec;

    // Keep track of the FFT 
    std::unique_ptr<icarus_signal_processing::ICARUSFFT<float>> fFFT; ///< Object to handle thread safe FFT
};

/**
* \class high pass Butterworth filter
*  
* \ingroup icarus_signal_processing
* 
* \brief Implements a low pass Butterworth filter
*/
class HighPassButterworthFilter : virtual public IFFTFilterFunction
{
public:
    /**
     *  @brief  Constructor
     */
    explicit HighPassButterworthFilter(unsigned int, unsigned int, unsigned int);

    /**
     *  @brief  Destructor
     */
    ~HighPassButterworthFilter();
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(Waveform<float>&)  const override;
    void operator()(Waveform<double>&) const override;

private:
    void highPassButterworthFilter(Waveform<float>&) const;

    unsigned int                     fThreshold;
    unsigned int                     fOrder;
    unsigned int                     fFrequencySize;
    std::vector<std::complex<float>> fFilterVec;

    // Keep track of the FFT 
    std::unique_ptr<icarus_signal_processing::ICARUSFFT<float>> fFFT; ///< Object to handle thread safe FFT
};

}

#endif
/** @} */ // end of doxygen group 
