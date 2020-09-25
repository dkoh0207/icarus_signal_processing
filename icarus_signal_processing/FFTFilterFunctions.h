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
   virtual void operator()(Waveform<float>&,  int)  const = 0;
   virtual void operator()(Waveform<double>&, int)  const = 0;
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
    explicit HighPassFFTFilter(const std::vector<double>&, const std::vector<double>&);

    /**
     *  @brief  Destructor
     */
    ~HighPassFFTFilter();
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(Waveform<float>&,  int)  const override;
    void operator()(Waveform<double>&, int)  const override;

private:
    void highPassFilter(Waveform<float>&, int) const;

    using KernelVec = std::vector<std::complex<float>>;

    std::vector<KernelVec> fHighPassFilterKernels;

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
    explicit LowPassFFTFilter(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~LowPassFFTFilter() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(Waveform<float>&,  int)  const override;
    void operator()(Waveform<double>&, int)  const override;

private:
    void lowPassFilter(Waveform<float>&, int) const;

    unsigned int fStructuringElement;

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
    explicit WindowFFTFilter(const std::vector<std::pair<double,double>>&, const std::vector<std::pair<double,double>>&);

    /**
     *  @brief  Destructor
     */
    ~WindowFFTFilter();
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(Waveform<float>&,  int)  const override;
    void operator()(Waveform<double>&, int)  const override;

private:
    void windowFilter(Waveform<float>&, int) const;

    using KernelVec = std::vector<std::complex<float>>;

    std::vector<KernelVec> fWindowFilterKernels;

    // Keep track of the FFT 
    std::unique_ptr<icarus_signal_processing::ICARUSFFT<float>> fFFT; ///< Object to handle thread safe FFT
};

}

#endif
/** @} */ // end of doxygen group 