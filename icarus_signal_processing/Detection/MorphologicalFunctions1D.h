/**
 * \file MorphologicalFunctions1D.h
 *
 * \ingroup icarus_signal_processing
 * 
 * \brief Class def header for objects to perform various morphological filter operations
 *
 * @author koh0207
 */

/** \addtogroup icarus_signal_processing

    @{*/
#ifndef __SIGPROC_TOOLS_MORPHOLOGICALFUNCTIONS1D_H__
#define __SIGPROC_TOOLS_MORPHOLOGICALFUNCTIONS1D_H__

#include <vector>
#include <memory>

#include "icarus_signal_processing/ICARUSSigProcDefs.h"

namespace icarus_signal_processing 
{
/**
 * \class IMorphologicalFunctions1D
 * 
 * \ingroup icarus_signal_processing
 * 
 * \brief Interface class for objects
 */
 
template <class T> using Waveform = std::vector<T>;

class IMorphologicalFunctions1D
{
public:
   /**
   *  @brief  Virtual Destructor
   */
   virtual ~IMorphologicalFunctions1D() noexcept = default;

   //@{ 
   /**
   *  @brief Interface functions which provided templated access
   */
   virtual void operator()(const Waveform<bool>&,   Waveform<bool>&)   const = 0;
   virtual void operator()(const Waveform<short>&,  Waveform<short>&)  const = 0;
   virtual void operator()(const Waveform<float>&,  Waveform<float>&)  const = 0;
   virtual void operator()(const Waveform<double>&, Waveform<double>&) const = 0;
   //@}
};

using FilterFunctionVec = std::vector<std::unique_ptr<IMorphologicalFunctions1D>>; 

/**
* \class Dilation1D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the dilation of a waveform
*/
class Dilation1D : virtual public IMorphologicalFunctions1D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Dilation1D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Dilation1D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<bool>&,   Waveform<bool>&)   const override;
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getDilation(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
};

/**
* \class Erosion1D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the erosion of a waveform
*/
class Erosion1D : virtual public IMorphologicalFunctions1D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Erosion1D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Erosion1D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<bool>&,   Waveform<bool>&)   const override;
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getErosion(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
};

/**
* \class Gradient1D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the gradient of a wabeform
*/
class Gradient1D : virtual public IMorphologicalFunctions1D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Gradient1D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Gradient1D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<bool>&,   Waveform<bool>&)   const override;
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getGradient(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
};

/**
* \class Average1D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the average of a waveform
*/
class Average1D : virtual public IMorphologicalFunctions1D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Average1D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Average1D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<bool>&,   Waveform<bool>&)   const override;
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getAverage(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
};

/**
* \class Median1D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the median of a waveform
*/
class Median1D : virtual public IMorphologicalFunctions1D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Median1D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Median1D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<bool>&,   Waveform<bool>&)   const override;
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getMedian(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
};

/**
* \class Opening1D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the opening of a waveform
*/
class Opening1D : virtual public IMorphologicalFunctions1D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Opening1D() {}

    /**
     *  @brief  Destructor
     */
    ~Opening1D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<bool>&,   Waveform<bool>&)   const override;
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getOpening(const Waveform<T>&, Waveform<T>&) const;

};

/**
* \class Closing1D
* 
* \ingroup icarus_signal_processing
* 
* \brief Return the closing of a waveform
*/
class Closing1D : virtual public IMorphologicalFunctions1D
{
public:
    /**
     *  @brief  Constructor
     */
    explicit Closing1D() {}

    /**
     *  @brief  Destructor
     */
    ~Closing1D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<bool>&,   Waveform<bool>&)   const override;
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getClosing(const Waveform<T>&, Waveform<T>&) const;

};

}

#endif
/** @} */ // end of doxygen group 
