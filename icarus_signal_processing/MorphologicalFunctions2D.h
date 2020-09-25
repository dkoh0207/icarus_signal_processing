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

#include <vector>
#include <memory>

namespace icarus_signal_processing 
{
/**
 * \class IMorphologicalFunctions2D
 * 
 * \ingroup icarus_signal_processing
 * 
 * \brief Interface class for objects
 */
 
template <class T> using Waveform = std::vector<T>;

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
   virtual void operator()(const Waveform<short>&,  Waveform<short>&)  const = 0;
   virtual void operator()(const Waveform<float>&,  Waveform<float>&)  const = 0;
   virtual void operator()(const Waveform<double>&, Waveform<double>&) const = 0;
};

using FilterFunctionVec = std::vector<std::unique_ptr<IMorphologicalFunctions2D>>; 

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
    explicit Dilation2D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Dilation2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getDilation(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
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
    explicit Erosion2D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Erosion2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getErosion(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
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
    explicit Gradient2D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Gradient2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getGradient(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
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
    explicit Average2D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Average2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getAverage(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
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
    explicit Median2D(const unsigned int structuringElement) : fStructuringElement(structuringElement) {}

    /**
     *  @brief  Destructor
     */
    ~Median2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getMedian(const Waveform<T>&, Waveform<T>&) const;

    unsigned int fStructuringElement;
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
    explicit Opening2D() {}

    /**
     *  @brief  Destructor
     */
    ~Opening2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getOpening(const Waveform<T>&, Waveform<T>&) const;

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
    explicit Closing2D() {}

    /**
     *  @brief  Destructor
     */
    ~Closing2D() {}
 
    /**
    *  @brief Interface functions which provided templated access
    */
    void operator()(const Waveform<short>&,  Waveform<short>&)  const override;
    void operator()(const Waveform<float>&,  Waveform<float>&)  const override;
    void operator()(const Waveform<double>&, Waveform<double>&) const override;

private:
    template <typename T> void getClosing(const Waveform<T>&, Waveform<T>&) const;

};

}

#endif
/** @} */ // end of doxygen group 