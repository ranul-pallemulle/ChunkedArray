#ifndef CHUNKEDARRAY_H
#define CHUNKEDARRAY_H

#include "Assertions.hpp"
#include "SharedArray.hpp"
#include <boost/preprocessor.hpp>
#include <cstdlib>
#include <memory>

// ChunkArray
// An efficient mapping between a SharedArray and a more
// cache-friendly, vectorization-friendly datastructure.
// The implementation is a hybrid approach to
// Structure-of-Arrays/Array-of-Structures.

#define ALIGNMENT 32

#ifndef ALIGNMENT
static_assert(false, "Data alignment unspecified for aligned memory allocation "
                     "for chunked array.");
#endif // ALIGNMENT

#ifndef CHUNKTYPE
#define CHUNKTYPE NekDouble
#endif // CHUNKTYPE

#if defined(__INTEL_COMPILER)
#define PRAGMA_VECTORIZE_IVDEP _Pragma("ivdep")
#define PRAGMA_VECTORIZE_ALIGNED _Pragma("vector aligned")
#elif defined(__clang__)
#define PRAGMA_VECTORIZE_IVDEP _Pragma("clang loop vectorize(enable)")
#define PRAGMA_VECTORIZE_ALIGNED
#elif defined(__GNUC__)
#define PRAGMA_VECTORIZE_IVDEP _Pragma("GCC ivdep")
#define PRAGMA_VECTORIZE_ALIGNED
#elif defined(_MSC_VER)
#define PRAGMA_VECTORIZE_IVDEP __pragma("loop (ivdep)")
#define PRAGMA_VECTORIZE_ALIGNED
#endif

#ifndef NUM_IN_OUT_VARS
static_assert(false, "Number of input and output variables for cell model not "
                     "defined. Cannot create ChunkUnit data structure.");
#endif // NUM_IN_OUT_VARS

#ifdef HAS_GATES
#ifndef NUM_GATE_VARS
static_assert(false, "Number of gating variables for cell model not defined. "
                     "Cannot create ChunkUnit data structure.");
#endif // NUM_GATE_VARS
#endif // HAS_GATES

#ifdef HAS_CONCENTRATIONS
#ifdef HAS_GATES
#define CONC_START_IDX BOOST_PP_ADD(NUM_GATE_VARS, 1)
#else
#define CONC_START_IDX 1
#endif // HAS_GATES
#define CONC_END_IDX BOOST_PP_SUB(NUM_IN_OUT_VARS, 1)
#endif // HAS_CONCENTRATIONS

#define APPEND_TEMP(x, y) x##y
#define APPEND(x, y) APPEND_TEMP(x, y)

#define GEN_INVAR(num) T in##num[chunk_len + pad];
#define GEN_OUTVAR(num) T out##num[chunk_len + pad];
#define GEN_GATEVAR(num) T gate##num[chunk_len + pad];

template <typename T, int chunk_len, int pad> class ChunkArray
{
public:
    // Unit of data is a set of short arrays
    struct ChunkUnit
    {
        // number of elements
        size_t size{0};
        // generate input variables
#define BOOST_PP_LOCAL_MACRO(n) GEN_INVAR(n)
#define BOOST_PP_LOCAL_LIMITS (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
        // generate output variables
#define BOOST_PP_LOCAL_MACRO(n) GEN_OUTVAR(n)
#define BOOST_PP_LOCAL_LIMITS (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()

#ifdef HAS_GATES
        // generate gating variables
#define BOOST_PP_LOCAL_MACRO(n) GEN_GATEVAR(n)
#define BOOST_PP_LOCAL_LIMITS (0, NUM_GATE_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
#endif // HAS_GATES
    };

    ChunkArray() = default;
    ChunkArray(size_t num_chunks);
    ChunkArray(size_t num_elements,
               std::initializer_list<double> initial_values);
    ChunkArray(const Array<OneD, const Array<OneD, double>> &inarray);
    ChunkArray(const ChunkArray &rhs);
    ChunkArray(ChunkArray &&rhs);
    ~ChunkArray();
    ChunkArray &operator=(const ChunkArray &rhs);
    ChunkArray &operator=(ChunkArray &&rhs);
    ChunkUnit &operator[](unsigned int i) const
    {
        return m_data[i];
    }
    ChunkUnit *begin()
    {
        return m_data;
    }
    ChunkUnit *end()
    {
        return m_data + m_num_chunks;
    }
    static void fromInArray(const Array<OneD, double> &inarray,
                            ChunkUnit &chunk, int chunk_num);
    static void toOutArray(const ChunkUnit &chunk, Array<OneD, double> &inarray,
                           int chunk_num);
    size_t num_elements()
    {
        return m_num_chunks;
    }

private:
    void allocate_aligned_memory(size_t size, size_t alignment);
    void release_aligned_memory();
    size_t m_num_chunks{0};
    size_t m_remainder{0};
    ChunkUnit *m_data{nullptr}; // heap allocated chunks
};

template <typename T, int chunk_len, int pad>
ChunkArray<T, chunk_len, pad>::ChunkArray(size_t num_chunks)
{
    m_num_chunks = num_chunks;
    m_remainder  = 0;
    allocate_aligned_memory(sizeof(ChunkUnit) * num_chunks, ALIGNMENT);
}

// Construct with initial values
template <typename T, int chunk_len, int pad>
ChunkArray<T, chunk_len, pad>::ChunkArray(
    size_t num_elements, std::initializer_list<double> initial_values)
{
    m_remainder = num_elements % chunk_len;
    if (m_remainder)
    {
        m_num_chunks = num_elements / chunk_len + 1;
    }
    else
    {
        m_num_chunks = num_elements / chunk_len;
    }
    allocate_aligned_memory(sizeof(ChunkUnit) * m_num_chunks, ALIGNMENT);
    int n_iter = m_remainder > 0 ? m_num_chunks - 1 : m_num_chunks;

    for (int i = 0; i < n_iter; ++i)
    {
        m_data[i].size = chunk_len;
        PRAGMA_VECTORIZE_IVDEP
        for (int j = 0; j < chunk_len; ++j)
        {
#define BOOST_PP_LOCAL_MACRO(n) m_data[i].in##n[j] = initial_values.begin()[n];
#define BOOST_PP_LOCAL_LIMITS (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
        }
    }

    if (m_remainder)
    {
        m_data[m_num_chunks - 1].size = m_remainder;
        PRAGMA_VECTORIZE_IVDEP
        for (size_t j = 0; j < m_remainder; ++j)
        {
#define BOOST_PP_LOCAL_MACRO(n)                                                \
    m_data[m_num_chunks - 1].in##n[j] = initial_values.begin()[n];
#define BOOST_PP_LOCAL_LIMITS (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
        }
    }
}

// Construct from Shared Array
template <typename T, int chunk_len, int pad>
ChunkArray<T, chunk_len, pad>::ChunkArray(
    const Array<OneD, const Array<OneD, double>> &inarray)
{
    int np      = inarray[0].num_elements();
    m_remainder = np % chunk_len;
    if (m_remainder)
    {
        m_num_chunks = np / chunk_len + 1;
    }
    else
    {
        m_num_chunks = np / chunk_len;
    }
    allocate_aligned_memory(sizeof(ChunkUnit) * m_num_chunks, ALIGNMENT);
    int n_iter = m_remainder > 0 ? m_num_chunks - 1 : m_num_chunks;

    for (int i = 0; i < n_iter; ++i)
    {
        // set chunk size
        m_data[i].size = chunk_len;

        PRAGMA_VECTORIZE_IVDEP
        for (int j = 0; j < chunk_len; ++j)
        {
            // generate assignments to m_data[i].inK[j] for all K
            // e.g: m_data[i].in2[j] = inarray[2][i*chunk_len + j];
#define BOOST_PP_LOCAL_MACRO(n)                                                \
    m_data[i].in##n[j] = inarray[n][i * chunk_len + j];
#define BOOST_PP_LOCAL_LIMITS (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
        }
    }

    if (m_remainder)
    {
        m_data[m_num_chunks - 1].size = m_remainder;
        // remainder loop
        PRAGMA_VECTORIZE_IVDEP
        for (unsigned int j = 0; j < m_remainder; ++j)
        {
            // generate assignments to m_data[i].inK[j] for all K
#define BOOST_PP_LOCAL_MACRO(n)                                                \
    m_data[m_num_chunks - 1].in##n[j] =                                        \
        inarray[n][(m_num_chunks - 1) * chunk_len + j];
#define BOOST_PP_LOCAL_LIMITS (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
        }
    }
}

// Copy constructor
template <typename T, int chunk_len, int pad>
ChunkArray<T, chunk_len, pad>::ChunkArray(const ChunkArray &rhs)
{
    m_num_chunks = rhs.m_num_chunks;
    m_remainder  = rhs.m_remainder;
    allocate_aligned_memory(sizeof(ChunkUnit) * rhs.m_num_chunks, ALIGNMENT);
    std::copy(rhs.m_data, rhs.m_data + rhs.m_num_chunks, m_data);
}

// Move constructor
template <typename T, int chunk_len, int pad>
ChunkArray<T, chunk_len, pad>::ChunkArray(ChunkArray &&rhs)
{
    m_num_chunks = rhs.m_num_chunks;
    m_remainder  = rhs.m_remainder;
    m_data       = rhs.m_data;
    rhs.m_data   = nullptr;
}

// Destructor
template <typename T, int chunk_len, int pad>
ChunkArray<T, chunk_len, pad>::~ChunkArray()
{
    release_aligned_memory();
    m_num_chunks = 0;
    m_remainder  = 0;
}

// Copy assignment operator
template <typename T, int chunk_len, int pad>
ChunkArray<T, chunk_len, pad> &ChunkArray<T, chunk_len, pad>::operator=(
    const ChunkArray &rhs)
{
    if (this == &rhs)
    {
        return *this;
    }
    if (m_num_chunks != rhs.m_num_chunks)
    {
        release_aligned_memory();
        allocate_aligned_memory(sizeof(ChunkUnit) * rhs.m_num_chunks,
                                ALIGNMENT);
    }
    m_num_chunks = rhs.m_num_chunks;
    m_remainder  = rhs.m_remainder;
    std::copy(rhs.m_data, rhs.m_data + rhs.m_num_chunks, m_data);
    return *this;
}

// Move assignment operator
template <typename T, int chunk_len, int pad>
ChunkArray<T, chunk_len, pad> &ChunkArray<T, chunk_len, pad>::operator=(
    ChunkArray &&rhs)
{
    ASSERTL0((this != &rhs), "Move assignment to self.");
    release_aligned_memory();
    m_num_chunks     = rhs.m_num_chunks;
    m_remainder      = rhs.m_remainder;
    m_data           = rhs.m_data;
    rhs.m_num_chunks = 0;
    rhs.m_remainder  = 0;
    rhs.m_data       = nullptr;
    return *this;
}

// Copy input voltage from SharedArray (single chunk)
template <typename T, int chunk_len, int pad>
void ChunkArray<T, chunk_len, pad>::fromInArray(
    const Array<OneD, double> &inarray, ChunkUnit &chunk, int chunk_num)
{
    PRAGMA_VECTORIZE_IVDEP
#if defined(__INTEL_COMPILER)
#pragma code_align 32
#endif
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
        chunk.in0[i] = inarray[chunk_num * chunk_len + i];
    }
}

// Copy output voltage into SharedArray (single chunk)
template <typename T, int chunk_len, int pad>
void ChunkArray<T, chunk_len, pad>::toOutArray(const ChunkUnit &chunk,
                                               Array<OneD, double> &outarray,
                                               int chunk_num)
{
    PRAGMA_VECTORIZE_IVDEP
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
        outarray[chunk_num * chunk_len + i] = chunk.out0[i];
    }
}

// Aligned memory allocation
template <typename T, int chunk_len, int pad>
void ChunkArray<T, chunk_len, pad>::allocate_aligned_memory(size_t size,
                                                            size_t alignment)
{
#if defined(__INTEL_COMPILER)
    m_data = static_cast<ChunkUnit *>(_mm_malloc(size, alignment));
    if (!m_data)
    {
        ASSERTL0(false, "_mm_malloc failed to allocate aligned memory.");
    }
#elif defined(__GNUC__)
    if (posix_memalign((void **)&m_data, alignment, size))
    {
        ASSERTL0(false, "posix_memalign failed to allocate aligned memory.");
    }
#elif defined(_MSC_VER)
    m_data = static_cast<ChunkUnit *>(_aligned_malloc(size, alignment));
    if (!m_data)
    {
        ASSERTL0(false, "_aligned_malloc failed to allocate aligned memory.");
    }
#else
    ASSERTL0(false, "Unsupported compiler for allocating aligned memory.");
#endif
}

// Aligned memory deallocation
template <typename T, int chunk_len, int pad>
void ChunkArray<T, chunk_len, pad>::release_aligned_memory()
{
#if defined(__INTEL_COMPILER)
    _mm_free(m_data);
#elif defined(__GNUC__)
    free(m_data);
#elif defined(_MSC_VER)
    _aligned_free(m_data);
#else
    ASSERTL0(false, "Unsupported compiler for allocating aligned memory.");
#endif
}

using NekChunkArray = ChunkArray<double, 512, 8>;

#endif /* CHUNKEDARRAY_H */
