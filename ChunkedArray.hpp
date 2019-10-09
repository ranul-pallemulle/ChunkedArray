#ifndef CHUNKEDARRAY_H
#define CHUNKEDARRAY_H

#include <iostream>
#include <cstdlib>
#include <string>
#include <memory>
#include <boost/preprocessor.hpp>

#include "Assertions.hpp"
#include "SharedArray.hpp"

// ChunkArray
// An efficient mapping between a SharedArray and a more
// cache-friendly, vectorization-friendly datastructure.
// The implementation is a hybrid approach to
// Structure-of-Arrays/Array-of-Structures.

#define ALIGNMENT 32

#ifndef CHUNKTYPE
#define CHUNKTYPE NekDouble
#endif	// CHUNKTYPE

// #define STRINGIFY(a) #a
// #define PRAGMA_VECTORIZE_PRIVATE(a,b) _Pragma(STRINGIFY("omp simd
// private("## a ## "," ## b ##")"))

#if defined(__INTEL_COMPILER)
#define PRAGMA_VECTORIZE_IVDEP _Pragma("ivdep")
#elif defined (__clang__)
#define PRAGMA_VECTORIZE_IVDEP _Pragma("clang loop vectorize(enable)")
#elif defined(__GNUC__)
#define PRAGMA_VECTORIZE_IVDEP _Pragma("GCC ivdep")
#elif defined (_MSC_VER)
#define PRAGMA_VECTORIZE_IVDEP __pragma("loop (ivdep)")
#endif


#ifndef NUM_IN_OUT_VARS
static_assert(false, "Number of input and output variables for cell model not defined. Cannot create ChunkUnit data structure.");
#endif	// NUM_IN_OUT_VARS

#ifndef NUM_GATE_VARS
static_assert(false, "Number of gating variables for cell model not defined. Cannot create ChunkUnit data structure.");
#endif	// NUM_GATE_VARS

#ifdef HAS_CONCENTRATIONS
#define CONC_START_IDX BOOST_PP_ADD(NUM_GATE_VARS, 1)
#define CONC_END_IDX BOOST_PP_SUB(NUM_IN_OUT_VARS, 1)
#endif

#define APPEND_V(x,y) x ## y
#define APPEND(x,y) APPEND_V(x,y)

#define GEN_INVAR(num) DataType in##num[chunk_len + pad];
#define GEN_OUTVAR(num) DataType out##num[chunk_len + pad];
#define GEN_GATEVAR(num) DataType gate##num[chunk_len + pad];

namespace Nektar
{
    template <typename DataType, int chunk_len, int pad>
    class ChunkArray {
    public:
	// Unit of data is a set of short arrays
	struct ChunkUnit {
	    // number of elements
	    size_t size{0};
	    // generate input variables
#define BOOST_PP_LOCAL_MACRO(n)   GEN_INVAR(n)
#define BOOST_PP_LOCAL_LIMITS     (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
	    // generate output variables
#define BOOST_PP_LOCAL_MACRO(n)   GEN_OUTVAR(n)
#define BOOST_PP_LOCAL_LIMITS     (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
	    // generate gating variables
#define BOOST_PP_LOCAL_MACRO(n)   GEN_GATEVAR(n)
#define BOOST_PP_LOCAL_LIMITS     (0, NUM_GATE_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
	};
	
	ChunkArray() = default;
	ChunkArray(size_t num_chunks);
	ChunkArray(const Array<OneD, const Array<OneD, double>>& inarray);
	ChunkArray(const ChunkArray& rhs);
	ChunkArray(ChunkArray&& rhs);
	~ChunkArray();
	ChunkArray& operator=(const ChunkArray& rhs);
	ChunkArray& operator=(ChunkArray&& rhs);
	ChunkUnit& operator[](unsigned int i) const {return m_data[i];}
	ChunkUnit* begin() {return m_data;}
	ChunkUnit* end() {return m_data + m_num_chunks;}
	void fromInArray(const Array<OneD, double>& inarray);
	void toOutArray(Array<OneD, double>& outarray);
	size_t num_elements() {return m_num_chunks;}
	// int get_chunk_length() {return chunk_len;}
    private:
	void allocate_aligned_memory(size_t size, size_t alignment);
	void release_aligned_memory();
	size_t m_num_chunks{0};
	size_t m_remainder{0};
	ChunkUnit* m_data = nullptr; // heap allocated chunks
    };

    template<typename DataType, int chunk_len, int pad>
    ChunkArray<DataType, chunk_len, pad>::ChunkArray(size_t num_chunks)
    {
	m_num_chunks = num_chunks;
	m_remainder = 0;
	allocate_aligned_memory(sizeof(ChunkUnit)*num_chunks, ALIGNMENT);
    }

    // Construct from Shared Array
    template<typename DataType, int chunk_len, int pad>
    ChunkArray<DataType, chunk_len, pad>::ChunkArray(
			  const Array<OneD, const Array<OneD, double>>& inarray)
    {
	int np = inarray[0].num_elements();
	m_remainder = np % chunk_len;
	if (m_remainder)
	{
	    m_num_chunks = np / chunk_len + 1;
	}
	else
	{
	    m_num_chunks = np / chunk_len;
	}
	allocate_aligned_memory(sizeof(ChunkUnit)*m_num_chunks, ALIGNMENT);
	int n_iter = m_remainder > 0 ? m_num_chunks-1 : m_num_chunks;

	for (int i = 0; i < n_iter; ++i)
	{
	    PRAGMA_VECTORIZE_IVDEP
	    for (int j = 0; j < chunk_len; ++j)
	    {
		// set chunk size
		m_data[i].size = chunk_len;
		// generate assignments to m_data[i].inK[j] for all K
		// e.g: m_data[i].in2[j] = inarray[2][i*chunk_len + j];
#define BOOST_PP_LOCAL_MACRO(n) m_data[i].in##n[j] = inarray[n][i*chunk_len + j];
#define BOOST_PP_LOCAL_LIMITS (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
	    }
	}
	// remainder loop
	PRAGMA_VECTORIZE_IVDEP
	for (unsigned int j = 0; j < m_remainder; ++j)
	{
	    m_data[m_num_chunks - 1].size = m_remainder;
	    // generate assignments to m_data[i].inK[j] for all K
#define BOOST_PP_LOCAL_MACRO(n) m_data[m_num_chunks - 1].in##n[j] = \
		inarray[n][(m_num_chunks-1)*chunk_len + j];
#define BOOST_PP_LOCAL_LIMITS (0, NUM_IN_OUT_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
	}
    }

    // Copy constructor
    template<typename DataType, int chunk_len, int pad>
    ChunkArray<DataType, chunk_len, pad>::ChunkArray(const ChunkArray& rhs)
    {
	m_num_chunks = rhs.m_num_chunks;
	m_remainder = rhs.m_remainder;
	allocate_aligned_memory(sizeof(ChunkUnit)*rhs.m_num_chunks, ALIGNMENT);
	std::copy(rhs.m_data, rhs.m_data + rhs.m_num_chunks, m_data);
    }

    // Move constructor
    template<typename DataType, int chunk_len, int pad>
    ChunkArray<DataType, chunk_len, pad>::ChunkArray(ChunkArray&& rhs)
    {
	m_num_chunks = rhs.m_num_chunks;
	m_remainder = rhs.m_remainder;
	m_data = rhs.m_data;

	rhs.m_data = nullptr;
    }

    // Destructor
    template<typename DataType, int chunk_len, int pad>
    ChunkArray<DataType, chunk_len, pad>::~ChunkArray()
    {
	release_aligned_memory();
	m_num_chunks = 0;
	m_remainder = 0;
    }

    // Copy assignment operator
    template<typename DataType, int chunk_len, int pad>
    ChunkArray<DataType, chunk_len, pad>&
    ChunkArray<DataType, chunk_len, pad>::operator=(const ChunkArray& rhs)
    {
	if (this == &rhs) {
	    return *this;
	}
	if (m_num_chunks != rhs.m_num_chunks) {
	    release_aligned_memory();
	    allocate_aligned_memory(sizeof(ChunkUnit)*rhs.m_num_chunks, ALIGNMENT);
	}
	m_num_chunks = rhs.m_num_chunks;
	m_remainder = rhs.m_remainder;
	std::copy(rhs.m_data, rhs.m_data + rhs.m_num_chunks, m_data);
    }

    // Move assignment operator
    template<typename DataType, int chunk_len, int pad>
    ChunkArray<DataType, chunk_len, pad>&
    ChunkArray<DataType, chunk_len, pad>::operator=(ChunkArray&& rhs)
    {
	ASSERTL0((this != &rhs), "Move assignment to self.");
	release_aligned_memory();
	m_num_chunks = rhs.m_num_chunks;
	m_remainder = rhs.m_remainder;
	m_data = rhs.m_data;
	rhs.m_num_chunks = 0;
	rhs.m_remainder = 0;
	rhs.m_data = nullptr;
	return *this;
    }

    // Copy input voltage from SharedArray
    template<typename DataType, int chunk_len, int pad>
    void ChunkArray<DataType, chunk_len, pad>::fromInArray(
			       const Array<OneD, double>& inarray)
    {
	int n_iter = m_remainder > 0 ? m_num_chunks-1 : m_num_chunks;
	for (int i = 0; i < n_iter; ++i)
	{
	    PRAGMA_VECTORIZE_IVDEP
	    for (int j = 0; j < chunk_len; ++j)
	    {
		m_data[i].in0[j] = inarray[i*chunk_len + j];
	    }
	}
	PRAGMA_VECTORIZE_IVDEP
	for (unsigned int j = 0; j < m_remainder; ++j)
	{
	    m_data[m_num_chunks - 1].in0[j] = inarray[(m_num_chunks-1)*chunk_len + j];
	}
    }


    // Copy output voltage into SharedArray
    template<typename DataType, int chunk_len, int pad>
    void ChunkArray<DataType, chunk_len, pad>::toOutArray(
			             Array<OneD, double>& outarray)
    {
	int n_iter = m_remainder > 0 ? m_num_chunks-1 : m_num_chunks;
	for (int i = 0; i < n_iter; ++i)
	{
	    PRAGMA_VECTORIZE_IVDEP
	    for (int j = 0; j < chunk_len; ++j)
	    {
		outarray[i*chunk_len + j] = m_data[i].out0[j];
	    }
	}
	PRAGMA_VECTORIZE_IVDEP
	for (unsigned int j = 0; j < m_remainder; ++j)
	{
	    outarray[(m_num_chunks-1) * chunk_len + j] = \
		m_data[m_num_chunks-1].out0[j];
	}
    }

    // Aligned memory allocation
    template<typename DataType, int chunk_len, int pad>
    void ChunkArray<DataType, chunk_len, pad>::allocate_aligned_memory(size_t size,
								       size_t alignment)
    {
#if defined (__INTEL_COMPILER)
	m_data = static_cast<ChunkUnit*>(_mm_malloc(size, alignment));
	if (!m_data) {
	    ASSERTL0(false, "_mm_malloc failed to allocate aligned memory.");
	}
#elif defined (__GNUC__)
	if (posix_memalign((void**)&m_data, alignment, size)) {
	    ASSERTL0(false, "posix_memalign failed to allocate aligned memory.");
	}
#elif defined (_MSC_VER)
	m_data = static_cast<ChunkUnit*>(_aligned_malloc(size, alignment));
	if (!m_data) {
	    ASSERTL0(false, "_aligned_malloc failed to allocate aligned memory.");
	}
#else
	static_assert(false, "Unsupported compiler for allocating aligned memory.");
#endif
    }
    
    // Aligned memory deallocation
    template<typename DataType, int chunk_len, int pad>
    void ChunkArray<DataType, chunk_len, pad>::release_aligned_memory()
    {
#if defined (__INTEL_COMPILER)
	_mm_free(m_data);
#elif defined (__GNUC__)
	free(m_data);
#elif defined (_MSC_VER)
	_aligned_free(m_data);
#else
	static_assert(false, "Unsupported compiler for allocating aligned memory.");
#endif	
    }
    
}

#endif /* CHUNKEDARRAY_H */
