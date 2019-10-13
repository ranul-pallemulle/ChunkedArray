#include <cmath>
#include "SharedArray.hpp"
#include "ChunkedArray.hpp"

#ifndef NekDouble
#define NekDouble double
#endif

extern NekDouble lastTime;
extern unsigned int substeps;
extern void v_Update(NekChunkArray::ChunkUnit& chunk);
extern NekChunkArray m_data;

void TimeIntegrate(const Array<OneD, const Array<OneD, NekDouble> >& inarray,
		         Array<OneD,       Array<OneD, NekDouble> >& outarray,
		   const NekDouble time)
{
    NekDouble delta_t = (time - lastTime)/substeps;

    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (unsigned int k = 0; k < m_data.num_elements(); ++k)
    {
	NekChunkArray::ChunkUnit& chunk = m_data[k];

	// Copy in voltage
	NekChunkArray::fromInArray(inarray[0], chunk, k);
	// Substep
	for (unsigned int i = 0; i < substeps - 1; ++i)
	{
	    v_Update(chunk);

	    PRAGMA_VECTORIZE_IVDEP
	    PRAGMA_VECTORIZE_ALIGNED
	    for (unsigned int j = 0; j < chunk.size; ++j)
	    {
		// Forward Euler on voltage
		chunk.in0[j] = delta_t * chunk.out0[j] + chunk.in0[j];

		// Forward Euler on concentrations
#ifdef HAS_CONCENTRATIONS
#define BOOST_PP_LOCAL_MACRO(n)						\
		chunk.in##n[j] = delta_t * chunk.out##n[j] + chunk.in##n[j];
#define BOOST_PP_LOCAL_LIMITS (CONC_START_IDX, CONC_END_IDX)
#include BOOST_PP_LOCAL_ITERATE()
#endif
		// Rush-Larsen integration on gating variables
#ifdef HAS_GATES		
#define BOOST_PP_LOCAL_MACRO(n)					\
		chunk.gate##n[j] = -delta_t / chunk.gate##n[j]; \
								\
		chunk.gate##n[j] = exp(chunk.gate##n[j]);	\
								\
		APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) =	\
		    APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) -	\
		    APPEND(chunk.out, BOOST_PP_ADD(n,1)[j]);	\
								\
		APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) =	\
		    APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) *	\
		    APPEND(chunk.gate, n[j]) +			\
		    APPEND(chunk.out, BOOST_PP_ADD(n,1)[j]);
		
#define BOOST_PP_LOCAL_LIMITS (0, NUM_GATE_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
#endif		
	    }
	}
	
	v_Update(chunk);

	PRAGMA_VECTORIZE_IVDEP
	PRAGMA_VECTORIZE_ALIGNED
	for (unsigned int j = 0; j < chunk.size; ++j)
	{
	    // Forward Euler on concentrations
#ifdef HAS_CONCENTRATIONS
#define BOOST_PP_LOCAL_MACRO(n)						\
	    chunk.in##n[j] = delta_t * chunk.out##n[j] + chunk.in##n[j];
#define BOOST_PP_LOCAL_LIMITS (CONC_START_IDX, CONC_END_IDX)
#include BOOST_PP_LOCAL_ITERATE()
#endif
	    // Rush-Larsen integration on gating variables
#ifdef HAS_GATES	    
#define BOOST_PP_LOCAL_MACRO(n)					\
	    chunk.gate##n[j] = -delta_t / chunk.gate##n[j];	\
								\
	    chunk.gate##n[j] = exp(chunk.gate##n[j]);		\
								\
	    APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) =		\
		APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) -	\
		APPEND(chunk.out, BOOST_PP_ADD(n,1)[j]);	\
								\
	    APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) =		\
		APPEND(chunk.in, BOOST_PP_ADD(n,1)[j]) *	\
		APPEND(chunk.gate, n[j]) +			\
		APPEND(chunk.out, BOOST_PP_ADD(n,1)[j]);
		
#define BOOST_PP_LOCAL_LIMITS (0, NUM_GATE_VARS - 1)
#include BOOST_PP_LOCAL_ITERATE()
#endif	    
	}
	// Copy out new voltage
	NekChunkArray::toOutArray(chunk, outarray[0], k);
    }
    lastTime = time;
}
