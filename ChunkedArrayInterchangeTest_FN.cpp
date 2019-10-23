#include <vector>
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants_FN.hpp"
#include "ChunkedArray.hpp"

NekChunkArray m_data;

extern NekDouble lastTime;
extern unsigned int substeps;

void init_test(int n)
{
    int nq = n;
    m_data = NekChunkArray(nq, {0.0, 0.0});
}

void v_Update(NekChunkArray::ChunkUnit& chunk)
{
    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& u_old = chunk.in0[i];
	const NekDouble& v_old = chunk.in1[i];
	NekDouble& u_new = chunk.out0[i];
	NekDouble& v_new = chunk.out1[i];

	u_new = u_old - (1.0/3.0)*u_old*u_old*u_old;
	u_new = v_old - u_new;
	u_new = (-1.0/m_epsilon)*u_new;

	v_new = (-1.0*m_gamma)*v_old + u_old;
	v_new = v_new + m_beta;
	v_new = m_epsilon*v_new;
    }
}
