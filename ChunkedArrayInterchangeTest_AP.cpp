#include <vector>
#include "VmathArray.hpp"
#include "Assertions.hpp"
#include "CellModelConstants_AP.hpp"
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
    NekDouble tmp1, tmp2;

    PRAGMA_VECTORIZE_IVDEP
    PRAGMA_VECTORIZE_ALIGNED
    for (unsigned int i = 0; i < chunk.size; ++i)
    {
	const NekDouble& u_old = chunk.in0[i];
	const NekDouble& v_old = chunk.in1[i];
	NekDouble& u_new = chunk.out0[i];
	NekDouble& v_new = chunk.out1[i];

	// Ru = au
	tmp1 = m_a * u_old;
	// Ru = -(1+a)uu + au
	tmp1 = (-1.0 - m_a) * u_old * u_old + tmp1;
	// Ru = uuu - (1+a)uu + au
	tmp1 = u_old * u_old * u_old + tmp1;
	// Ru = k(uuu - (1+a)uu + au)
	tmp1 = m_k * tmp1;
	// Ru = k(uuu - (1+a)uu + au) + I_stim
	u_new = u_new + tmp1;
	u_new = u_old * v_old + tmp1;
	u_new = -u_new;

	tmp2 = m_mu2 + u_old;
	tmp2 = v_old/tmp2;
	tmp2 = m_mu1*tmp2;
	tmp2 = tmp2 + m_eps;
	tmp1 = (-m_a - 1) + u_old;
	tmp1 = m_k * tmp1;
	tmp1 = u_old * tmp1 + v_old;
	tmp1 = -tmp1;
	v_new = tmp1 * tmp2;
    }
}
