// TXMLPGPU network trained with NNO NetworkTrainer at Sat Sep 19 00:33:45 2026
// Input parameters  mom:acos(theta):svt:emc:drc:dch:ifr:ifrExp:ifrAdd
// Output parameters abs(pid)==1
// Training files:
//../Data/PidTuple.root

#include "RhoNNO/TXMLPGPU.h"

double* Recall(double *invec)
{
	static TXMLPGPU net("TXMLPGPU.net");
	float x[7];
	x[0] 	= 1.02292	*	invec[0];	// mom
	x[1] 	= 2.84358	*	invec[1];	// acos(theta)
	x[2] 	= 0.267283	*	invec[2];	// svt
	x[3] 	= 0.763076	*	invec[3];	// emc
	x[4] 	= 1.75257	*	invec[4];	// drc
	x[5] 	= 0.00179238	*	invec[5];	// dch
	x[6] 	= 0.331141	*	invec[6];	// ifr
	return net.Recall(x);
}
