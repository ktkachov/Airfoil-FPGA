package airfoil;

import static utils.Utils.array2_t;
import static utils.Utils.array4_t;
import static utils.Utils.arith_t;

import com.maxeler.maxcompiler.v1.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v1.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWVar;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.composite.KArray;

public class BResCalcKernel extends Kernel {

	public BResCalcKernel(KernelParameters params) {
		super(params);
		KArray<HWVar> x1 = io.input("x1", array2_t);
		KArray<HWVar> x2 = io.input("x2", array2_t);
		KArray<HWVar> q1 = io.input("q1", array4_t);
		HWVar adt1 = io.input("adt1", arith_t);
		HWVar bound = io.input("bound", hwBool());
		HWVar dx = x1[0] - x2[0];
		HWVar dy = x1[1] - x2[1];
		HWVar ri = 1.0f / q1[0];

		HWVar gm1 = io.scalarInput("gm1", arith_t);

		HWVar p1 = gm1*(q1[3]-0.5f*ri*(q1[1]*q1[1]+q1[2]*q1[2]));

		KArray<HWVar> res1 = array4_t.newInstance(this);

		KArray<HWVar> qinf = io.scalarInput("qinf", array4_t);

		HWVar vol1 =  ri*(q1[1]*dy - q1[2]*dx);
		HWVar newRi =  1.0f/qinf[0];
		HWVar p2 =  gm1*(qinf[3]-0.5f*ri*(qinf[1]*qinf[1]+qinf[2]*qinf[2]));
		HWVar vol2 =  newRi*(qinf[1]*dy - qinf[2]*dx);
		HWVar eps = io.scalarInput("eps", arith_t);
		HWVar mu = adt1*eps;

		res1[0] <== bound ? 0.0f   : 0.5f*(vol1* q1[0] + vol2* qinf[0]) + mu*(q1[0]-qinf[0]);
		res1[1] <== bound ? p1*dy  : 0.5f*(vol1* q1[1] + p1*dy + vol2* qinf[1] + p2*dy) + mu*(q1[1]-qinf[1]);
		res1[2] <== bound ? -p1*dx : 0.5f*(vol1* q1[2] - p1*dx + vol2* qinf[2] - p2*dx) + mu*(q1[2]-qinf[2]);
		res1[3] <== bound ? 0.0f   : 0.5f*(vol1*(q1[3]+p1) + vol2*(qinf[3]+p2)) + mu*(q1[3]-qinf[3]);

		io.output("res1", res1.getType()) <== res1;
	}

}
