package airfoil;
import static utils.Utils.array4_t;
import static utils.Utils.arith_t;

import java.util.ArrayList;
import java.util.List;

import utils.Utils;

import com.maxeler.maxcompiler.v1.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v1.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.CounterChain;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWVar;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.composite.KArray;

public class UpdateKernel extends Kernel {

	public UpdateKernel(KernelParameters parameters) {
		super(parameters);
		KArray<HWVar> qold = io.input("qold", array4_t);
		KArray<HWVar> res = io.input("res", array4_t);

		HWVar adt = io.input("adt", arith_t);
		HWVar adti = 1.0f / adt;
		KArray<HWVar> q = array4_t.newInstance(this);
		List<HWVar> rmsToSum = new ArrayList<HWVar>(4);
		for (int i = 0; i < qold.getSize(); ++i) {
			HWVar del = adti * res[i];
			q[i] <== qold[i] - del;
			rmsToSum.add(del * del);
		}

		final int loopLength = 13;
		CounterChain chain = control.count.makeCounterChain();
		HWVar x = chain.addCounter(32, 1);
		HWVar loopCounter = chain.addCounter(loopLength, 1);


		HWVar rmsSum = Utils.adderTree(rmsToSum);
		HWVar count = control.count.simpleCounter(32);
		HWVar carriedSum = arith_t.newInstance(this);
		HWVar rms = count.eq(0) ? 0.0 : carriedSum;
		HWVar newSum = rms + rmsSum;
		carriedSum <== stream.offset(newSum, -loopLength);
		io.scalarOutput("rms", rms.getType()) <== rms;
		io.output("q", q.getType()) <== q;
	}

}
