package airfoil_backup;

import java.util.ArrayList;
import java.util.List;

import utils.Utils;

import com.maxeler.maxcompiler.v1.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v1.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.KernelMath;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWType;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWVar;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.composite.KArray;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.composite.KArrayType;

public class ADTKernel extends Kernel {

	private final HWType floatType = hwFloat(8, 24);

	private final KArrayType<HWVar> fourArray = new KArrayType<HWVar>(floatType, 4);
	private final KArrayType<HWVar> twoArray = new KArrayType<HWVar>(floatType, 2);

	protected ADTKernel(KernelParameters params) {
		super(params);

		//adt_calc kernel
		KArray<HWVar> q = io.input("q", fourArray);


		KArray<HWVar> x1 = io.input("x1", twoArray);
		KArray<HWVar> x2 = io.input("x2", twoArray);
		KArray<HWVar> x3 = io.input("x3", twoArray);
		KArray<HWVar> x4 = io.input("x4", twoArray);

		HWVar gam = io.scalarInput("gam", floatType);
		HWVar gml = io.scalarInput("gml", floatType);

		HWVar ri = 1.0f / q[0];
		HWVar u = ri * q[1];
		HWVar v = ri * q[2];
		HWVar c = KernelMath.sqrt(gam * gml * (ri * q[3] - 0.5f * (u*u + v*v)));

		List<HWVar> toSum = new ArrayList<HWVar>(4);

		HWVar dx = x2[0] - x1[0];
		HWVar dy = x2[1] - x1[1];
		toSum.add(KernelMath.abs(u*dy - v*dx) + c*KernelMath.sqrt(dx*dx + dy*dy));

		dx = x3[0] - x2[0];
		dy = x3[1] - x2[1];
		toSum.add(KernelMath.abs(u*dy - v*dx) + c*KernelMath.sqrt(dx*dx + dy*dy));

		dx = x4[0] - x3[0];
		dy = x4[1] - x3[1];
		toSum.add(KernelMath.abs(u*dy - v*dx) + c*KernelMath.sqrt(dx*dx + dy*dy));

		dx = x1[0] - x4[0];
		dy = x1[1] - x4[1];
		toSum.add(KernelMath.abs(u*dy - v*dx) + c*KernelMath.sqrt(dx*dx + dy*dy));

		HWVar cfl = io.scalarInput("cfl", floatType);
		HWVar adt = Utils.adderTree(toSum) / cfl;

		io.output("adt", adt.getType()) <== adt;
		//end adt_calc kernel

	}


}
