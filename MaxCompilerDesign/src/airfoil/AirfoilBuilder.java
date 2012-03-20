package airfoil;


import com.maxeler.maxcompiler.v1.managers.MAX3BoardModel;

public class AirfoilBuilder {

	public static void main(String[] args) {
		AirfoilManager m = new AirfoilManager(MAX3BoardModel.MAX3448A, "AirfoilResCalc");
		m.build();

//
//		Manager sim = new Manager(true, "ResCalcSim", MAX3BoardModel.MAX3448A);
//		Kernel k = new ResCalcKernel(sim.makeKernelParameters("AirfoilResCalc"));
//		sim.setKernel(k);
//		sim.setIO(IOType.ALL_PCIE);
//		sim.build();
	}

}
