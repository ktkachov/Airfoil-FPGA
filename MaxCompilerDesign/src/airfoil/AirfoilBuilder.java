package airfoil;


import com.maxeler.maxcompiler.v1.managers.MAX3BoardModel;
import com.maxeler.maxcompiler.v1.managers.custom.CustomManager;

public class AirfoilBuilder {

	public static void main(String[] args) {
		CustomManager sim_manager = new AirfoilSimManager(MAX3BoardModel.MAX3424A, "ResCalcSim", CustomManager.Target.MAXFILE_FOR_SIMULATION);
		sim_manager.build();
//		CustomManager hw_manager = new AirfoilManager(MAX3BoardModel.MAX3424A, "AirfoilResCalc", CustomManager.Target.MAXFILE_FOR_HARDWARE);
//		hw_manager.build();


	}

}
