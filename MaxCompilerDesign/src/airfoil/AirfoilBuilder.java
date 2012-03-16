package airfoil;


import com.maxeler.maxcompiler.v1.managers.MAX3BoardModel;

public class AirfoilBuilder {

	public static void main(String[] args) {
		AirfoilManager m = new AirfoilManager(MAX3BoardModel.MAX3448A, "AirfoilResCalc");
		m.build();
	}

}
