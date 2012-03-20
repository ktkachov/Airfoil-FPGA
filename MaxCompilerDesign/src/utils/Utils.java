package utils;

import java.util.ArrayList;
import java.util.List;

import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWFix.SignMode;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWType;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWTypeFactory;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWVar;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.composite.KArrayType;

public class Utils {

	public static final HWType arith_t = HWTypeFactory.hwFix(8, 24, SignMode.TWOSCOMPLEMENT);
	public static final KArrayType<HWVar> array2_t = new KArrayType<HWVar>(arith_t, 2);
	public static final KArrayType<HWVar> array4_t = new KArrayType<HWVar>(arith_t, 4);

	public static HWVar adderTree(List<HWVar> elems) {
		if (elems.size() == 1) {
			return elems.get(0);
		} else {
			int midpoint = elems.size()/2;
			List<HWVar> leftSubTree = new ArrayList<HWVar>(elems.subList(0, midpoint));
			List<HWVar> rightSubTree = new ArrayList<HWVar>(elems.subList(midpoint, elems.size()));
			return adderTree(leftSubTree) + adderTree(rightSubTree);
		}
	}

}
