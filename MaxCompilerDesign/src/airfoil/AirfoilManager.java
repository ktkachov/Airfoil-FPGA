package airfoil;

import com.maxeler.maxcompiler.v1.managers.MAXBoardModel;
import com.maxeler.maxcompiler.v1.managers.custom.CustomManager;
import com.maxeler.maxcompiler.v1.managers.custom.Stream;
import com.maxeler.maxcompiler.v1.managers.custom.blocks.KernelBlock;
import com.maxeler.maxcompiler.v1.managers.custom.stdlib.MemoryControlGroup.MemoryAccessPattern;

public class AirfoilManager extends CustomManager {

	public AirfoilManager(MAXBoardModel board_model, String name) {
		super(board_model, name, CustomManager.Target.MAXFILE_FOR_HARDWARE);


		KernelBlock resCalc = addKernel(new ResCalcKernel(makeKernelParameters("ResCalcKernel")));
		Stream in_host = addStreamFromHost("halo_cells");
		Stream nodes_dram = addStreamFromOnCardMemory("nodes_from_fram", MemoryAccessPattern.LINEAR_1D);
		Stream cells_dram = addStreamFromOnCardMemory("cells_from_dram", MemoryAccessPattern.LINEAR_1D);
		Stream addresses_dram = addStreamFromOnCardMemory("addresses_from_dram", MemoryAccessPattern.LINEAR_1D);
		Stream sizes_dram = addStreamFromOnCardMemory("sizes", MemoryAccessPattern.LINEAR_1D);


		resCalc.getInput("node_input_dram") <== nodes_dram;
		resCalc.getInput("cell_input_dram") <== cells_dram;
		resCalc.getInput("addresses") <== addresses_dram;
		resCalc.getInput("input_host") <== in_host;
		resCalc.getInput("sizes") <== sizes_dram;

		Stream to_host = addStreamToHost("res");
		to_host <== resCalc.getOutput("result_pcie");
		Stream to_dram = addStreamToOnCardMemory("to_dram", MemoryAccessPattern.LINEAR_1D);
		to_dram <== resCalc.getOutput("result_dram");


		config.setAllowNonMultipleTransitions(true);//FIXME: Must remove later!!!!!!
	}


}
