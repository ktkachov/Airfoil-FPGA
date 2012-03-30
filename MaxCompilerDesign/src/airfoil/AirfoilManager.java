package airfoil;

import com.maxeler.maxcompiler.v1.managers.MAXBoardModel;
import com.maxeler.maxcompiler.v1.managers.custom.CustomManager;
import com.maxeler.maxcompiler.v1.managers.custom.Stream;
import com.maxeler.maxcompiler.v1.managers.custom.blocks.KernelBlock;
import com.maxeler.maxcompiler.v1.managers.custom.stdlib.MemoryControlGroup.MemoryAccessPattern;

public class AirfoilManager extends CustomManager {

	public AirfoilManager(MAXBoardModel board_model, String name, CustomManager.Target target) {
		super(board_model, name, target);


		KernelBlock resCalc = addKernel(new ResCalcKernel(makeKernelParameters("ResCalcKernel")));
		Stream in_host_cell = addStreamFromHost("halo_cells");
		Stream in_host_node = addStreamFromHost("halo_nodes");

		Stream nodes_dram		= addStreamFromOnCardMemory("nodes_from_fram", MemoryAccessPattern.LINEAR_1D);
		Stream cells_dram 		= addStreamFromOnCardMemory("cells_from_dram", MemoryAccessPattern.LINEAR_1D);
		Stream addresses_dram 	= addStreamFromOnCardMemory("addresses_from_dram", MemoryAccessPattern.LINEAR_1D);
		Stream sizes_dram 		= addStreamFromOnCardMemory("sizes", MemoryAccessPattern.LINEAR_1D);


		resCalc.getInput("node_input_dram") <== nodes_dram;
		resCalc.getInput("cell_input_dram") <== cells_dram;
		resCalc.getInput("addresses") 		<== addresses_dram;
		resCalc.getInput("input_host_cell") <== in_host_cell;
		resCalc.getInput("input_host_node") <== in_host_node;

		resCalc.getInput("sizes") <== sizes_dram;

		Stream to_host = addStreamToHost("res");
		to_host <== resCalc.getOutput("result_pcie");
		Stream to_dram = addStreamToOnCardMemory("to_dram", MemoryAccessPattern.LINEAR_1D);
		to_dram <== resCalc.getOutput("result_dram");

		Stream fromHost = addStreamFromHost("host_to_dram");
		Stream toDram = addStreamToOnCardMemory("write_dram", MemoryAccessPattern.LINEAR_1D);
		toDram <== fromHost;

		Stream toHost 	= addStreamToHost("dram_to_host");
		Stream fromDram = addStreamFromOnCardMemory("read_dram", MemoryAccessPattern.LINEAR_1D);
		toHost <== fromDram;


//		config.setAllowNonMultipleTransitions(true);//FIXME: Must remove later!!!!!!
	}


}
