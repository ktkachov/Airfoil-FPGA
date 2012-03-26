package tests;

import com.maxeler.maxcompiler.v1.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v1.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v1.kernelcompiler.SMIO;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Count;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Count.Counter;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWType;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWVar;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.composite.KStructType;
import com.maxeler.maxcompiler.v1.utils.MathUtils;

public class ResCalcTestKernel extends Kernel {
	private final int max_partition_size = 1<<14;
	private final int halo_size = max_partition_size / 10;
	private final int arithmeticPipelineLatency = 6;

	private final int size_width = MathUtils.bitsToAddress(max_partition_size);
	private final HWType size_width_t = hwUInt(size_width);
	private final int sizes_padding = 74;
	private final KStructType size_struct_t
		= new KStructType(
			KStructType.sft("nodes", size_width_t),
			KStructType.sft("cells", size_width_t),
			KStructType.sft("edges", size_width_t),
			KStructType.sft("halo_cells", size_width_t),
			KStructType.sft("halo_nodes", size_width_t),

			KStructType.sft("iph_cells", size_width_t),
			KStructType.sft("iph_nodes", size_width_t),

			KStructType.sft("nhd1_cells", size_width_t),
			KStructType.sft("nhd1_nodes", size_width_t),
			KStructType.sft("nhd1_edges", size_width_t),

			KStructType.sft("nhd2_cells", size_width_t),
			KStructType.sft("nhd2_nodes", size_width_t),
			KStructType.sft("nhd2_edges", size_width_t),

			KStructType.sft("padding", hwUInt(sizes_padding)) //FIXME: REMOVE later!!!
		);

	public ResCalcTestKernel(KernelParameters params) {
		super(params);
		HWVar total_count = control.count.simpleCounter(48);
		int width = 10;


		int lat = 10;


		HWType width_t = hwUInt(width);

		SMIO sm1 = addStateMachine("sm1", new TestStateMachine(this, width));
		SMIO sm2 = addStateMachine("sm2", new TestStateMachine(this, width));

		HWVar en1 = sm1.getOutput("read_sizes");
		HWVar en2 = sm2.getOutput("read_sizes");
		HWVar in_en = (total_count === 0) | (en1 | en2);

		HWVar nElems = io.scalarInput("nElems", width_t);
		Count.Params parts_count_params = control.count.makeParams(width).withEnable(in_en);
		Counter parts_count = control.count.makeCounter(parts_count_params);


		HWVar real_en = in_en & (parts_count.getCount() < nElems);
		Count.Params which_sm_params = control.count.makeParams(1).withEnable(real_en);
		Counter which_sm = control.count.makeCounter(which_sm_params);
		HWVar partitionSize = io.input("input1", width_t, real_en);
		HWVar halfPartitionSize = io.input("input2", width_t, real_en);

		int slat = 1;
//		HWVar which_sm_offset = total_count < 8 ? 0 : stream.offset(which_sm.getCount(), -8);
		sm1.connectInput("size", 		real_en & which_sm.getCount() 	? stream.offset(partitionSize, -lat) 		: stream.offset(sm1.getOutput("size_out"), -slat));
		sm1.connectInput("half_size", 	real_en & which_sm.getCount() 	? stream.offset(halfPartitionSize, -lat) 	: stream.offset(sm1.getOutput("half_size_out"), -slat));
		sm2.connectInput("size", 		real_en & which_sm.getWrap() 	? stream.offset(partitionSize, -lat) 		: stream.offset(sm2.getOutput("size_out"), -slat));
		sm2.connectInput("half_size", 	real_en & which_sm.getWrap() 	? stream.offset(halfPartitionSize, -lat) 	: stream.offset(sm2.getOutput("half_size_out"), -slat));

		debug.printf("cycle: %d, read_size1: %d, size1:%d, read_size2: %d, size2: %d\n", total_count, en1, sm1.getOutput("size_out"), en2, sm2.getOutput("size_out"));

		io.output("output", width_t) <== sm1.getOutput("size_out");

	}

}
