package tests;

import com.maxeler.maxcompiler.v1.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v1.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v1.kernelcompiler.SMIO;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Count;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Count.Counter;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWType;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWVar;

public class ResCalcTestKernel extends Kernel {
//	private final int max_partition_size = 1<<14;
//	private final int halo_size = max_partition_size / 10;
//	private final int arithmeticPipelineLatency = 6;
//
//	private final int size_width = MathUtils.bitsToAddress(max_partition_size);
//	private final HWType size_width_t = hwUInt(size_width);
//	private final int sizes_padding = 74;
//	private final KStructType size_struct_t
//		= new KStructType(
//			KStructType.sft("nodes", size_width_t),
//			KStructType.sft("cells", size_width_t),
//			KStructType.sft("edges", size_width_t),
//			KStructType.sft("halo_cells", size_width_t),
//			KStructType.sft("halo_nodes", size_width_t),
//
//			KStructType.sft("iph_cells", size_width_t),
//			KStructType.sft("iph_nodes", size_width_t),
//
//			KStructType.sft("nhd1_cells", size_width_t),
//			KStructType.sft("nhd1_nodes", size_width_t),
//			KStructType.sft("nhd1_edges", size_width_t),
//
//			KStructType.sft("nhd2_cells", size_width_t),
//			KStructType.sft("nhd2_nodes", size_width_t),
//			KStructType.sft("nhd2_edges", size_width_t),
//
//			KStructType.sft("padding", hwUInt(sizes_padding)) //FIXME: REMOVE later!!!
//		);

	public ResCalcTestKernel(KernelParameters params) {
		super(params);
		HWVar total_count = control.count.simpleCounter(48);
		int width = 10;


		int lat = 5;


		HWType width_t = hwUInt(width);

		SMIO sm = addStateMachine("ResSM", new TestStateMachine(this, width, lat));
		HWVar in_en = sm.getOutput("read_sizes");

		HWVar nElems = io.scalarInput("nElems", width_t);
//		Count.Params parts_count_params = control.count.makeParams(width).withEnable(in_en).withWrapMode(WrapMode.COUNT_LT_MAX_THEN_WRAP);
//		Counter parts_count = control.count.makeCounter(parts_count_params);

		HWVar real_en = in_en;

		Count.Params which_sm_params = control.count.makeParams(1).withEnable(real_en);
		Counter which_sm = control.count.makeCounter(which_sm_params);

		HWVar partitionSize 	= io.input("input1", width_t, real_en);
		HWVar halfPartitionSize = io.input("input2", width_t, real_en);

		sm.connectInput("size", stream.offset(partitionSize, -lat));
		sm.connectInput("half_size", stream.offset(halfPartitionSize, -lat));

//		debug.printf("cycle: %d, real_en:%d read_size1: %d, size1:%d, read_size2: %d, size2: %d, size:%d, half_size:%d, which_sm:%d\n",
//				total_count,real_en, en1, sm1.getOutput("size_out"), en2, sm2.getOutput("size_out"),
//				partitionSize, halfPartitionSize, which_sm.getCount()
//				);

		io.output("output", width_t) <== partitionSize;

	}

}
