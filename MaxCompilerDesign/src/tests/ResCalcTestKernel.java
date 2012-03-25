package tests;

import com.maxeler.maxcompiler.v1.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v1.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Count;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Count.Counter;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Count.WrapMode;
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


		int lat = 12;


		HWType width_t = hwUInt(width);

		HWVar delay = control.count.simpleCounter(1);

		HWVar num1 = width_t.newInstance(this);
		HWVar half1 = width_t.newInstance(this);
		Count.Params nCountParams1 = control.count.makeParams(width).withMax(num1).withEnable(delay).withWrapMode(WrapMode.COUNT_LT_MAX_THEN_WRAP);
		Counter nCount1 = control.count.makeCounter(nCountParams1);
		HWVar halfPartitionDone1 = (nCount1.getCount() === half1) & half1 !== 0;

		HWVar num2 = width_t.newInstance(this);
		HWVar half2 = width_t.newInstance(this);
		Count.Params nCountParams2 = control.count.makeParams(width).withMax(num2).withEnable(delay).withWrapMode(WrapMode.COUNT_LT_MAX_THEN_WRAP);
		Counter nCount2 = control.count.makeCounter(nCountParams2);
		HWVar halfPartitionDone2 = (nCount2.getCount() === half2) & half2 !== 0;


		HWVar in_en = (total_count === 0) | ((halfPartitionDone1 | halfPartitionDone2));

		HWVar nElems = io.scalarInput("nElems", width_t);
		Count.Params parts_count_params = control.count.makeParams(width).withEnable(in_en);
		Counter parts_count = control.count.makeCounter(parts_count_params);


		HWVar real_en = in_en & (parts_count.getCount() < nElems);
		Count.Params which_sm_params = control.count.makeParams(1).withEnable(real_en);
		Counter which_sm = control.count.makeCounter(which_sm_params);

		HWVar partitionSize = io.input("input1", width_t, real_en);
		HWVar halfPartitionSize = io.input("input2", width_t, real_en);

		int slat = 1;
		HWVar which_sm_offset = total_count < 8 ? 0 : stream.offset(which_sm.getCount(), -8);

		num1 <== which_sm_offset ? stream.offset(partitionSize, -lat) : stream.offset(num1, -slat);
		half1 <== which_sm_offset ? stream.offset(halfPartitionSize, -lat) : stream.offset(half1, -slat);

		num2 <== ~which_sm_offset ? stream.offset(partitionSize, -lat) : stream.offset(num2, -slat);
		half2<== ~which_sm_offset ? stream.offset(halfPartitionSize, -lat) : stream.offset(half2, -slat);

		debug.printf("cycle: %d, input_en: %d, num1: %d, input %d, count1 = %d, half1: %d, num2: %d, count2 = %d, half2: %d, partsCount: %d, which_sm:%d\n",
				total_count, in_en, num1, partitionSize, nCount1.getCount(), half1,  num2, nCount2.getCount(), half2, parts_count.getCount(), which_sm_offset);

		io.output("output", width_t) <== num1;

	}

}
