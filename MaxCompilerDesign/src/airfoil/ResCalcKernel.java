/*
 * ResCalc kernel specification
 * */


package airfoil;

import static utils.Utils.array2_t;
import static utils.Utils.array4_t;
import static utils.Utils.float_t;

import com.maxeler.maxcompiler.v1.kernelcompiler.Kernel;
import com.maxeler.maxcompiler.v1.kernelcompiler.KernelParameters;
import com.maxeler.maxcompiler.v1.kernelcompiler.SMIO;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Count;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Mem;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Mem.DualPortMemOutputs;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Mem.RamPortMode;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Mem.RamPortParams;
import com.maxeler.maxcompiler.v1.kernelcompiler.stdlib.core.Mem.RamWriteMode;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWType;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.base.HWVar;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.composite.KArray;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.composite.KStruct;
import com.maxeler.maxcompiler.v1.kernelcompiler.types.composite.KStructType;
import com.maxeler.maxcompiler.v1.utils.MathUtils;


public class ResCalcKernel extends Kernel {

	private final int max_partition_size = 1<<13;
	private final int halo_size = 1<<8;

	private final KStructType node_struct_t
		= new KStructType(
				KStructType.sft("x", array2_t)
			);

	private final KStructType cell_struct_t
		= new KStructType(
				KStructType.sft("q", array4_t),
				KStructType.sft("adt", float_t),
				KStructType.sft("padding", hwUInt(96))
		);

	final int addr_width = MathUtils.bitsToAddress(max_partition_size);
	final HWType addr_t = hwUInt(addr_width);
	final int halo_addr_width = MathUtils.bitsToAddress(halo_size);
	final HWType halo_addr_t = hwUInt(halo_addr_width);


	private final KStructType address_struct_t
		= new KStructType(
				KStructType.sft("node1", addr_t),
				KStructType.sft("node2", addr_t),
				KStructType.sft("cell1", hwUInt(addr_width + 1)),
				KStructType.sft("cell2", hwUInt(addr_width + 1)),
				KStructType.sft("padding", hwUInt(10))
		);

	private final KStructType res_struct_t
		= new KStructType(
				KStructType.sft("res1", array4_t),
				KStructType.sft("res2", array4_t)
			);

	private final KStructType size_struct_t
		= new KStructType(
			KStructType.sft("nodes", addr_t),
			KStructType.sft("cells", addr_t),
			KStructType.sft("edges", addr_t),
			KStructType.sft("halo_cells", addr_t),

			KStructType.sft("iph_cells", addr_t),
			KStructType.sft("iph_nodes", addr_t),
			KStructType.sft("iph_edges", addr_t),

			KStructType.sft("nhd1_cells", addr_t),
			KStructType.sft("nhd1_nodes", addr_t),
			KStructType.sft("nhd1_edges", addr_t),

			KStructType.sft("nhd2_cells", addr_t),
			KStructType.sft("nhd2_nodes", addr_t),
			KStructType.sft("nhd2_edges", addr_t),

			KStructType.sft("padding", hwUInt(23)) //FIXME: REMOVE later!!!
		);


	public ResCalcKernel(KernelParameters params) {
		super(params);

		HWVar total_count = control.count.simpleCounter(48);

		KStruct sizes = size_struct_t.newInstance(this);

		SMIO control_sm = addStateMachine("io_control_sm", new ResControlSM(this, addr_width, 10));
		control_sm.connectInput("cells", (HWVar) sizes.get("cells"));
		control_sm.connectInput("edges", (HWVar) sizes.get("edges"));
		control_sm.connectInput("nodes", (HWVar) sizes.get("nodes"));
		control_sm.connectInput("halo_cells", (HWVar) sizes.get("halo_cells"));
		control_sm.connectInput("nhd1_cells", (HWVar) sizes.get("nhd1_cells"));
		control_sm.connectInput("nhd1_nodes", (HWVar) sizes.get("nhd1_nodes"));
		control_sm.connectInput("nhd1_edges", (HWVar) sizes.get("nhd1_edges"));
		control_sm.connectInput("nhd2_cells", (HWVar) sizes.get("nhd2_cells"));
		control_sm.connectInput("nhd2_nodes", (HWVar) sizes.get("nhd2_nodes"));
		control_sm.connectInput("nhd2_edges", (HWVar) sizes.get("nhd2_edges"));
		control_sm.connectInput("iph_cells", (HWVar) sizes.get("iph_cells"));
		control_sm.connectInput("iph_nodes", (HWVar) sizes.get("iph_nodes"));
		control_sm.connectInput("iph_edges", (HWVar) sizes.get("iph_edges"));

		HWVar read_cell = control_sm.getOutput("read_cell");
		HWVar read_node = control_sm.getOutput("read_node");
		HWVar read_edge = control_sm.getOutput("read_edge");
		HWVar read_sizes = control_sm.getOutput("read_sizes");

		HWVar processing = control_sm.getOutput("processing");
		HWVar output_data = control_sm.getOutput("writing");
		HWVar output_halo = control_sm.getOutput("writing_halo");
		HWVar read_host_halo = control_sm.getOutput("halo_read");


		KStruct zero_sizes = size_struct_t.newInstance(this);
		for (String field : size_struct_t.getFieldNames()) {
			zero_sizes[field] = field.equals("padding") ? hwUInt(23).newInstance(this, 0) : addr_t.newInstance(this, 0);
		}

		KStruct size_input = io.input("sizes", size_struct_t, read_sizes);
		sizes <== total_count < 6 ? zero_sizes : stream.offset(size_input, -6);


//		HWVar nhd1Size = io.scalarInput("nhd1Size", addr_t);
//		HWVar nhd2Size = io.scalarInput("nhd2Size", addr_t);
//		HWVar intraHaloSize = io.scalarInput("intraHaloSize", addr_t);
//		HWVar haloDataSize = io.scalarInput("halo_size", addr_t);

		KStruct node_input_dram = io.input("node_input_dram", node_struct_t, read_cell);
		KStruct cell_input_dram = io.input("cell_input_dram", cell_struct_t, read_cell);
		KStruct address_struct = io.input("addresses", address_struct_t, read_edge);
		KStruct cell_data_host = io.input("input_host", cell_struct_t, read_host_halo);

		HWVar gm1 = io.scalarInput("gm1", float_t);
		HWVar eps = io.scalarInput("eps", float_t);


		HWVar ram_write_counter = control.count.simpleCounter(addr_width, max_partition_size);

		Count.Params hwc_params = control.count.makeParams(halo_addr_width)
			.withEnable(read_host_halo)
			.withMax(halo_size)
			;
		HWVar halo_write_count = control.count.makeCounter(hwc_params).getCount();

		HWVar cell1_addr = address_struct["cell1"];
		HWVar isCell1Halo = cell1_addr > max_partition_size;
		cell1_addr = isCell1Halo ? cell1_addr - max_partition_size : cell1_addr;

		HWVar cell2_addr = address_struct["cell2"];
		HWVar isCell2Halo = cell2_addr > max_partition_size;
		cell2_addr = isCell2Halo ? cell2_addr - max_partition_size : cell2_addr;



		HWVar read_halo_ram = isCell1Halo | isCell2Halo;
		HWVar halo_ram_addr = (read_halo_ram ? (isCell1Halo ? cell1_addr : cell2_addr) : 0).cast(halo_addr_t);

		//RAMs for halo data
		Mem.RamPortParams<KStruct> halo_cell_portA_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_write_count, cell_struct_t)
				.withDataIn(cell_data_host)
				.withWriteEnable(read_host_halo)
				;

		Mem.RamPortParams<KStruct> halo_cell_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, halo_ram_addr, cell_struct_t);

		Mem.DualPortMemOutputs<KStruct> halo_ram_output = mem.ramDualPort(halo_size, RamWriteMode.READ_FIRST, halo_cell_portA_params, halo_cell_portB_params);


		// RAMs for the node data
		Mem.RamPortParams<KStruct> node_ram_portA_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, ram_write_counter, node_struct_t)
				.withDataIn(node_input_dram)
				.withWriteEnable(read_cell)
				;
		HWVar node1_addr = address_struct["node1"];
		RamPortParams<KStruct> node_ram1_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, node1_addr, node_struct_t);
		Mem.DualPortMemOutputs<KStruct> node_ram1_output
			= mem.ramDualPort(max_partition_size, RamWriteMode.WRITE_FIRST, node_ram_portA_params, node_ram1_portB_params);


		HWVar node2_addr = address_struct["node2"];
		RamPortParams<KStruct> node_ram2_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, node2_addr, node_struct_t);
		Mem.DualPortMemOutputs<KStruct> node_ram2_output
			= mem.ramDualPort(max_partition_size, RamWriteMode.WRITE_FIRST, node_ram_portA_params, node_ram2_portB_params);

		//RAMs for cell data
		Mem.RamPortParams<KStruct> cell_ram_portA_params
		= mem.makeRamPortParams(RamPortMode.READ_WRITE, ram_write_counter, cell_struct_t)
			.withDataIn(cell_input_dram)
			.withWriteEnable(read_cell)
			;

		RamPortParams<KStruct> cell_ram1_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, cell1_addr.cast(addr_t), cell_struct_t);
		Mem.DualPortMemOutputs<KStruct> cell_ram1_output
			= mem.ramDualPort(max_partition_size, RamWriteMode.WRITE_FIRST, cell_ram_portA_params, cell_ram1_portB_params);

		RamPortParams<KStruct> cell_ram2_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, cell2_addr.cast(addr_t), cell_struct_t);
		Mem.DualPortMemOutputs<KStruct> cell_ram2_output
			= mem.ramDualPort(max_partition_size, RamWriteMode.WRITE_FIRST, cell_ram_portA_params, cell_ram2_portB_params);


		KStruct cell1 = isCell1Halo ? halo_ram_output.getOutputB() : cell_ram1_output.getOutputB();
		KStruct cell2 = isCell2Halo ? halo_ram_output.getOutputB() : cell_ram2_output.getOutputB();

		//The arithmetic pipeline
		KStruct res_increments = doResMath(
					node_ram1_output.getOutputB(),
					node_ram2_output.getOutputB(),
					cell1,
					cell2,
					eps,
					gm1
				);


		KArray<HWVar> previous_res_value_cell1 = array4_t.newInstance(this);
		KArray<HWVar> previous_res_value_cell2 = array4_t.newInstance(this);

		KArray<HWVar> new_res_value_cell1 = array4_t.newInstance(this);
		KArray<HWVar> new_res_value_cell2 = array4_t.newInstance(this);
		KArray<HWVar> zeroes = array4_t.newInstance(this);

		for (int i = 0; i < array4_t.getSize(); ++i) {
			KArray<HWVar> res1 = res_increments.get("res1");
			new_res_value_cell1[i] <== previous_res_value_cell1[i] + res1[i];

			KArray<HWVar> res2 = res_increments.get("res2");
			new_res_value_cell2[i] <== previous_res_value_cell2[i] + res2[i];

			zeroes[i] <== float_t.newInstance(this, 0.0);

		}


		//RAMs for halo res data

		Count.Params hroc_params = control.count.makeParams(halo_addr_width)
				.withEnable(output_halo)
			;
		HWVar halo_ram_output_count = control.count.makeCounter(hroc_params).getCount();

		KArray<HWVar> halo_res_ram_input = read_halo_ram ? (isCell1Halo ? new_res_value_cell1 : new_res_value_cell2) : zeroes;
		Mem.RamPortParams<KArray<HWVar>> halo_res_ram_portA_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_ram_addr, array4_t)
				.withWriteEnable(read_halo_ram)
				.withDataIn(halo_res_ram_input)
				;

		Mem.RamPortParams<KArray<HWVar>> halo_res_ram_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_ram_output_count, array4_t)
				.withDataIn(zeroes)
				.withWriteEnable(output_halo)
				;

		Mem.DualPortMemOutputs<KArray<HWVar>> halo_res_ram_output
			= mem.ramDualPort(halo_size, RamWriteMode.READ_FIRST, halo_res_ram_portA_params, halo_res_ram_portB_params);


		//RAMs for partition res data

		Count.Params res_output_count_params = control.count.makeParams(addr_width)
			.withEnable(output_data)
			;
		HWVar res_output_count = control.count.makeCounter(res_output_count_params).getCount();


		Mem.RamPortParams<KArray<HWVar>> res_ram1_portA_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, cell1_addr.cast(addr_t), array4_t)
				.withDataIn(new_res_value_cell1)
				.withWriteEnable(~isCell1Halo)
			;

		Mem.RamPortParams<KArray<HWVar>> res_ram1_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, res_output_count, array4_t)
				.withDataIn(zeroes)
				.withWriteEnable(output_data)
				;
		DualPortMemOutputs<KArray<HWVar>> res_ram1_output = mem.ramDualPort(max_partition_size, RamWriteMode.READ_FIRST, res_ram1_portA_params, res_ram1_portB_params);


		Mem.RamPortParams<KArray<HWVar>> res_ram2_portA_params
		= mem.makeRamPortParams(RamPortMode.READ_WRITE, cell2_addr.cast(addr_t), array4_t)
			.withDataIn(new_res_value_cell2)
			.withWriteEnable(~isCell2Halo)
		;
		Mem.RamPortParams<KArray<HWVar>> res_ram2_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, res_output_count, array4_t)
				.withDataIn(zeroes)
				.withWriteEnable(output_data)
				;
		DualPortMemOutputs<KArray<HWVar>> res_ram2_output = mem.ramDualPort(max_partition_size, RamWriteMode.READ_FIRST, res_ram2_portA_params, res_ram2_portB_params);



		//Connect stream offsets to create the loops
		//FIXME: Why 17?
		previous_res_value_cell1 <== stream.offset(isCell1Halo ? halo_res_ram_output.getOutputA() : res_ram1_output.getOutputA(), -17);
		previous_res_value_cell2 <== stream.offset(isCell2Halo ? halo_res_ram_output.getOutputA() : res_ram2_output.getOutputA(), -17);


		KArray<HWVar> res_output = array4_t.newInstance(this);
		for (int i = 0; i < res_output.getSize(); ++i) {
			res_output[i] <== res_ram1_output.getOutputB()[i] + res_ram2_output.getOutputB()[i];
		}

		io.output("result_dram", res_output.getType(), output_data) <== res_output;
		io.output("result_pcie", res_output.getType(), output_halo) <== halo_res_ram_output.getOutputB();
	}


	// The math that produces the res1 and res2 vectors
	private KStruct doResMath(KStruct node1, KStruct node2, KStruct cell1, KStruct cell2, HWVar eps, HWVar gm1){

		KArray<HWVar> x1 = node1["x"];
		KArray<HWVar> x2 = node2["x"];
		KArray<HWVar> q1 = cell1["q"];
		KArray<HWVar> q2 = cell2["q"];
		HWVar adt1 = cell1["adt"];
		HWVar adt2 = cell2["adt"];
		HWVar mu = 0.5f*(adt1+adt2)*eps;

		HWVar dx = x1[0] - x2[0];
		HWVar dy = x1[1] - x2[1];
		HWVar ri = 1.0f / q1[0];
		HWVar p1 = gm1 * (q1[3] - 0.5f*ri*( q1[1] * q1[1] + q1[2] * q1[2]) );
		HWVar vol1 = ri * (q1[1]*dy - q1[2]*dx);

		ri = 1.0f / q1[0];
		HWVar p2 = gm1*(q2[3]-0.5f*ri*(q2[1]*q2[1]+q2[2]*q2[2]));
		HWVar vol2 = ri*(q2[1]*dy - q2[2]*dx);

		KStruct result = res_struct_t.newInstance(this);
		KArray<HWVar> res1 = result["res1"];
		KArray<HWVar> res2 = result["res2"];

		HWVar f = 0.5f*(vol1* q1[0] + vol2* q2[0]) + mu*(q1[0]-q2[0]);
		res1[0] <==  f;
		res2[0] <== -f;

		f = 0.5f*(vol1* q1[1] + p1*dy + vol2* q2[1] + p2*dy) + mu*(q1[1]-q2[1]);
		res1[1] <==  f;
		res2[1] <== -f;

		f = 0.5f*(vol1* q1[2] - p1*dx + vol2* q2[2] - p2*dx) + mu*(q1[2]-q2[2]);
		res1[2] <==  f;
		res2[2] <== -f;

		f = 0.5f*(vol1*(q1[3]+p1)     + vol2*(q2[3]+p2)    ) + mu*(q1[3]-q2[3]);
		res1[3] <==  f;
		res2[3] <== -f;

		return result;

	}

}
