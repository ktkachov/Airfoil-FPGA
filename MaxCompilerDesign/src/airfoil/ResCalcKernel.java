/*
 * ResCalc kernel specification
 * */


package airfoil;

import static utils.Utils.arith_t;
import static utils.Utils.array2_t;
import static utils.Utils.array4_t;

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

	private final int max_partition_size = 1<<14;
	private final int halo_size = max_partition_size / 10;
	private final int arithmeticPipelineLatency = 6;

	private final int size_width = MathUtils.bitsToAddress(max_partition_size);
	private final HWType size_width_t = hwUInt(size_width);

	private final KStructType node_struct_t
		= new KStructType(
				KStructType.sft("x", array2_t),
				KStructType.sft("padding", hwUInt(64))
			);

	private final int cell_padding = 96;
	private final KStructType cell_struct_t
		= new KStructType(
				KStructType.sft("q", array4_t),
				KStructType.sft("adt", arith_t),
				KStructType.sft("padding", hwUInt(cell_padding))
		);

	final int addr_width = MathUtils.bitsToAddress(max_partition_size);
	final HWType addr_t = hwUInt(addr_width);
	final int halo_addr_width = MathUtils.bitsToAddress(halo_size);
	final HWType halo_addr_t = hwUInt(halo_addr_width);


	private final KStructType address_struct_t
		= new KStructType(
				KStructType.sft("node1", hwUInt(addr_width)),
				KStructType.sft("node2", hwUInt(addr_width)),
				KStructType.sft("cell1", hwUInt(addr_width)),
				KStructType.sft("cell2", hwUInt(addr_width)),
				KStructType.sft("padding", hwUInt(8))
		);

	private final KStructType res_struct_t
		= new KStructType(
				KStructType.sft("res1", array4_t),
				KStructType.sft("res2", array4_t)
			);


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


	public ResCalcKernel(KernelParameters params) {
		super(params);

		System.err.println("addr_t width = " + addr_t);
		System.err.println("cell_struct width = " + cell_struct_t.getTotalBits());
		System.err.println("size_struct width = " + size_struct_t.getTotalBits());
		System.err.println("address_struct type width = " + address_struct_t.getTotalBits());

		HWVar total_count = control.count.simpleCounter(48);
		debug.printf("cycle %d\n", total_count);
		KStruct sizes = size_struct_t.newInstance(this);

		SMIO control_sm = addStateMachine("io_control_sm", new ResControlSM(this, addr_width, 10));
		control_sm.connectInput("cells", (HWVar) sizes.get("cells"));
		control_sm.connectInput("edges", (HWVar) sizes.get("edges"));
		control_sm.connectInput("nodes", (HWVar) sizes.get("nodes"));
		control_sm.connectInput("halo_cells", (HWVar) sizes.get("halo_cells"));
		control_sm.connectInput("halo_nodes", (HWVar) sizes.get("halo_nodes"));
		control_sm.connectInput("nhd1_cells", (HWVar) sizes.get("nhd1_cells"));
		control_sm.connectInput("nhd1_nodes", (HWVar) sizes.get("nhd1_nodes"));
		control_sm.connectInput("nhd1_edges", (HWVar) sizes.get("nhd1_edges"));
		control_sm.connectInput("nhd2_cells", (HWVar) sizes.get("nhd2_cells"));
		control_sm.connectInput("nhd2_nodes", (HWVar) sizes.get("nhd2_nodes"));
		control_sm.connectInput("nhd2_edges", (HWVar) sizes.get("nhd2_edges"));
		control_sm.connectInput("iph_cells", (HWVar) sizes.get("iph_cells"));
		control_sm.connectInput("iph_nodes", (HWVar) sizes.get("iph_nodes"));

		HWVar read_cell = control_sm.getOutput("read_cell");
		HWVar read_node = control_sm.getOutput("read_node");
		HWVar read_edge = control_sm.getOutput("read_edge");
		HWVar read_sizes = control_sm.getOutput("read_sizes");
		debug.printf("read_cell:%d, read_node:%d, read_edge:%d, read_sizes:%d\n", read_cell, read_node, read_edge, read_sizes);

		HWVar processing = control_sm.getOutput("processing");
		HWVar output_data = control_sm.getOutput("writing");
		HWVar output_halo = control_sm.getOutput("writing_halo");
		HWVar read_host_halo_cell = control_sm.getOutput("halo_read_cell");
		HWVar read_host_halo_node = control_sm.getOutput("halo_read_node");


		KStruct zero_sizes = size_struct_t.newInstance(this);
		for (String field : size_struct_t.getFieldNames()) {
			zero_sizes[field] = field.equals("padding") ? hwUInt(sizes_padding).newInstance(this, 0) : size_width_t.newInstance(this, 0);
		}

		KStruct size_input = io.input("sizes", size_struct_t, read_sizes);
//		sizes <== total_count < 6 ? zero_sizes : stream.offset(size_input, -7);
		sizes <== stream.offset(size_input, -6);

		debug.printf("sizes.nodes=%d\n", sizes["nodes"]);

//		HWVar nhd1Size = io.scalarInput("nhd1Size", addr_t);
//		HWVar nhd2Size = io.scalarInput("nhd2Size", addr_t);
//		HWVar intraHaloSize = io.scalarInput("intraHaloSize", addr_t);
//		HWVar haloDataSize = io.scalarInput("halo_size", addr_t);

		KStruct node_input_dram = io.input("node_input_dram", node_struct_t, read_node);
		KStruct cell_input_dram = io.input("cell_input_dram", cell_struct_t, read_cell);
		KStruct edge = io.input("addresses", address_struct_t, read_edge);
		KStruct cell_data_host = io.input("input_host_cell", cell_struct_t, read_host_halo_cell);
		KStruct node_data_host = io.input("input_host_node", node_struct_t, read_host_halo_node);

		HWVar gm1 = io.scalarInput("gm1", arith_t);
		HWVar eps = io.scalarInput("eps", arith_t);


		HWVar ram_write_counter = control.count.simpleCounter(addr_width, max_partition_size);

		Count.Params hwc_params = control.count.makeParams(halo_addr_width)
			.withEnable(read_host_halo_cell)
			.withMax(halo_size)
			;
		HWVar halo_write_count = control.count.makeCounter(hwc_params).getCount();

		HWVar isNoopEdge = isEdgeNoop(edge);

		HWVar cell1_addr = edge["cell1"];
		HWVar isCell1Halo = cell1_addr > (HWVar)sizes["cells"].cast(cell1_addr.getType());
		cell1_addr = isCell1Halo ? cell1_addr - (HWVar)sizes["cells"].cast(cell1_addr.getType()) : cell1_addr;

		HWVar cell2_addr = edge["cell2"];
		HWVar isCell2Halo = cell2_addr > (HWVar)sizes["cells"].cast(cell2_addr.getType());
		cell2_addr = isCell2Halo ? cell2_addr - (HWVar)sizes["cells"].cast(cell2_addr.getType()) : cell2_addr;

		HWVar node1_addr = edge["node1"];
		HWVar isNode1Halo = node1_addr > (HWVar)sizes["nodes"];
		node1_addr = isNode1Halo ? node1_addr - (HWVar)sizes["nodes"] : node1_addr;

		HWVar node2_addr = edge["node2"];
		HWVar isNode2Halo = node2_addr > (HWVar)sizes["nodes"];
		node2_addr = isNode2Halo ? node2_addr - (HWVar)sizes["nodes"] : node2_addr;


		//RAMs for halo data


		HWVar halo_cell_ram1_addr = (isCell1Halo ? cell1_addr : 0).cast(halo_addr_t);
		Mem.RamPortParams<KStruct> halo_cell_ram1_portA_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_write_count, cell_struct_t)
				.withDataIn(cell_data_host)
				.withWriteEnable(read_host_halo_cell)
				;

		Mem.RamPortParams<KStruct> halo_cell_ram1_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, halo_cell_ram1_addr, cell_struct_t);

		Mem.DualPortMemOutputs<KStruct> halo_cell_ram1_output = mem.ramDualPort(halo_size, RamWriteMode.READ_FIRST, halo_cell_ram1_portA_params, halo_cell_ram1_portB_params);


		HWVar halo_cell_ram2_addr = (isCell2Halo ? cell2_addr : 0).cast(halo_addr_t);
		Mem.RamPortParams<KStruct> halo_cell_ram2_portA_params
		= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_write_count, cell_struct_t)
			.withDataIn(cell_data_host)
			.withWriteEnable(read_host_halo_cell)
			;

		Mem.RamPortParams<KStruct> halo_cell_ram2_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, halo_cell_ram2_addr, cell_struct_t);

		Mem.DualPortMemOutputs<KStruct> halo_cell_ram2_output = mem.ramDualPort(halo_size, RamWriteMode.READ_FIRST, halo_cell_ram2_portA_params, halo_cell_ram2_portB_params);



		HWVar halo_node_ram1_addr = (isNode1Halo ? node1_addr : 0).cast(halo_addr_t);
		Mem.RamPortParams<KStruct> halo_node_ram1_portA_params
		= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_write_count, node_data_host.getType())
			.withDataIn(node_data_host)
			.withWriteEnable(read_host_halo_node)
			;

		Mem.RamPortParams<KStruct> halo_node_ram1_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, halo_node_ram1_addr, node_data_host.getType());

		Mem.DualPortMemOutputs<KStruct> halo_node_ram1_output = mem.ramDualPort(halo_size, RamWriteMode.READ_FIRST, halo_node_ram1_portA_params, halo_node_ram1_portB_params);


		HWVar halo_node_ram2_addr = (isNode2Halo ? node2_addr : 0).cast(halo_addr_t);
		Mem.RamPortParams<KStruct> halo_node_ram2_portA_params
		= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_write_count, node_data_host.getType())
			.withDataIn(node_data_host)
			.withWriteEnable(read_host_halo_node)
			;

		Mem.RamPortParams<KStruct> halo_node_ram2_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, halo_node_ram2_addr, node_data_host.getType());

		Mem.DualPortMemOutputs<KStruct> halo_node_ram2_output = mem.ramDualPort(halo_size, RamWriteMode.READ_FIRST, halo_node_ram2_portA_params, halo_node_ram2_portB_params);



		// RAMs for the node data
		Mem.RamPortParams<KStruct> node_ram_portA_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, ram_write_counter, node_struct_t)
				.withDataIn(node_input_dram)
				.withWriteEnable(read_node)
				;
		RamPortParams<KStruct> node_ram1_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_ONLY, node1_addr, node_struct_t);
		Mem.DualPortMemOutputs<KStruct> node_ram1_output
			= mem.ramDualPort(max_partition_size, RamWriteMode.WRITE_FIRST, node_ram_portA_params, node_ram1_portB_params);


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


		KStruct cell1 = isCell1Halo ? halo_cell_ram1_output.getOutputB() : cell_ram1_output.getOutputB();
		KStruct cell2 = isCell2Halo ? halo_cell_ram2_output.getOutputB() : cell_ram2_output.getOutputB();
		KStruct node1 = isNode1Halo ? halo_node_ram1_output.getOutputB() : node_ram1_output.getOutputB();
		KStruct node2 = isNode2Halo ? halo_node_ram2_output.getOutputB() : node_ram2_output.getOutputB();

		cell1 = isNoopEdge ? cell1 : noOpCell();
		cell2 = isNoopEdge ? cell2 : noOpCell();
		node1 = isNoopEdge ? node1 : noOpNode();
		node2 = isNoopEdge ? node2 : noOpNode();

		//The arithmetic pipeline
		KStruct res_increments = doResMath(
					node1,
					node2,
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

			zeroes[i] <== arith_t.newInstance(this, 0.0);

		}


		//RAMs for halo res data

		Count.Params hroc_params = control.count.makeParams(halo_addr_width)
				.withEnable(output_halo)
			;
		HWVar halo_ram_output_count = control.count.makeCounter(hroc_params).getCount();

		KArray<HWVar> halo_res_ram1_input = isCell1Halo ? new_res_value_cell1 : zeroes;
		KArray<HWVar> halo_res_ram2_input = isCell2Halo ? new_res_value_cell2 : zeroes;

		Mem.RamPortParams<KArray<HWVar>> halo_res_ram1_portA_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_cell_ram1_addr, array4_t)
				.withWriteEnable(processing & isCell1Halo)
				.withDataIn(halo_res_ram1_input)
				;

		Mem.RamPortParams<KArray<HWVar>> halo_res_ram1_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_ram_output_count, array4_t)
				.withDataIn(zeroes)
				.withWriteEnable(output_halo)
				;

		Mem.DualPortMemOutputs<KArray<HWVar>> halo_res_ram1_output
			= mem.ramDualPort(halo_size, RamWriteMode.READ_FIRST, halo_res_ram1_portA_params, halo_res_ram1_portB_params);



		Mem.RamPortParams<KArray<HWVar>> halo_res_ram2_portA_params
		= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_cell_ram2_addr, array4_t)
			.withWriteEnable(processing & isCell2Halo)
			.withDataIn(halo_res_ram2_input)
			;

		Mem.RamPortParams<KArray<HWVar>> halo_res_ram2_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, halo_ram_output_count, array4_t)
			.withDataIn(zeroes)
			.withWriteEnable(output_halo)
			;

		Mem.DualPortMemOutputs<KArray<HWVar>> halo_res_ram2_output
			= mem.ramDualPort(halo_size, RamWriteMode.READ_FIRST, halo_res_ram2_portA_params, halo_res_ram2_portB_params);




		//RAMs for partition res data

		Count.Params res_output_count_params = control.count.makeParams(addr_width)
			.withEnable(output_data)
			;
		HWVar res_output_count = control.count.makeCounter(res_output_count_params).getCount();


		Mem.RamPortParams<KArray<HWVar>> res_ram1_portA_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, cell1_addr.cast(addr_t), array4_t)
				.withDataIn(new_res_value_cell1)
				.withWriteEnable(~isCell1Halo & processing)
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
			.withWriteEnable(~isCell2Halo & processing)
		;
		Mem.RamPortParams<KArray<HWVar>> res_ram2_portB_params
			= mem.makeRamPortParams(RamPortMode.READ_WRITE, res_output_count, array4_t)
				.withDataIn(zeroes)
				.withWriteEnable(output_data)
				;
		DualPortMemOutputs<KArray<HWVar>> res_ram2_output = mem.ramDualPort(max_partition_size, RamWriteMode.READ_FIRST, res_ram2_portA_params, res_ram2_portB_params);



		//Connect stream offsets to create the loops
		previous_res_value_cell1 <== stream.offset(isCell1Halo ? halo_res_ram1_output.getOutputA() : res_ram1_output.getOutputA(), -arithmeticPipelineLatency);
		previous_res_value_cell2 <== stream.offset(isCell2Halo ? halo_res_ram2_output.getOutputA() : res_ram2_output.getOutputA(), -arithmeticPipelineLatency);


		KArray<HWVar> res_dram_output = array4_t.newInstance(this);
		KArray<HWVar> res_host_output = array4_t.newInstance(this);

		for (int i = 0; i < res_dram_output.getSize(); ++i) {
			res_dram_output[i] <== res_ram1_output.getOutputB()[i] + res_ram2_output.getOutputB()[i];
			res_host_output[i] <== halo_res_ram1_output.getOutputB()[i] + halo_res_ram2_output.getOutputB()[i];
		}
		debug.printf("----------------------------\n");
		io.output("result_dram", res_dram_output.getType(), output_data) <== res_dram_output;
		io.output("result_pcie", res_dram_output.getType(), output_halo) <== res_host_output;
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


	KStruct noOpNode() {
		KStruct res = node_struct_t.newInstance(this);
		KArray<HWVar> x = array2_t.newInstance(this);
		for (int i = 0; i < x.getSize(); ++i) {
			x[i] <== x[i].getType().newInstance(this, i);
		}
		KArray<HWVar> n = res["x"];
		n <== x;
		HWVar r_pad = res["padding"];
		r_pad <== r_pad.getType().newInstance(this, 0);
		return res;
	}

	@SuppressWarnings("unchecked")
	KStruct noOpCell() {
		KStruct res = cell_struct_t.newInstance(this);
		KArray<HWVar> q = (KArray<HWVar>) res["q"].getType().newInstance(this);

		for (int i = 0; i < q.getSize(); ++i) {
			q[i] <== q[i].getType().newInstance(this, i);
		}
		KArray<HWVar> r_q = res["q"];
		r_q <== q;
		HWVar r_adt = res["adt"];
		r_adt <== r_adt.getType().newInstance(this, 0);
		HWVar r_pad = res["padding"];
		r_pad <== r_pad.getType().newInstance(this, 0);
		return res;
	}

	HWVar isEdgeNoop(KStruct edge) {
		return ((HWVar)edge["node1"] === 0) &
		((HWVar)edge["node2"] === 0) &
		((HWVar)edge["cell1"] === 0) &
		((HWVar)edge["cell2"] === 0)
		;

	}

}
