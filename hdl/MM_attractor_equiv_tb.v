// ============================================================================
// MM_attractor_equiv_tb.v  —  self-checking equivalence testbench.
// Reads test vectors written by src/MM_verilog_equiv.py (Python reference) and
// checks that one clock of the Verilog FSM gives the same next state.
// Vector layout (96 bits, MSB first):
//   [95:88] Isa Dara Btz Dex Len PTEN CRBN BIM_del
//   [87:66] state   [65:44] sa_en   [43:22] sa_val   [21:0] expected next state
// Run: python3 src/MM_verilog_equiv.py   (builds and runs this testbench)
// ============================================================================
`timescale 1ns/1ps
`ifndef NVEC
`define NVEC 20000
`endif
module tb_equiv;
    reg clk=0, load=0;
    reg [21:0] init=0, sa_en=0, sa_val=0, expect_q=0;
    reg Isa=0,Dara=0,Btz=0,Dex=0,Len=0,PTEN=1,CRBN=1,BIM_del=0;
    wire [21:0] Q; wire APOP_out, SURVIVAL_out;
    reg [95:0] vec [0:`NVEC-1];
    integer k, mism;

    mm_attractor dut(.clk(clk),.load(load),.init(init),.Isa(Isa),.Dara(Dara),
        .Btz(Btz),.Dex(Dex),.Len(Len),.PTEN(PTEN),.CRBN(CRBN),.BIM_del(BIM_del),
        .sa_en(sa_en),.sa_val(sa_val),.Q(Q),.APOP_out(APOP_out),.SURVIVAL_out(SURVIVAL_out));

    initial begin
        $readmemb("equiv_vectors.txt", vec);
        mism = 0;
        for (k = 0; k < `NVEC; k = k + 1) begin
            {Isa,Dara,Btz,Dex,Len,PTEN,CRBN,BIM_del} = vec[k][95:88];
            init = vec[k][87:66]; sa_en = vec[k][65:44];
            sa_val = vec[k][43:22]; expect_q = vec[k][21:0];
            load = 1; #1 clk = 1; #1 clk = 0;      // load state
            load = 0; #1 clk = 1; #1 clk = 0;      // one transition
            if (Q !== expect_q) begin
                mism = mism + 1;
                if (mism <= 5) $display("mismatch @%0d: got %b expected %b", k, Q, expect_q);
            end
        end
        $display("checked %0d cases, mismatches = %0d", `NVEC, mism);
        $finish;
    end
endmodule
