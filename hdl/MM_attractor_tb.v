// ============================================================================
// MM_attractor_tb.v  —  testbench: load an initial state, clock until the state
// register stops changing (fixed-point attractor), report the result.
// Run with: iverilog -o sim MM_attractor.v MM_attractor_tb.v && vvp sim
// ============================================================================
`timescale 1ns/1ps
module tb;
    reg clk=0, load=0;
    reg [21:0] init=0, sa_en=0, sa_val=0;
    reg Isa=0,Dara=0,Btz=0,Dex=0,Len=0,PTEN=1,CRBN=1,BIM_del=0;
    wire [21:0] Q; wire APOP_out, SURVIVAL_out;
    integer i; reg [21:0] prev;

    mm_attractor dut(.clk(clk),.load(load),.init(init),.Isa(Isa),.Dara(Dara),
        .Btz(Btz),.Dex(Dex),.Len(Len),.PTEN(PTEN),.CRBN(CRBN),.BIM_del(BIM_del),
        .sa_en(sa_en),.sa_val(sa_val),.Q(Q),.APOP_out(APOP_out),.SURVIVAL_out(SURVIVAL_out));

    always #5 clk = ~clk;

    task run_to_attractor(input [127:0] name);
        begin
            @(negedge clk); load=1; @(negedge clk); load=0;   // load init
            prev = 22'hxxxxxx;
            for (i=0;i<40 && Q!==prev;i=i+1) begin prev=Q; @(negedge clk); end
            $display("%0s : settled Q=%b  APOP=%b SURVIVAL=%b  (%0d cycles)",
                     name, Q, APOP_out, SURVIVAL_out, i);
        end
    endtask

    initial begin
        // Case 1: ERK-driver fault, untreated  (bit4=ERK pinned to 1)
        sa_en=22'd0; sa_val=22'd0; sa_en[4]=1; sa_val[4]=1;
        Isa=0;Dara=0;Btz=0;Dex=0;Len=0; PTEN=1;CRBN=1;BIM_del=0;
        init=22'd0; init[11]=1; init[12]=1;            // start IRF4=MYC=1
        run_to_attractor("ERK fault, untreated, init=11");
        init=22'd0;                                    // start IRF4=MYC=0
        run_to_attractor("ERK fault, untreated, init=00");

        // Case 2: NF-kB fault (bit10) + VRd (Btz,Dex,Len)
        sa_en=22'd0; sa_val=22'd0; sa_en[10]=1; sa_val[10]=1;
        Btz=1;Dex=1;Len=1; init=22'd0;
        run_to_attractor("NFKB fault + VRd");

        // Case 3: MCL1-high fault (bit19) + Bortezomib
        sa_en=22'd0; sa_val=22'd0; sa_en[19]=1; sa_val[19]=1;
        Btz=1;Dex=0;Len=0; init=22'd0;
        run_to_attractor("MCL1-high fault + Btz");
        $finish;
    end
endmodule
