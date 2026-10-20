// ============================================================================
// MM_attractor.v  —  Multiple-Myeloma Boolean attractor model as a synthesizable
// SEQUENTIAL circuit (finite-state machine).
//   * 22 state nodes -> one 22-bit state register Q (D flip-flops), clocked.
//   * Next-state logic D = combinational AND/OR/NOT (Section 8 netlist).
//   * Mutation = STUCK-AT fault: sa_en[i]=1 pins node i to sa_val[i].
//   * Drugs = input pins; APOP is an absorbing latch (death).
// Logic is a line-by-line translation of MM_attractor._next (validated).
//
// Bit order of Q / D / sa_* (index : node):
//   0:STAT3 1:PI3K 2:AKT 3:MTOR 4:ERK 5:BCAT 6:IGF1R 7:RANK 8:NIK 9:TNF
//   10:NFKB 11:IRF4 12:MYC 13:CCND2 14:MDM2 15:P53 16:P21 17:CDK46 18:E2F
//   19:MCL1 20:BCLXL 21:APOP
// ============================================================================
module mm_attractor (
    input  wire        clk,
    input  wire        load,        // 1: load `init` into the register (set initial condition)
    input  wire [21:0] init,        // initial state (use init[21]=APOP=0 to start alive)
    // drug input pins
    input  wire        Isa, Dara, Btz, Dex, Len,
    // microenvironment flags (1 = functional / present; fault drives these low)
    input  wire        PTEN,        // PTEN loss fault -> PTEN=0
    input  wire        CRBN,        // CRBN loss (lenalidomide resistance) fault -> CRBN=0
    input  wire        BIM_del,     // BCL2L11 (BIM) deletion fault -> 1
    // stuck-at fault control over the 22 state nodes
    input  wire [21:0] sa_en,       // 1 = node is pinned (stuck-at)
    input  wire [21:0] sa_val,      // value the pinned node is held at
    output wire [21:0] Q,           // current state (flip-flop outputs)
    output wire        APOP_out,    // committed to apoptosis (death)
    output wire        SURVIVAL_out // anti-apoptotic and not dead
);
    reg [21:0] q;
    assign Q = q;

    // ---- unpack current state (flip-flop outputs feed back here) ----
    wire STAT3=q[0],PI3K=q[1],AKT=q[2],MTOR=q[3],ERK=q[4],BCAT=q[5],IGF1R=q[6],
         RANK=q[7],NIK=q[8],TNF=q[9],NFKB=q[10],IRF4=q[11],MYC=q[12],CCND2=q[13],
         MDM2=q[14],P53=q[15],P21=q[16],CDK46=q[17],E2F=q[18],MCL1=q[19],
         BCLXL=q[20],APOP=q[21];

    // ---- derived combinational signals ----
    wire Len_active = Len & CRBN;
    wire antiCD38   = Isa | Dara;
    wire NOXA       = Btz;
    wire BIM        = (Dex & ~STAT3 & ~AKT) & ~BIM_del;
    wire free_MCL1  = MCL1  & ~NOXA & ~BIM;
    wire free_BCLXL = BCLXL & ~BIM;
    wire guarded    = free_MCL1 | free_BCLXL;
    wire pro        = NOXA | BIM | antiCD38;

    // ---- next-state logic (D), before stuck-at faults ----
    wire [21:0] Draw;
    assign Draw[0]  = 1'b0;                                   // STAT3' (IL6=0)
    assign Draw[1]  = IGF1R | ~PTEN;                          // PI3K'
    assign Draw[2]  = PI3K;                                   // AKT'
    assign Draw[3]  = AKT;                                    // MTOR'
    assign Draw[4]  = IGF1R;                                  // ERK'
    assign Draw[5]  = 1'b0;                                   // BCAT' (WNT=0)
    assign Draw[6]  = 1'b0;                                   // IGF1R' (IGF1=0)
    assign Draw[7]  = 1'b0;                                   // RANK' (RANKL=0)
    assign Draw[8]  = RANK;                                   // NIK'
    assign Draw[9]  = 1'b0;                                   // TNF' (TNF_in=0)
    assign Draw[10] = (TNF | RANK | AKT | NIK) & ~Btz & ~Dex; // NFKB'
    assign Draw[11] = ~Len_active & (NFKB | MYC);             // IRF4'  (<- MYC feedback)
    assign Draw[12] = IRF4 & (ERK | BCAT | STAT3 | NFKB | MTOR); // MYC' (<- IRF4 feedback)
    assign Draw[13] = MYC | ERK | BCAT | NFKB;               // CCND2'
    assign Draw[14] = AKT | ERK;                             // MDM2'
    assign Draw[15] = ~MDM2;                                 // P53'
    assign Draw[16] = P53 | Len_active;                     // P21'
    assign Draw[17] = CCND2 & ~P21;                         // CDK46'
    assign Draw[18] = CDK46;                                // E2F'
    assign Draw[19] = STAT3 | NFKB | IRF4;                  // MCL1'
    assign Draw[20] = NFKB | STAT3;                         // BCLXL'
    assign Draw[21] = APOP | (pro & ~guarded);              // APOP' (absorbing latch)

    // ---- apply stuck-at faults: D = en ? val : Draw ----
    wire [21:0] D = (Draw & ~sa_en) | (sa_val & sa_en);

    // ---- state register (D flip-flops) ----
    always @(posedge clk) begin
        if (load) q <= (init & ~sa_en) | (sa_val & sa_en);
        else      q <= D;
    end

    // ---- decoded phenotype outputs ----
    wire ANTIAPOP = MCL1 | BCLXL;
    assign APOP_out     = APOP;
    assign SURVIVAL_out = ANTIAPOP & ~APOP;
endmodule
