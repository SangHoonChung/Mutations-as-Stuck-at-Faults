# 다발성골수종(Multiple Myeloma) 10개 KEGG Pathway 통합 매핑
### Master Crosstalk Edge List + 5개 약물 오버레이

업로드된 10개 맵을 기준으로, **맵끼리 이어지는 노드-노드 연결(crosstalk)** 을 hub 중심으로 정리하고, 그 위에 isatuximab / daratumumab / bortezomib / dexamethasone / lenalidomide의 작용점을 얹었습니다. 통합 figure를 그릴 때 그대로 edge로 쓸 수 있도록 만든 문서입니다.

---

## 0. 표기법 (notation)

- **`[맵ID] 노드`** : 어느 맵의 어떤 노드인지
- `→` 활성화 · `⊣` 억제 · `⇢` 전사/유전자발현(DNA) · `⊗` 분해(proteasomal degradation)
- **★** = 약물 작용점
- 신뢰도 태그: **(KEGG)** = 해당 맵에 실제로 그려진 crosstalk · **(확립)** = 문헌상 확립된 연결 · **(추정)** = 내가 추론한 연결(검증 권장)

업로드 맵 ID 대조:
`hsa04020` Calcium · `hsa05202` Transcriptional misregulation in cancer · `hsa04110` Cell cycle · `hsa04210` Apoptosis · `map03050` Proteasome · `hsa04310` Wnt · `hsa04014` Ras · `hsa04630` JAK-STAT · `hsa04151` PI3K-Akt · `hsa04064` NF-κB

---

## 1. 맵 인벤토리 — 각 맵의 핵심 노드 & 맵이 선언한 crosstalk

각 맵에서 통합에 쓰이는 핵심 노드와, KEGG가 그림 안에 박스로 선언해 둔 외부 pathway 연결만 추렸습니다.

| 맵 | 핵심 노드 (통합 관련) | 맵에 그려진 crosstalk 박스 |
|---|---|---|
| **Calcium** `hsa04020` | **ADPRcyclase(=CD38)** → cADPR/NAADP, RYR, IP3R, TPC/TRPML, Ca²⁺, CALM→CaN→NFAT, RTK→PLCγ | → MAPK, → Apoptosis, → Phosphatidylinositol |
| **Transcr. misreg.** `hsa05202` | **IgH-MAF→CCND2 / ITGB7**, IgH-MMSET→H3→CCND2 (MM 블록), c-Myc, p53→p21/Bax/Bak | → p53 signaling |
| **Cell cycle** `hsa04110` | GSK3β, CycD-CDK4/6, Rb-E2F, p21(Cip1)/p27(Kip1), p53, APC/C, SCF-Skp2 | → MAPK, Ubiquitin mediated proteolysis, → Apoptosis |
| **Apoptosis** `hsa04210` | Bcl-2/Bcl-XL/**Mcl-1** ⊣ Bax/Bak/**Bim**/Bid/Bad/**NOXA**/PUMA, MOMP→CytC→CASP9→CASP3, ER stress(PERK→CHOP), IKK→IκBα→NF-κB | NF-κB, PI3K-Akt, MAPK, Calcium, TNF, p53, Protein processing in ER |
| **Proteasome** `map03050` | 26S = 19S(PA700) + 20S(α1-7, **β5**=bortezomib 표적), immunoproteasome β1i/β2i/β5i | Ubiquitin mediated proteolysis |
| **Wnt** `hsa04310` | Frizzled/LRP→Dvl, 파괴복합체(Axin/APC/**GSK-3β**/CKIα), β-catenin→TCF/LEF⇢c-myc/cycD, β-TrCP⊗, Wnt/Ca²⁺ | → Cell cycle, Ubiquitin mediated proteolysis, p53, TGF-β, MAPK |
| **Ras** `hsa04014` | RTK→Ras, Ras→Raf→MEK→ERK, Ras→PI3K→Akt, Akt→IKK→NFκB, Akt→BAD⊣BCL-X | MAPK, PI3K-Akt, Calcium, Rap1 |
| **JAK-STAT** `hsa04630` | IL-6family→JAK1/2→**STAT3**⇢**MCL1**/Bcl-XL/c-Myc/CycD/p21, SHP2→Ras, PI3K→AKT→mTOR, SOCS⊗ | MAPK, PI3K-AKT, Ubiquitin mediated proteolysis |
| **PI3K-Akt** `hsa04151` | RTK→PI3K→PIP3→PDK1→**AKT** → (GSK3, IKK→NFκB, BAD, FOXO→Bim, Casp9, MDM2⊣p53, mTOR, CCND1/Myc) | mTOR, NFκB, MAPK, FoxO, p53, Cell cycle, Apoptosis, JAK/STAT |
| **NF-κB** `hsa04064` | IKKα/β/NEMO→IκBα **⊗(Proteasome)**→p50/p65⇢생존유전자(Bcl-XL/Bcl-2/XIAP/cIAP), 비정규: NIK→p100/RelB→p52, RANKL→RANK→NIK | **Proteasome**, Ubiquitin mediated proteolysis, Calcium, TCR/BCR, Osteoclast diff. |

---

## 2. 마스터 Crosstalk 엣지 리스트 (Hub 중심)

10개 맵을 실제로 묶어주는 "접점(hub)"은 아래 10개입니다. 각 hub가 어떤 맵의 어떤 노드를 잇는지를 edge로 적었습니다.

### Hub 1 — IκB / NF-κB (정규 생존축)
- `[map03050] 20S proteasome` ⊗ `[hsa04064] IκBα` → NF-κB(p50/p65) 핵 이행 **(KEGG)**
- `[hsa04064] NF-κB` ⇢ `[hsa04210] Bcl-XL / Bcl-2 / XIAP / cIAP / FLIP` (pro-survival 유전자 공유) **(KEGG)**
- `[hsa04151] AKT` → `[hsa04064] IKK` → NF-κB **(KEGG)**
- `[hsa04014] Akt/PKB` → `[hsa04064] IKK` → NF-κB **(KEGG)**
- `[hsa04064] NIK→p52` ⇢ SDF-1α / BAFF / ICAM (비정규, MM 생존+골질환) **(KEGG)**
- ★ **Bortezomib** = `[map03050] β5` ⊗ 차단 → IκBα 안정화 → NF-κB OFF *(단, 구성적 NF-κB는 역설적 활성 가능 — caveat 참고)*
- ★ **Dexamethasone** = NF-κB transrepression (off-map GR → 이 축 억제) **(확립)**

### Hub 2 — Ubiquitin-proteasome (분해 기계)
- `[hsa04110] SCF-Skp2 / APC/C` ⊗ `p21(Cip1) / p27(Kip1) / cyclins` **(KEGG)**
- `[hsa04310] β-TrCP/Skp1/Cul1` ⊗ `β-catenin` **(KEGG)**
- `[hsa04630] SOCS→ECS complex` ⊗ STAT 신호 음성피드백 **(KEGG)**
- 위 세 분해 반응이 모두 `[map03050] 26S proteasome`으로 수렴 **(확립)**
- ★ **Bortezomib** = `[map03050] β5` (전 분해 반응 상류 차단)
- ★ **Lenalidomide** = off-map **CRL4-CRBN E3 ligase**가 IKZF1/3에 폴리유비퀴틴 부착 → 이 proteasome 기계로 투입 **(확립)** *(CRBN/IKZF는 이 맵들에 노드로 없음)*

### Hub 3 — GSK3β (3개 맵 공유 노드)
- `[hsa04151] AKT` ⊣ `GSK3` **(KEGG)**
- `[hsa04310] GSK-3β` = β-catenin 파괴복합체 구성 **(KEGG)**
- `[hsa04110] GSK3β` ⊣ `CycD` **(KEGG)**
- `[hsa04151] GSK3` ⊣ `CCND1 / Myc` 안정성 조절 **(KEGG)** → Hub 6로 연결

### Hub 4 — RTK / Ras / PI3K / MAPK (증식 backbone)
- `[hsa04014] RTK→Ras→Raf→MEK→ERK` = `[hsa04151]`·`[hsa04210]`에도 동일 ERK 축 **(KEGG)**
- `[hsa04014] Ras` → `[hsa04151] PI3K→AKT` **(KEGG)**
- `[hsa04014/hsa04151] RTK` → `[hsa04020] PLCγ → Ca²⁺` **(KEGG)**
- `[hsa04630] SHP2→Ras→Raf` (IL-6 → MAPK) **(KEGG)**
- `[MAPK] → [hsa04110] GSK3β→CycD` (증식 → 세포주기) **(KEGG)**
- 약물 직접 표적 아님. RAS 돌연변이는 MM에서 가장 흔(40–55%) → 약물이 맞서는 증식 입력축으로 표시

### Hub 5 — Bcl-2 family / 미토콘드리아 apoptosis ★수렴 핵심★
- 항-apoptotic 수렴: `[hsa04630] STAT3⇢MCL1`, `[hsa04151] CREB⇢Mcl-1 / AKT⊣BAD`, `[hsa04064] NF-κB⇢Bcl-XL/Bcl-2`, `[hsa04014] Akt→BCL-X` → 모두 `[hsa04210] Mcl-1/Bcl-XL/Bcl-2` **(KEGG)**
- 친-apoptotic 투입: `[hsa04151] FOXO→Bim`, `p53→PUMA/NOXA`, `Akt⊣BAD(인산화)` → `[hsa04210]` **(KEGG)**
- `[hsa04210] Bax/Bak → MOMP → CytC → CASP9 → CASP3` (최종 출력) **(KEGG)**
- ★ **Bortezomib** → `[hsa04210] NOXA` 안정화 ⊣ Mcl-1 → 미토콘드리아 apoptosis **(확립)**
- ★ **Dexamethasone** → `[hsa04210] Bim` 유도 **(확립)**
- ★ **Isatuximab/Daratumumab** → 직접 apoptosis 입력(면역 layer, off-map) → effector caspase **(확립)**

### Hub 6 — c-Myc / Cyclin D / p21 (전사→세포주기)
- c-Myc 수렴: `[hsa04310] β-catenin⇢c-myc`, `[hsa04630] STAT⇢c-Myc`, `[hsa04151] GSK3→Myc`, `[hsa05202] IRF4/MYC` **(KEGG)**
- CycD 수렴: `[hsa04310]`, `[hsa04151] GSK3→CCND1`, `[hsa04630] CycD`, **`[hsa05202] IgH-MAF⇢CCND2` (MM 특이!)** → `[hsa04110] CycD-CDK4/6` **(KEGG)**
- p21 수렴: `p53`, `[hsa04630]`, `[hsa04151]` → `[hsa04110] p21(Cip1)` **(KEGG)**
- ★ **Lenalidomide** 종착점 = `[hsa05202] MYC↓ (via IRF4↓)` + `[hsa04110] p21↑ → G1 정지` **(확립)**

### Hub 7 — p53
- `[hsa04110] p53→p21/GADD45/14-3-3σ` · `[hsa04210] p53→Bax/PUMA/NOXA` · `[hsa04151] AKT→MDM2⊣p53` · `[hsa05202] p53 블록` **(KEGG)**

### Hub 8 — CD38 / Calcium / Apoptosis ★약물 특이 hub★
- `[hsa04020] ADPRcyclase(=CD38)` → cADPR → `RYR` → Ca²⁺ / NAADP → `TPC/TRPML` → Ca²⁺ **(KEGG)**
- `[hsa04020] Ca²⁺` → **Apoptosis 출력** → `[hsa04210]` 연결 **(KEGG)**
- ★ **Isatuximab / Daratumumab** = `[hsa04020] ADPRcyclase(CD38)` 노드에 결합 / ectoenzyme 활성 억제 + (주기전) 면역 매개 살상은 off-map

### Hub 9 — IL-6 / BMM (골수 미세환경 생존 입력)
- `[hsa04630] IL-6family→JAK→STAT3⇢MCL1/Bcl-XL` (BMSC 유래 IL-6) **(KEGG)** → Hub 5
- `[hsa04630] →Ras/MAPK` & `→PI3K-Akt` 분기 **(KEGG)**
- `[hsa04064] RANKL→RANK→NIK→p52` (파골세포 분화 → MM 골질환) **(KEGG)**
- IL-6/PI3K(SHP2) ⊣ **Dexamethasone**-유도 apoptosis → 병용 근거 **(확립)**

### Hub 10 — MM 특이 전사 root
- `[hsa05202] IgH-MAF ⇢ CCND2`(myeloma proliferation) + `ITGB7`(adhesion to bone marrow stroma); `IgH-MMSET→H3⇢CCND2` **(KEGG)**
- 이 블록이 전체 그림의 **유전적 출발점(root)** — CCND2를 통해 `[hsa04110]`로 내려감

---

## 3. 약물 오버레이 — 어디에 ★를 찍나

| 약물 | 분자 표적 | ★ 찍을 위치 | 노드 상태 | 부호/효과 |
|---|---|---|---|---|
| **Isatuximab** | CD38 | `[hsa04020] ADPRcyclase` + `[hsa04210]` 직접 apoptosis 입력 | ✅ 실제 노드 (Calcium) | 결합/억제; 직접 apoptosis(가교 불필요) + ADCC/ADCP/CDC(off-map) |
| **Daratumumab** | CD38 | `[hsa04020] ADPRcyclase` + `[hsa04210]` 직접 apoptosis 입력 | ✅ 실제 노드 (Calcium) | 결합; CDC/ADCC/ADCP 우위, 직접 apoptosis는 FcγR 가교 의존 |
| **Bortezomib** | 26S proteasome β5 | `[map03050] β5` (주표적) → `[hsa04064] IκBα⊗ 차단` + `[hsa04210] ER stress / NOXA↑` + `[hsa04110] p21/p27 안정화` | ✅ 실제 노드 | 억제(⊣); NF-κB OFF + 내인성 apoptosis |
| **Dexamethasone** | GR (NR3C1) | `[hsa04210] Bim↑` + `[hsa04064] NF-κB transrepression` | ⚠️ 표적 미표시, 하류로 진입 | GR은 off-map 주석; 입력 화살표로 Bim·NF-κB에 |
| **Lenalidomide** | CRBN→IKZF1/3 ⊗ | `[hsa05202] MYC↓(IRF4↓)` + `[hsa04110] p21↑/G1` (+ Hub 2 분해기계 투입) | ⚠️ 표적 미표시, 하류로 진입 | CRBN/IKZF는 off-map 주석; 입력 화살표로 MYC·p21에 |

**핵심 정리**
- **실제 KEGG 노드를 가진 약**: 항-CD38 2종(Calcium의 ADPRcyclase 덕분), Bortezomib(Proteasome β5). → 맵 위에 ★를 직접 찍을 수 있음.
- **표적이 맵에 없는 약**: Dexamethasone(GR), Lenalidomide(CRBN→IKZF1/3). → 표적은 맵 바깥(off-map) 박스로 그리고, **하류 효과 노드(Bim / NF-κB / MYC / p21)로 들어가는 입력 화살표**로 표현.

---

## 4. 통합 Figure 레이어 제안 (text mapping → 그림 골격)

5개 약의 작용 계층이 서로 다르므로, 평면 네트워크보다 **4-레이어 + 2개 수렴 hub** 구조가 가장 읽기 좋습니다.

1. **Root 레이어** — `[hsa05202]` MM 블록 (IgH-MAF→CCND2/ITGB7, IgH-MMSET→CCND2): 유전적 출발점.
2. **표면/면역 레이어** — `[hsa04020] CD38(ADPRcyclase)`에 ★Isatuximab/★Daratumumab; 옆에 면역 effector(NK/대식/보체) off-map 박스.
3. **세포질 단백질항상성 레이어** — `[map03050] β5`에 ★Bortezomib; off-map **CRL4-CRBN**에 ★Lenalidomide(→IKZF1/3 ⊗).
4. **전사/핵 레이어** — off-map **GR**에 ★Dexamethasone(→Bim, ⊣NF-κB); Lenalidomide 종착 `[hsa05202] MYC↓` + `[hsa04110] p21↑`.

**Spine(중심 척추)로 그릴 2개 수렴 hub**
- **Hub 5 (Bcl-2 family / 미토콘드리아 apoptosis)** = 모든 약효가 모이는 최종 출력. NOXA(Bortezomib)·Bim(Dex)·CD38 직접 apoptosis·MYC↓(Lena) → Mcl-1/Bcl-XL ⊣ → MOMP → CASP9/3.
- **Hub 1 (NF-κB 생존축)** = Bortezomib·Dexamethasone이 동시에 끄는 축. 골수미세환경(Hub 9 IL-6/RANKL)이 이 축과 Hub 5를 떠받침.

그림에서 부호는 ⊣(억제) / →(활성) / ⊗(분해) / 점선(off-map·면역·synergy)으로 구분하고, synergy 엣지 3개는 점선으로 명시:
- Lenalidomide → CD38↑ ⟹ 항-CD38 효과 강화
- Lenalidomide → NK 활성↑ ⟹ ADCC 강화
- IL-6/PI3K(SHP2) ⊣ Dexamethasone apoptosis (병용 근거)

---

## 5. Caveats (검증/주의)

- **Bortezomib–NF-κB는 단순하지 않음.** *유도형* NF-κB는 확실히 억제하지만, MM에서 *구성적* 정규 NF-κB는 역설적으로 활성될 수 있음(IKKβ 경유). 현재는 ER stress/UPR + NOXA-미토콘드리아 축이 주 세포독성 기전으로 봄. figure 범례에 "NF-κB 억제는 한 기전(논쟁적)"으로 명시 권장.
- **CD38는 정규 세포내 신호 노드가 아님.** Calcium 맵의 ADPRcyclase는 구조적으로 옳은 anchor지만, MM에서 항체의 주기전은 면역 매개 + 직접 apoptosis이고 효소(Ca²⁺ 축) 억제는 부차적. 범례에 표기.
- **Dexamethasone의 apoptosis는 GRE transactivation(Bim 유도) 의존.** NF-κB transrepression은 기여하나 필수는 아님 → NF-κB 화살표를 과대평가하지 말 것.
- **Dex(GR=NR3C1), Lenalidomide(CRBN, IKZF1/3) 표적 노드는 이 10개 맵에 없음** → off-map 박스 + 하류 입력 화살표로 처리(위 3·4절).
- **(추정) 태그가 붙은 연결은 없음** — 본 문서의 엣지는 모두 맵에 그려졌거나(KEGG) 문헌상 확립(확립). 단, off-map 약물 2종의 하류 진입 노드(Bim/MYC/p21)는 효과의 대표 노드를 고른 것이므로, 출판 전 표적-하류 경로를 한 번 더 확인 권장.
- KEGG 맵 버전: 업로드본 기준(예: NF-κB 04064 '24, Wnt 04310 '26, Apoptosis 04210 '25). 노드 좌표를 KEGG에서 직접 받아 KEGG-overlay 버전을 만들 때는 동일 버전 사용 권장.
