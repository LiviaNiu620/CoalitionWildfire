# Coalition Design 论文审阅与下一步计划

审阅对象：*Coalition Design under Heterogeneous AV Objectives: Regimes and Low\-Regret Selection in Mixed\-Autonomy Networks*（2026\-09\-11 版，35 页）

## 总体判断

论文内部逻辑严谨，数值认证（VI gap、omitted\-route slack）做得非常细，这是优点。但以 TS / TAC 的标准看，目前有五个主要风险：

1. **模型本质是势博弈**：当前的 coalition 目标函数使整个 Nash–Wardrop 均衡等价于一个凸优化问题（已符号验证，见 §2.1），而文中写"not assumed to possess a common potential"，与事实不符。
2. **理论偏局部**：核心结果是固定 support、仿射成本下的灵敏度分析，加一个两路径 bottleneck 定理；缺少全局性或网络类层面的结论。
3. **理论和算法脱节**：merge\-direction score M\_Π 是全文理论核心，但 one\-step / two\-step 策略并没有用到它，而是对每个候选合并都重新求解均衡。
4. **实验规模和状态空间偏小**：完整审计只有 n\=4（15 个划分）；Sioux Falls 的主结论只基于 4 个状态；Braess 和 3×3 grid 属于玩具网络；效应量很小（组成对比 ≤0.07% J^UE）。
5. **写作硬伤**：PDF 里所有引用显示为 "?"，References 为空；符号冲突多；免责声明重复出现。

* * *

## 0\. 投稿前必须先修的硬伤

- **参考文献没有编译出来**：正文所有 `\cite` 都显示为 "(?)"，第 35 页 References 为空。检查 `refs.bib` 路径、`\input{response_paper/...}` 的相对路径，以及 bibtex 是否真的跑过。
- **符号冲突**（审稿人一定会指出）：
  - η 同时表示 BPR 指数（η\=4）和 merge path 参数 η∈\[0,1\]；
  - B 同时表示 BPR 常数 0.15 和仿射斜率对角矩阵；
  - N 同时表示公司集合和 bottleneck 总需求；
  - H 同时表示 HDV 下标、bottleneck 中的固定流量标量、以及 H\_C 曲率矩阵；
  - K 同时表示联盟数 |Π|、矩阵 K\_Π 和 K\_H。
- **重复内容**：设计问题 (1)、(41)、(47) 写了三遍；bottleneck 公式 (37) 和 (46) 重复；5.3 节与 6.9 节重复。
- **过度防御性写作**："not endogenous formation"、"not a universal guarantee"、"benchmark evidence, not a theorem" 在摘要、引言、正文、讨论、结论和附录中反复出现十余次。建议集中到一个 Limitations 小节，其他位置删掉。这类句式反复出现也容易让人觉得是 AI 生成的文本。
- **数字精度过高**：18.695、9.9998479×10⁻⁷、1.0308357×10⁻⁶ pp 等，统一保留 2–3 位有效数字。
- **摘要过长、数字过密**（约 330 词）。IEEE 期刊摘要一般限制在 250 词以内；TS 也更偏好讲清楚 insight 的摘要。
- **"480 certified profiles" 的表述**：这其实是 2 网络 × 16 状态 × 15 划分的穷举计数，不是 480 个独立样本。摘要里这样写会显得在夸大证据量。

* * *

## 1\. TODO 1：是否假设所有公司完全知道其他公司和全局信息？

**结论：是的，而且是隐含假设，全文没有明确写出来。** 可以分四个层次说明：

| 层次 | 当前隐含的信息假设 | 文中位置 |
| --- | --- | --- |
| 联盟（玩家）之间 | Nash 均衡要求每个联盟知道或能观察到每条边上其他人的流量 x^{\-C} 和物理成本函数 c\_e。若按一次性博弈解释，还需要博弈结构为共同知识（common knowledge）。 | Def. 1，式 (8) |
| 参考曲率 | b̄\_e \= c'\_e(x\_e^UE)（附录 C），即每家公司的目标函数用到了**全 HDV 用户均衡下的全网斜率**，这等于假设公司知道全局 UE 状态。 | App. C |
| 联盟内部 | Ω\_{C,e} 需要成员之间共享 τ\_i 以及逐边流量 f\_e^i，即联盟内部信息完全共享。 | 式 (12)–(13) |
| HDV | 确定性 Wardrop：完全知道路径时间，没有感知误差。 | Def. 1 |
| **设计者（选 Π 的人）** | **信息要求最强**：必须知道所有 τ\_i、各公司 OD 需求 d^i、HDV 需求、网络与成本函数、HDV active support（用于 A\_H），才能算 x\*(Π) 和 M\_Π。 | 式 (41)、7.6 节 |
| 需求 | 需求批次确定、无随机性。 | 3\.2 节 |

**建议：**

1. **加一个明确的 Assumption（Information structure）**，把均衡解释为**稳态（steady\-state）均衡**，而不是一次性完全信息博弈。由于博弈是单调的，而且实际上是势博弈（见 §2.1），均衡可以通过分散式动态达到：每个联盟只需知道自己的需求、观察到的路段时间和 b̄\_e，不需要知道别家的 τ\_j。可以补一个命题，或一个"分散式投影梯度 / best\-response 动态收敛到同一均衡"的数值演示。这对 TAC 方向很加分。
2. **把玩家的信息和设计者的信息分开写**。玩家层面可以弱化；设计者层面在实践中不可能精确知道 τ，这是审稿人最可能追问的地方。
3. **加一个低成本的鲁棒性实验**：设计者用带噪声或错误分类的 τ̂ 选 Π，再在真实 τ 下计算 regret；或者做一个"最小化最坏情况 regret"的鲁棒划分选择。这个实验工作量小，但能直接回应"信息完全"的质疑。
4. **不完全信息的 Bayesian Nash–Wardrop**（公司之间私有类型）属于独立的大扩展，放进 Future Work 即可，不建议在这篇论文里做。
5. 可选鲁棒性检查：HDV 改用 logit SUE。

* * *

## 2\. TODO 2：模型和实验是否太简单？

**直接回答：实验在"严谨性"上不简单，但在"规模、广度、以及对理论的检验"上偏薄。模型最大的问题不是太简单，而是它的数学结构比文中宣称的更简单。**

### 2\.1 模型层面

**M1. 当前模型是精确势博弈（最重要）。**
由于式 (8) 的物理项用积分 ∫c 表示，而 Ψ\_C 只依赖本联盟的流量，堆叠算子 F\_Π 恰好是下面这个函数的梯度：

> P(z) \= Σ\_e ∫₀^{x\_e} c\_e(s) ds \+ Σ\_C Ψ\_C(y^C)

也就是说，均衡等价于"Beckmann 目标 \+ 各联盟二次项"的凸规划最小化，相当于**带公司特定线性附加费（∇Ψ \= Ω f）的多类 UE**。我用 sympy 检查过：在 BPR、b̄ 冻结的设定下，F 的 Jacobian 严格对称；而**真正的 atomic splittable**（玩家内部化 c'\_e(x)·F\_e^C）的 Jacobian 不对称，非对角差为 (9/5)(f₁−f₂)x²。

这带来三个后果：

- Prop. 3/4（存在性、唯一性）几乎是平凡的，审稿人会认为贡献偏弱；
- 在 BPR 实验中，模型**不是**文中所说的 atomic splittable 博弈，只在仿射成本下两者才一致（附录 A.1）；
- 冻结斜率 b̄\_e \= c'\_e(x^UE) 在 α\=0.9 时可能严重失准，因为此时均衡流量离 UE 很远。**"高渗透率下由 grand 转为 partial"这一主结论，有可能是冻结斜率加 √τ 加权造成的伪影。**

建议：(a) 在正文里明确写出势函数，并加以利用：更快的求解器（Frank–Wolfe、Algorithm B），可以扩展到大网络，还能把设计问题写成下层为凸的双层规划，与收费设计文献对接；(b) 做一组鲁棒性实验，使用状态相关的真实 atomic 内部化，或把 b̄ 更新到均衡点的不动点，检查 regime 转变是否仍然存在。

**M2. "异质目标"实际上是一维的。** τ\_i 只是同一曲率的缩放系数，q\_C\=0，而且只有 0.5 和 1.5 两种类型。标题里的 heterogeneous objectives 有夸大之嫌。建议至少加入：带 q\_C 的能耗/距离/路线偏好；一个追求 SO 的"利他型"运营商（同时内部化对 HDV 的外部性）；以及异质路线权限。

**M3. "Mixed autonomy"的 AV 特征较弱。** 模型中 AV 与 HDV 的唯一区别是 AV 被集中调度，没有容量不对称（比如 Lazar et al. 的 platoon 容量模型）。审稿人可能会说这其实是"车队运营商 \+ 背景交通"的模型。有两种处理方式：加一个 c\_e(x\_H, x\_A) 的扩展；或者在定位上改为 fleet operators（AV、网约车、导航 app）。

**M4. 为什么用"划分"作为调控手段？** 如果设计者能指定 Π 和 γ，为什么不直接收费，或直接调 γ？需要更强的动机，例如反垄断合并审查（正好对应 HHI 这条线）、互操作联盟、数据共享联盟。建议增加 (Π, γ) 的联合设计，并与"最优 AV 收费"做基准比较，说明划分设计能恢复多少差距。

**M5. 成员激励。** 文中完全没有报告各公司自身的成本和 HDV 的成本。建议补充 individual rationality 检查（每家公司在推荐划分下是否比 singleton 好）、core 或 Shapley 分配，以及 AV 与 HDV 之间的分配效应。这部分工作量小，而且能提前回应"企业为什么会接受"这个问题。

**M6.** 治理成本 C\_gov \= γ²×(配对比例) 是人为设定的，Prop. 9 只是下包络，贡献很薄。可以保留，但不宜作为主要卖点。

### 2\.2 理论层面

**T1. 理论偏局部。** Thm 1 基于固定 support 和仿射成本；Prop. 8 本质上是链式法则 dJ/dη；Thm 2 只针对单瓶颈两路径；Prop. 6 的"极小性"比较直接。TS / TAC 审稿人可能会认为这是"灵敏度分析的整理"。**至少需要一个全局性或结构性结果**，可选方向：

- 平行网络或 series\-parallel 网络加仿射成本下，J 在划分细化/粗化上的单调性（在 Huang 2013 "well\-designed network \+ affine ⇒ collusion 有益"的基础上加入 HDV 与异质 τ）；
- 以 K、质量分布、τ 离散度、α 为参数的 PoA 上界；
- 贪心合并的最优性或近似保证，例如在 merge score 符号一致或某种次模性条件下。

**T2. 理论与算法脱节（投入产出比最高的修改）。** 目前的策略没有用到 M\_Π。建议实现 **score\-guided 合并**：在当前均衡处用灵敏度矩阵计算所有候选合并的 M\_Π 并据此排序，只对前几个候选求解均衡。然后报告三件事：(i) M\_Π 的符号对 BPR 下有限合并 ΔJ 的预测准确率；(ii) regret；(iii) 均衡求解次数和计算时间的节省。这样理论就变成了可用的工具。

**T3.** 把 Thm 2 的 Σ\* 和真实网络联系起来。目前"consistent with Theorem 2"只是叙述性的。需要在 Sioux Falls 上构造一个网络层面的有效响应指标，并证明它能预测转变发生的位置。

**T4.** 必须和合谋文献对比定位：Hayrapetyan–Tardos–Wexler (2006)、Huang (2013)、Cominetti–Correa–Stier\-Moses (2009)。"合谋可能有害"已是已知结论，需要讲清楚本文的新意在哪里（HDV 屏蔽、异质目标、成员 OD 约束、设计策略）。

### 2\.3 实验层面

| 问题 | 现状 | 建议 |
| --- | --- | --- |
| 网络规模 | SF（76 链路）\+ 5 边 Braess \+ 3×3 grid | 增加 Anaheim / Chicago Sketch / Winnipeg 中的至少一个。借助势函数结构，用 FW / Algorithm B 求解很容易 |
| 公司数 | 完整审计只有 n\=4 | n\=4 时 one\-step 最多求解约 10 个均衡，穷举也只要 15 个，**看不出计算优势**。需要做 n\=6（203 个划分）、n\=8（4140 个）的小网络穷举，以及 n\=10–20 的启发式比较，并报告求解次数和耗时 |
| 状态空间 | SF 主结论只有 α∈{0.5,0.9}×γ∈{0.5,1}，共 4 个状态 | α∈{0.1,…,0.9}、γ∈{0,…,1}、需求倍数 {0.8,1,1.2,1.5}，画出最优 K 的相图 |
| 效应量 | SF 的 UE–SO 差距只占 TSTT 的 3.8%，组成对比 ≤0.07% J^UE | 同时报告绝对车·小时；增加更拥堵的需求情景 |
| 策略价值 | "always grand" 在 32 个状态中 27 个是精确最优，one\-step 是 29 个 | 需要设计非 grand 最优频繁出现的实例，才能体现选择策略的价值 |
| 异质性 | 两种 τ、等质量、比例 OD | 随机 τ、不等质量、空间专门化组合，与划分审计一起做 |
| 输出指标 | 只有 TSTT 和 recovery | 增加 AV/HDV 分项成本和各公司成本 |
| 鲁棒性 | 无 | 真实 atomic 斜率、SUE、路线集大小 |

### 2\.4 下一步 TODO（按优先级排序）

**Phase A（1–2 周，必须做）**

- [ ] 修复 bib，统一符号，删除重复内容和冗余免责声明，把摘要压到 250 词以内
- [ ] 增加 Information\-structure 假设（稳态解释，玩家信息与设计者信息分开写）
- [ ] 明确写出势函数命题，并修正"not assumed to possess a common potential"
- [ ] 补充合谋 / atomic splittable 文献（Hayrapetyan 2006、Huang 2013、Cominetti 2009、Battifarano & Qian 2023、Toso et al. 2024）

**Phase B（3–5 周，决定能否冲 TS 的关键）**

- [ ] 实现 score\-guided 合并策略，报告符号准确率、regret 和计算节省（T2）
- [ ] 真实 atomic 斜率或更新 b̄ 的鲁棒性实验，验证 regime 转变不是伪影（M1）
- [ ] n\=6/8 穷举与 n≥10 启发式的可扩展性实验
- [ ] SF 的稠密 (α, γ, 需求) 相图，外加一个大网络（Anaheim 或 Chicago Sketch）

**Phase C（2–4 周，加分项）**

- [ ] 更丰富的异质目标（q\_C、SO 型运营商）以及空间专门化与划分审计的组合
- [ ] 个体理性、core、Shapley 检查，以及 AV/HDV 分配效应
- [ ] 设计者 τ 设定错误时的 regret 与鲁棒划分
- [ ] 尝试一个全局定理（平行网络单调性、α 阈值，或贪心保证）

完成 Phase B 之后，再根据理论和计算哪一侧更强来决定投稿期刊。

* * *

## 3\. TODO 3：期刊选择与近年相关论文

### 3\.1 相关论文基准

**Transportation Science**

- **Battifarano & Qian (2023), "The Impact of Optimized Fleets in Transportation Networks", TS 57(4):1047–1068**。与本文最接近：研究一个优化车队加普通用户，提出 critical fleet size（CFS\-UE / CFS\-SO）的概念，实验用 Sioux Falls 加真实 Pittsburgh 网络。结论之一是小渗透率下车队可能让系统变差，这和本文"grand 不总是最优"的逻辑直接相关，**必须引用并对比**。
- **Xu, Chen, Yin & Ye (2021), "Equilibrium Analysis of Urban Traffic Networks with Ride\-Sourcing Services", TS 55(6)**：两个均衡模型加求解算法，并用滴滴实际数据验证。
- **Cummings, Vaze, Ergun & Barnhart, "Multimodal Transportation Pricing Alliance Design", TS**：属于"联盟设计"，使用大规模优化（两阶段分解、坐标下降、warm start）。这是 TS 接受"联盟/划分设计"类问题的先例，但计算规模很大。
- 对 TS 标准的判断：需要新的建模洞见、严谨的理论，**以及真实或中大型网络上的实质性数值研究**。

**IEEE Transactions on Automatic Control**

- **Lazar, Coogan & Pedarsani (2021), "Routing for Traffic Networks with Mixed Autonomy", TAC 66(6)**：用考虑容量不对称的宏观模型给出 PoA 界，核心是一般性定理。
- **Paccagnan, Gentile, Parise, Kamgarpour & Lygeros, "Nash and Wardrop Equilibria in Aggregative Games with Coupling Constraints", TAC**：Nash 与 Wardrop 之间的关系、收敛性和近似界。
- 同类控制方向的工作：Toso, Parise, Frasca & Kibangou (arXiv 2024)，研究协调车队规模对 PoA 的影响，只用 3–7 条链路的例子；Yi & Wei (arXiv 2026)，研究利他 AV 的 VI 存在唯一性，只用 Braess 例子；Mehr & Horowitz (IEEE TCNS 2020)。
- 对 TAC 标准的判断：**要求一般性定理（网络类、界、收敛性），数值部分可以很小**。本文目前的定理偏局部，要投 TAC 必须补 T1 中的全局结果，同时大幅压缩实验部分。

**IEEE Transactions on Intelligent Transportation Systems**

- 覆盖面广、偏应用，大量论文以仿真、RL、数据驱动为主。例如 Eom & Kwon (2025) 在 T\-ITS 发表了 mixed\-autonomy 网络中基于 RL 的"price of autonomous strategy"。
- 对理论深度要求较低，但要压缩到双栏约 10–14 页，并强调算法、规模和政策含义。本文的灵敏度理论在 T\-ITS 读者群中不太能得到充分认可。

### 3\.2 对比与建议

| 期刊 | 主题契合 | 当前稿件差距 | 需要补的工作 |
| --- | --- | --- | --- |
| **Transportation Science** | ★★★ | 网络规模、状态空间、理论与算法的连接、势函数定位 | Phase A \+ B（至少一个大网络，外加 score\-guided 策略） |
| TR Part B（替代选择） | ★★★ | 同 TS，对"模型加理论"更友好 | Phase A \+ B |
| **IEEE TAC** | ★★ | 缺全局定理，实验占篇幅过多 | 重写成理论论文：PoA 或单调性定理、分散式收敛性；实验压到约 2 页 |
| IEEE TCNS / Automatica | ★★ | 同 TAC，要求略低 | 同上 |
| **IEEE T\-ITS** | ★★ | 篇幅需大幅压缩，要突出应用和规模 | Phase A 加大网络实验；理论只保留核心部分 |

**我的建议：** 首选 **TS 或同档次的 TR\-B**，但要在完成 Phase A 和 Phase B 之后再投。目前的版本直接投 TS，大概率会遇到"势函数结构、理论偏局部、实验规模小"这三条意见。只有决定走纯理论路线时才考虑 TAC。T\-ITS 可以作为稳妥的备选。另外，合作者 Ruolin Li 组近期的 Stackelberg–Wardrop / SVO 工作（Wang, He & Li, arXiv 2026）可以在叙事上保持一致并相互引用。

* * *

## Sources

- Battifarano & Qian, TS 2023：https://ideas.repec.org/a/inm/ortrsc/v57y2023i4p1047\-1068.html
- Xu, Chen, Yin & Ye, TS 2021：https://ideas.repec.org/a/inm/ortrsc/v55y2021i6p1260\-1279.html
- Cummings et al., alliance design (TS)：https://arxiv.org/abs/2301.03414 ；https://pubsonline.informs.org/doi/10.1287/trsc.2023.0009
- Xie, Liu & Chen, TS mixed\-autonomy MoD：https://pubsonline.informs.org/doi/10.1287/trsc.2022.1188
- Lazar, Coogan & Pedarsani, TAC 2021：https://coogan.ece.gatech.edu/papers/lazar2018routing.html
- Paccagnan et al., TAC：https://arxiv.org/abs/1702.08789
- Toso, Parise, Frasca & Kibangou (arXiv 2024)：https://arxiv.org/abs/2408.15742
- Yi & Wei (arXiv 2026)：https://arxiv.org/abs/2605.23782
- Mehr & Horowitz (IEEE TCNS)：https://arxiv.org/abs/1901.05168
- Huang, Collusion in Atomic Splittable Routing Games (TOCS 2013)：https://link.springer.com/article/10.1007/s00224\-012\-9421\-4
- Eom & Kwon, T\-ITS 2025：https://api.openalex.org/works/doi:10.1109%2FTITS.2025.3623119
- Wang, He & Li (arXiv 2026)：https://arxiv.org/html/2604.21941
- Coalition formation in congestion games (arXiv)：https://arxiv.org/abs/2410.06797
