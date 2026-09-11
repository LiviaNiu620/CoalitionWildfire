# 野火扰动交通网络中的商业 AV Coalition Design：研究计划与执行路线

## 先看结论

本项目的核心仍然是 coalition design，不是 wildfire evacuation simulation，也不是 AV rescue dispatch。野火是外生的网络扰动和公共信息状态；商业 AV 继续服务已接受的 pickup-dropoff 订单；订单 OD 和公司归属保持不变；普通 HDV 仍然是能够获得公开道路信息、可以自行改道、但不受直接控制的交通参与者；监管者可以强制公司参加临时 operational coalition，但不逐车指定路线，也不重新分配订单。

当前原 coalition continuation model 已经可以作为可靠 baseline，但还不是完整 wildfire application model。现有代码支持静态网络、普通 HDV、固定 OD、identity-preserving coalition、symmetric PSD governance、BPR continuation 和 certified partition scans；仍缺少 hazard-dependent edge state、容量退化、风险属性、critical-corridor metrics、hazard-specific reference slope 和信息误差实验。

因此当前状态是：

```text
coalition baseline: ready
potential correction: complete
ordinary-HDVs wildfire scope: locked
wildfire smoke test: complete
Sioux Falls wildfire adapter: pending
formal wildfire evidence: not started
```

## 一、研究故事

野火发生后，道路网络可能出现关闭、容量下降、风险上升和公开道路信息快速变化。商业 AV 订单并不会在数学上自动消失；在一个短时间窗内，已经接受的订单仍然需要由原公司服务。与此同时，普通 HDV 根据公开道路信息自行改道，但监管者不能逐车控制它们。

AV 的特殊性不在于承担救援任务，而在于其交通流可以在公司层面集中调度。问题在于，这种可编程性分散在多个运营商之间。不同公司拥有不同的订单 OD、路线暴露和管理目标，并且各公司通常不会自动 internalize 其他公司的拥堵影响。

因此监管者面对的是一个组织设计问题：在不能接管订单、不能直接指定车辆路线的条件下，是否应当强制公司参与临时运营联盟？如果需要协调，应该是 singleton、partial coalition 还是 grand coalition？

核心政策问题是：

> When wildfire-induced network disruption changes road accessibility and public information, which regulator-mandated temporary operational coalition structure among commercial AV operators minimizes network congestion, emergency-critical corridor overload, order delay, risk-weighted traffic, and governance burden?

## 二、研究边界

### 本文研究

- 商业 AV operators 的 coalition structure；
- 已接受商业订单的固定 pickup-dropoff OD；
- ordinary HDV 的公开信息下 route response；
- 野火导致的道路关闭、容量退化和风险状态；
- regulator-mandated temporary operational coalition；
- coalition governance intensity；
- coalition response、HDV screening 和 activated response；
- TSTT、critical-corridor overload、订单服务损失和治理成本。

### 本文不研究

- 居民疏散需求本身；
- evacuee population 或 evacuation OD；
- AV rescue fleet；
- vulnerable-household pickup；
- shelter assignment；
- 911 dispatch；
- fire-spread PDE；
- 公司法律合并；
- 自愿 merge-and-split stability；
- 新订单接受、价格和乘客平台选择；
- inter-epoch empty-vehicle repositioning。

## 三、模型状态

### 3.1 野火状态

对 hazard scenario \(\omega\) 和时间窗 \(t\)，定义：

\[
s_t^\omega=(E_t^\omega,t_{e,t}^{0,\omega},\kappa_{e,t}^\omega,r_{e,t}^\omega,\mathcal I_t^\omega,E_{\mathrm{critical}}^\omega).
\]

其中：

- \(E_t^\omega\)：仍然开放的道路集合；
- \(t_{e,t}^{0,\omega}\)：野火状态下的自由流时间；
- \(\kappa_{e,t}^\omega\)：动态道路容量；
- \(r_{e,t}^\omega\)：风险或未来关闭风险；
- \(\mathcal I_t^\omega\)：公开道路和关闭信息；
- \(E_{\mathrm{critical}}^\omega\)：应急敏感道路集合。

第一版使用 static hazard snapshots。rolling horizon 是后续扩展，不进入第一版理论假设。

### 3.2 流量

\[
x_{e,t}^\omega
=u_{e,t}^\omega+h_{e,t}+\sum_{i=1}^{n}f_{e,t}^i.
\]

其中：

- \(u^\omega\)：外生背景流；
- \(h\)：普通 HDV 流；
- \(f^i\)：公司 \(i\) 的商业 AV 订单流。

HDV 不改名为 evacuee，也不额外建立居民疏散 population。

### 3.3 商业订单约束

公司 \(i\) 在当前短时间窗内接受的订单为 \(d_{w,t}^i\)，满足：

\[
\sum_{r\in\mathcal R_{w,t}^\omega}y_{r,t}^i=d_{w,t}^i.
\]

订单归属和 pickup-dropoff OD 不变。coalition 可以联合规划这些订单的路径，但不能把订单改成其他 OD，也不自动改变订单所有权。

### 3.4 Coalition feasibility

对 coalition \(C\)：

\[
\mathcal Y_{C,t}^{\mathrm{id}}
=\{y_t^C\geq0:\Gamma_{C,t}y_t^C=d_t^C\}.
\]

主模型使用 identity-preserving coordination。Operational pooling 只作为后续 extension，因为 pooling 会改变服务 authority。

### 3.5 Coalition objective and potential

coalition 目标保留原结构：

\[
\Phi_{C,t}(y^C;z_{-C})
=\sum_e\int_0^{F^C_{e,t}}c_{e,t}^\omega(x^{-C}_{e,t}+s)\,ds
+\Psi_{C,t}(y_t^C).
\]

固定参考曲率下：

\[
\Psi_{C,t}(y^C)
=\frac12(y^C)^\top H_C(\gamma)y^C+q_C^\top y^C.
\]

对称 PSD governance family 对应凸势函数：

\[
\mathcal P_{\Pi,t}^\omega
=\sum_e\int_0^{x_{e,t}^\omega}c_{e,t}^\omega(u)\,du
+\sum_{C\in\Pi}\Psi_{C,t}(y_t^C).
\]

因此需要区分：

- 一般 behavioral model：不预设 common potential；
- 当前 fixed-reference numerical model：是 convex potential-based continuation model；
- nonlinear BPR true own-travel-time atomic game：与 fixed-reference model 不自动相同。

### 3.6 HDV 信息结构

HDV：

- 可以看到公开道路旅行时间；
- 可以看到公开道路关闭状态；
- 可以获得公开风险信息；
- 可以自行改道；
- 不受监管者逐车控制；
- 不知道公司私有订单和私有 objective type。

基线使用 deterministic Wardrop。信息延迟、bounded rationality 和 logit response 作为 robustness。

### 3.7 监管者目标

监管者选择 temporary operational partition \(\Pi\) 和 coordination intensity \(\gamma\)，评价：

\[
L^\omega(\Pi)
=J_{\mathrm{TSTT}}^\omega(\Pi)
+\lambda_AJ_{\mathrm{order}}^\omega(\Pi)
+\lambda_RJ_{\mathrm{risk}}^\omega(\Pi)
+\lambda_CJ_{\mathrm{critical}}^\omega(\Pi)
+\lambda_GC_{\mathrm{gov}}(\Pi,\gamma).
\]

这里不使用 evacuee-specific objective。\(J_{\mathrm{critical}}\) 表示应急敏感道路的拥堵或容量占用，不等于居民疏散流。

## 四、研究假设

### H1：Coalition effect

在 hazard state、订单 OD 和 HDV demand 相同的情况下，改变 coalition structure 会改变商业 AV route response 和网络拥堵。

### H2：State dependence

最优 coalition count 取决于道路关闭、容量退化、订单 corridor exposure、HDV response 和治理成本。

### H3：Partial coordination

Grand coalition 不一定普遍最优；在订单空间异质、路线暴露不同或治理成本较高时，partial coalition 可能更好。

### H4：HHI insufficiency

相同公司质量和 HHI 下，不同订单 OD exposure、objective composition 和 coalition grouping 仍可产生不同结果。

### H5：Reference-slope sensitivity

高 AV load 下的 coalition transition 可能依赖 \(\bar b\) calibration，必须通过 state-dependent slope 和 reference-state robustness 检验。

### H6：Information robustness

公开道路信息延迟和 regulator-side state error 会增加 partition-selection regret，但不一定消除 coalition design 的价值。

### H7：Score-guided selection

固定-support merge score 可以作为 candidate screening 工具，但不能未经 finite-endpoint audit 就声称 global optimality。

## 五、实验阶段

### Phase 0：Baseline integrity（已完成）

目的：冻结原 coalition model。

执行：

- 运行所有现有 validators；
- 验证 potential gradient；
- 验证 Jacobian symmetry；
- 保存原始 JSON；
- 固定 Python、numpy、scipy、networkx、sympy 版本；
- 记录 VI 和 omitted-route tolerances。

数据：现有 Sioux Falls、Braess、directed grid 及其 committed results。

验收：

- 所有 validators 通过；
- potential residual 通过；
- no baseline JSON 被覆盖。

状态：已完成。

### Phase 1：Hazard adapter smoke test（已完成最小版本）

目的：确认野火只改变网络状态，不改变 HDV 身份和 coalition 定义。

执行：

- 三条并行道路；
- ordinary HDV；
- 四家公司商业订单；
- singleton versus grand coalition；
- no hazard；
- capacity degradation；
- single-edge closure。

验收：固定路线 VI gap 通过，订单 OD 不变，关闭道路不再被使用。

状态：最小 smoke test 已完成；尚未移植到 Sioux Falls。

### Phase 2：Sioux Falls static hazard adapter

目的：把 hazard state 接入真实 benchmark network。

执行：

1. 建立 `hazard_scenarios.json` schema；
2. 加入 open/closed edge；
3. 加入 capacity multiplier；
4. 加入 free-flow multiplier；
5. 加入 risk score；
6. 加入 public-information metadata；
7. 加入 critical-edge metadata；
8. 确认 shortest-path oracle 排除 closed edges。

首批 scenario：

- no hazard；
- mild capacity degradation；
- local closure；
- critical corridor closure；
- multiple corridor degradation；
- asymmetric regional degradation。

验收：

- no-hazard 输出与 baseline 一致；
- 所有 route 不含 closed edge；
- cost finite；
- VI gap <= `1e-6`；
- relative omitted-route slack <= `1e-6`。

### Phase 3：Pilot full partition audit

目的：确认 hazard state 是否真的改变 coalition ranking。

设计：

- n=4；
- 全部 15 个 partitions；
- ordinary HDV；
- fixed commercial AV OD；
- no hazard、capacity loss、critical closure、multiple degradation；
- 两个 AV load；
- 两个 gamma。

报告：

- best partition set；
- best K；
- TSTT；
- order delay；
- critical-edge overload；
- grand/singleton regret。

验收：如果没有 ranking change，也要如实报告 null result，不为了支持 partial regime 而调整 scenario。

### Phase 4：Main hazard regime scan

目的：构建 wildfire-state coalition regime map。

扫描因素：

- hazard severity；
- closure topology；
- capacity multiplier；
- AV load；
- order exposure to critical corridors；
- objective composition；
- gamma；
- governance burden。

主结果：

- hazard state x AV load 的 best-K heatmap；
- order exposure x coalition composition 的 recovery/TSTT comparison；
- critical-edge overload comparison；
- governance-cost lower envelope；
- fixed grand versus state-contingent coalition regret。

### Phase 5：Reference-curvature robustness

目的：判断高 AV load 下 partial regime 是否由固定参考斜率造成。

比较：

\[
\bar b_e^\omega
=\kappa_b c_e^{\omega\prime}(x_e^{\mathrm{ref},\omega}),
\qquad
\kappa_b\in\{0.25,0.5,1,2,4\}.
\]

参考状态：

- hazard-state HDV-only equilibrium；
- hazard-state no-coordination equilibrium；
- nominal baseline state。

同时使用代码中的 optional `management_slope_mode="state"` 做 nonlinear diagnostic。

报告：

- best K；
- transition location；
- ranking agreement；
- cross-model regret；
- slope mismatch \(c'(x^*)/\bar b\)。

### Phase 6：Information robustness

目的：验证 informed-but-uncontrolled HDV 和 regulator information error 的影响。

HDV 情形：

- zero delay；
- one-window delay；
- multiple-window delay；
- partial closure information；
- bounded-rational or logit sensitivity。

监管者误差：

- \(\widehat\tau_i\)；
- \(\widehat d_w^i\)；
- \(\widehat\kappa_e\)；
- \(\widehat r_e\)；
- \(\widehat G^\omega\)。

报告 selection regret、false-grand rate、false-partial rate 和 critical-edge overload。

### Phase 7：Merge-score policy

目的：把理论的 \(\mathcal M_\Pi\) 接入算法。

比较：

- full partition audit；
- TSTT-based greedy merge；
- score-guided greedy merge；
- score-guided two-step；
- noisy-state score policy。

报告：

- score sign accuracy；
- rank correlation；
- finite regret；
- equilibrium solves saved；
- support-change fallback rate；
- runtime。

### Phase 8：Scalability and LA grounding

目的：验证方法的计算价值和现实输入适配性。

规模：

- n=4：全 audit；
- n=6：203 partitions；
- n=8：selected-state audit；
- n>=10：beam search、score-guided search 或 local search。

网络：

- Sioux Falls：主机制与 hazard audit；
- Anaheim：中等规模城市网络迁移；
- Chicago Sketch：只做有限扩展性演示。

现实校准：

- fire perimeter 和 road closure；
- traffic speed/flow；
- critical corridor；
- AV service-area 和订单 OD 的可获得部分。

## 六、数据方案

### Level 1：仿真机制数据

- Sioux Falls；
- Braess；
- directed grid；
- synthetic hazard snapshots；
- synthetic commercial AV OD portfolios；
- controlled \(\tau_i\)；
- controlled \(\gamma\)。

用途：验证机制、理论、算法和稳健性。

### Level 2：公开网络和野火数据

- OpenStreetMap：路网拓扑和道路等级；
- Caltrans PeMS：速度和流量；
- LADOT data：城市道路状态；
- CAL FIRE：火灾 perimeter 和 incident；
- NIFC：历史火灾边界；
- NASA FIRMS：hotspot；
- NOAA：风和天气；
- Caltrans QuickMap：道路关闭；
- LA County emergency GIS：应急道路和关键设施。

用途：把真实 wildfire event 转换成 edge availability、capacity 和 risk scenario。

### Level 3：AV 和制度数据

理想数据：

- company-level accepted order OD；
- order timestamp；
- route and delay；
- cancellation；
- fleet size；
- service area；
- emergency operation policy。

这些数据大概率需要运营商或公共部门合作。如果无法取得，应把 AV OD、\(\tau_i\) 和 governance cost 标为 calibrated/synthetic，而不是 empirical estimates。

## 七、论文结果如何分析

每组结果至少报告三种量纲：

1. 绝对 TSTT 或 vehicle-time；
2. 相对于 \(J^{UE}\) 的百分比；
3. recovery points 或 partition regret。

重点分析：

- coalition 是否改变网络性能；
- effect 是否依赖 hazard topology；
- order corridor exposure 是否比 HHI 更有解释力；
- grand coordination 是否存在尾部失败；
- partial coalition 是否在治理成本下更有优势；
- state-dependent slope 是否改变 transition；
- information error 是否改变监管者推荐；
- merge score 是否减少实际 equilibrium solves。

任何结果如果只在某个 reference slope、某个 hazard scenario 或某个 objective design 下成立，都必须写成条件性结论。

## 八、论文 Introduction 结构

1. 野火造成道路网络扰动和应急敏感 corridor；
2. 商业 AV 订单仍然存在并继续占用道路容量；
3. 普通 HDV 能看到公开道路信息但不能被直接控制；
4. AV 路由具有公司级可编程性，但控制权分散；
5. 监管者可以强制临时 operational coalition，但不能直接指定路线；
6. grand、partial、singleton 的权衡来自订单 OD、目标异质性、HDV response 和治理成本；
7. 引出 coalition response、HDV screening 和 activated response；
8. 提出研究问题和四项贡献；
9. 明确不研究 evacuee demand、rescue 和 voluntary formation。

## 九、建议贡献

1. **Hazard-conditioned coalition model**：在野火扰动网络中保留商业 AV 订单 OD 和 identity-preserving coalition。
2. **Potential and response characterization**：证明固定参考对称 PSD governance 下 continuation equilibrium 是凸势最小点，并刻画 HDV screening 后的 activated response。
3. **State-dependent coalition regimes**：展示道路扰动、订单空间暴露、目标异质性和治理成本如何改变最优 coalition structure。
4. **Robust selection and computation**：比较 reference/state slope、信息误差和 merge-score policy 的 regret、准确率和求解成本。

## 十、当前执行状态与下一步

已完成：

- baseline validators；
- canonical regime check；
- potential structure certificate；
- ordinary HDV wildfire smoke test；
- proposal、model specification、experiment protocol 和 traceability plan。

下一步：

1. 实现 Sioux Falls static hazard adapter；
2. 运行 no-hazard regression；
3. 加入 capacity degradation 和 single critical-edge closure；
4. 对 n=4 全部 15 个 partitions 做 pilot audit；
5. 再决定是否启动完整 hazard regime scan。

在 Sioux Falls hazard adapter 和 pilot audit 通过之前，不应在论文中写入任何 wildfire-specific coalition ranking 或 grand-to-partial transition 结论。
