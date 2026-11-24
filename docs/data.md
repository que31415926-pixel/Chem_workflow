# 接入测量并决定下一批合作实验

## 先保留原始单位，再输入模型

1. `handoff` 输出下一批的 `next_batch.csv` 与空白 `wet_records_blank.csv`。记录瓶号、独立实验日、实际前体浓度/总钴、溶剂与体积、光路或顶空、标品、回收、定量限和原始数据路径。同一瓶多次进样保留原始记录，统计时合并为一个反应样。
2. 物料衡算填 `projects/{co,fe}/examples/material_blank.csv`。所有字段以 **mmol/整瓶** 填入；例如 `co_atoms_in` 以钴原子计，双核催化剂每摩尔计两摩尔钴。校准所得的气相浓度先折算成整瓶物质的量，溶解气单列。`other_b_atoms` 直接填未被其他列覆盖的硼原子 mmol，`aryl_di` 填双硼化产物分子 mmol。铁表的 `*_side_eq` 按该片段折算，不重复计入二聚物、盐或酸。只有经测量或由设计确定确实为零时才填 0；未测项留空，依赖该项的结果也留空。`blank_corrected_H2_per_target` 仅比较空白校正后的 H₂ 与目标产物，副反应和其他异构体仍会影响解释。
3. 动力学拟合表内 `value_M`、`sigma_M` 是 **mol/L 液体积**，时间为小时。输入真实反应前先把模型中对应物种的初始值、有效光通量、条件扰动倍率和检测误差改成实际值，保留改动记录。Fe 的 `Salt` 用真实有效浓度。COPASI 归一化模型容积为 1 L，即使反应实物是微量瓶，预测浓度一致。`Hgas` 同样以气体摩尔数除以反应液体积表示。
4. 可在本地私有副本增加实际条件、物种与候选副反应。`network.json` 中 `lhs`、`rhs` 为化学计量数，`k` 为质量作用系数。一级单位 h⁻¹，二级 M⁻¹ h⁻¹。自二聚写为 `v=k[R]²`，自由基消耗是 `2v`。初始值/模型/测量表必须对应同一反应及温度。引入新池时更新守恒映射，`simulate` 会检查每一步。导出的每组 SBML 与 SciPy 使用同一化学计量和条件。

## 试跑算例

下列输入名带 `synthetic` 的表都是程序验证用的合成数字。

```bash
python -m chemflow balance co --data projects/co/examples/synthetic_material.csv --out runs/co-balance
python -m chemflow balance fe --data projects/fe/examples/synthetic_material.csv --out runs/fe-balance
python -m chemflow fit fe --data projects/fe/examples/synthetic_observations.csv --out runs/fe-fit
python -m chemflow simulate fe --calibrated runs/fe-fit/calibrated.json --out runs/fe-after-fit
python -m chemflow paired-fe projects/fe/examples/synthetic_paired.csv
```

钴的拟合把 `fe` 换为 `co` 并使用对应目录中的数据。拟合输入列含 `batch, condition, time_h, species, value_M, sigma_M, kind, split, origin`。时间程可来自不同独立日；`kind=value` 为定量值，`kind=upper` 为经校准的浓度上界（含义和口径在原始记录注明）。超过上界才产生单侧损失，这是一种筛查评分，并非完整删失似然。至少留出一个**整个条件**作 `holdout`；不要把同一曲线的晚时间点移作外部验证。汇总不同文献平台、反应物、产率口径的记录很容易得到误导性排序。

`fit` 为每个备选解释仅拟合预设的一个速率，其余参数固定；默认搜索范围是初值的 0.01 至 100 倍。在 `scores.csv` 同时检查 `holdout` 的标准化均方根误差、是否碰到边界，在 `residuals.csv` 看观测种类、时间程和残差方向。样例真值由其中一个假说生成，恢复它只验证程序。真实数据下几个模型相近、均碰边界或系统性偏离时，应修改可检验的假设或引入新观测。给出置信区间、光强/气液误差与批次效应时，应另做参数范围扫描、重采样或正式的联合拟合。

`diagnostics.csv` 的每行对应一对假说。先计算干预减基线的变化，再取两种解释的差值，除以预设的**配对差标准差**。优先考虑重复性好、确实可测、能区分仍有可能的两种解释的项目。不能辨识的有效物种（如未指认钴峰）不应仅凭分数安排实验。评分没有对每种实验的操作成本、失效概率或后验概率建模。`provenance.json` 标明是否为默认演示或部分拟合。

## 把真实证据转换成实验清单

复制 `projects/co/state.json` 或 `projects/fe/state.json` 到本地 `private/`，保留原始模板。项目阶段状态由研究者根据 QC、独立日和配对结果填写。初始状态全为 `unknown`。可用值和阶段定义在 `wet_plan.json`。

| 项目 | 状态字段 | 示例含义 |
|---|---|---|
| 钴 | controls · parent · balance · signal · heldout | 阳性与定量通过；母体已满足目标/有缺口/待定；收支；偏向氢限制、C4 或 C–H；留出结果 |
| 铁 | controls · salt · target · supply · heldout | 阳性通过；可投原盐/只能改环/均不足；目标可行/低效/待定；绕过前体后是否恢复；留出结果 |

```bash
python -m chemflow handoff fe --state private/fe_state.json --out runs/fe-next
python -m chemflow handoff co --state private/co_state.json --out runs/co-next
```

例如，铁项目设 `controls=pass`、`salt=modified_only`，其余保持 unknown，得到改环同骨架的最小比较，并注明原环体系不可评价。钴项目设 `controls=pass`、`parent=shortfall`、`balance=pass`，在尚未定位原因时会得到一项配对干预，按模型及测量能力选择。单次演示输出或拟合优劣都不会自动推进实验阶段。

`paired-fe` 接收同日一级参照与二级目标、相同盐骨架/酸和反应条件下的**百分产率区间**，每个独立区组一行。定量值的上下界相等，未定量可写 [0,经标定的上界]。函数返回目标和“一级减二级”的均值范围、保守的双侧 95% t 区间包络。此算法基于独立配对区组及均值 t 近似，删失时以最大可能样本方差构造保守包络；小样本不稳。`works` 暂按目标下限 ≥30%；`low` 暂按目标上限 <30% 且配对差下限 >15 个百分点；其余为 `uncertain`。这两条线和至少三个独立区组是本项目的预设操作线，可在看结果前修改。只接受 QC 已通过的输入，失败批应单独归档并明确排除理由。钴的目标线和对照要求见钴项目说明。

