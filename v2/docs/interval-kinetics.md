# 区间动力学与端点回溯

工作任务使用 `engine=analysis`、`stage=mechanism`、`operation=interval_kinetics`。输入为一个 JSON 文件，输出为一个声明的 JSON 报告。计算节点需要 NumPy 和 SciPy。原始数据、审核过的机理模型和观测区间随研究包提供。机制输入可用现有 `$artifact` 绑定上游能量或原件记录。每个竞争机理单独提交，在同一组训练端点和条件下比较可行集。

## 模型和条件

| 字段 | 含义 |
|---|---|
| model_source_sha256 | 经审核的冻结模型的 SHA-256 |
| temperature_K, end_h | 本组各实验臂共用的温度和观测时长 |
| species_gibbs_kcal_mol | 物种标准溶液吉布斯能，共用同一参考和 1 M 标准态 |
| parameter_grid | 预先指定的能量参数及其递增取值。能量字段中的字符串引用同名参数。同一参数可出现在多个反应或物种中 |
| reactions | 每步给出 id、lhs、rhs 两侧的正整数计量数。热反应使用 type=thermal 和 ts_gibbs_kcal_mol，填写与物种使用同一参考的过渡态自由能。光激发使用 type=photon 和 quantum_yield，反应物为吸光体 |
| conserved | 守恒池到物种系数的映射，覆盖全部需核对的物料与金属，包括所有金属配合物 |
| henry_bar_L_mol | 可选，气体物种在指定温度下的 Henry 挥发度，单位 bar L mol⁻¹ |
| light.epsilon_L_mol_cm | 所有吸光体在入射波长的消光系数，须包含竞争吸收者 |
| solver_yield_tolerance | 普通与收紧容差积分得到的终点收率允许的绝对差，收率以 0 至 1 表示 |
| energy_coordinates | 可选，将参数线性组合为相图能量坐标。系数应依据冻结的能量定义推导。例如取代过渡态能量减去游离反应物能量，可表示取代势垒加结合自由能 |
| phase_contrast | 可选，primary、secondary 指向两套实验条件，minimum_difference 指定冻结的收率差阈值 |

每个 conditions 条目含有 role、initial_total_M 和 observable。初始值为每升液相所含的**总摩尔数**，包括平衡气相。observable 提供产物名 product 及定义最大理论收率的守恒池名 limiting_pool。无铁臂用 absent_pools 显式声明金属守恒池，其全部形态初始总量必须为零。气液平衡需要 liquid_L、headspace_L。光化学需要 incident_photons_M_h（按液相体积计的入射光子量）和 path_cm。

训练端点 role=training，要求 yield_interval（两个分数）和 observation_sha256。区间包括实验误差或可追溯的测量界限。无检出限的定性观察仅能在有依据的单独评价规则下裁决，不能任意换算为精确零收率。role=holdout 与 role=design 只携带实验条件，输入中出现观测区间会被拒绝。留出测量通过独立的 freeze_prediction 和 score_holdout 交接。

## 方法与输出解释

热反应的正逆速率由同一过渡态和物种能量计算。质量作用使用各物种溶液活度近似值和 1 M 标准态，浓度单位为 mol/L，时间单位为 h。每次求导以 Henry 定律和理想气体方程分配气液两相，守恒核对使用总量。光驱动步根据 Beer–Lambert 吸收和各吸光体贡献分配已吸收光子。该实现支持单波长、均匀液相、快速气液分配和固定温度，不包含空间光场、气液传质迟滞或激发态消光系数的自动预测。

引擎枚举全部网格组合，验证反应计量守恒，并对每种实际积分条件核查浓度、守恒量、物理收率和收紧求解器后的数值稳定性。任一数值错误都会中止分析。只有同时满足全部训练端点的参数组合才产生留出预测。输出保留参数组合、失败端点、可行参数支持、观测包络、每个参数在其他参数固定时的训练端点差分敏感度，以及预测分歧最大的后续实验。每个可行点另报告归一化参数轴上的训练端点差分雅可比矩阵奇异值、受积分容差限制的局部数值秩和未分辨维数。该指标只描述网格尺度上训练端点的局部敏感方向。笛卡尔网格的范围、密度及测量区间决定这些包络；它们不是统计置信区间或对连续空间的穷尽证明。空可行集以 indeterminate 阻止依赖推进。

最小可复核实例见 tests/test_interval_kinetics.py。测试包含有解析解的可逆网络、带顶空气相的平衡终点、留出隔离和有光竞争的无金属对照。真实研究使用审核过的原始数据、能量和照射条件。
