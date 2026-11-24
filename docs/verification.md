# 交付前核验

| 层级 | 实际执行 | 结果 |
|---|---|---|
| 有机物构象 | CREST 3.0.2/GFN2，单线程 quick，在两个项目的真实 CLI 入口分别对藜芦醚（ALPB-THF）和特戊酸（ALPB-DMSO）运行 | 两次正常结束，本次分别输出 2、4 个构象，见各自 `examples/*_crest.xyz`。随机采样在其他机器上可产生不同数量 |
| 起始结构 | 用 RDKit 2025.09.6 重建示例几何；检查已注册种子的元素、电荷/多重度奇偶性和作业生成 | 可为两个金属项目生成作业。初猜仍需化学审核与量化优化 |
| 外部文献坐标 | 从 Shimozono 原 SI 的 PDF 导入 m/p-TS1、m/p-TS3 并对其中一个生成文献方法的 ORCA 输入 | 各对分别为 96、89 个原子，组成与注册记录一致。原 PDF 由用户本地提供 |
| 化学动力学 | SciPy 1.17.0 / NumPy 2.3.5 积分两张网络，libSBML 5.20.5 检查模型，在 COPASI 4.48.311 Python 引擎逐一导入运行同一批 SBML | 共 24 个假说/条件组合；比较每个物种的终点浓度。容差为绝对 2×10⁻⁸ M、相对 2×10⁻⁵；各项通过 |
| 数据和实验选择 | 对示例数据做两个项目的拟合/留出、质量收支、配对判定；测试无钴/无铁、暗对照、未测值、错误自旋及不配平反应 | 七组回归测试通过。拟合测试使用标为 synthetic 的数据 |
| ORCA | 验证输入生成和结果读取逻辑，使用明确标记的模拟日志测试收敛、虚频和热力学检查 | 此环境没有授权 ORCA 可执行文件，未运行新的 DFT 作业、计算研究势垒或验证目标收率 |

实际测试使用 Linux / Python 3.12。Windows/WSL 的安装指令供使用者选择，没有在一台全新的 Windows 主机上重装验证。COPASI 核验使用数值引擎。图形界面可按 README 手工导入。示范的反应范围和基准记录见两个项目的说明。

## 在本地复验

```bash
python -m unittest discover -s tests -v
python -m chemflow demo --out runs/validation
python -m pip install python-copasi==4.48.311
python tools/check_copasi.py runs/validation
```

默认安装只需 requirements.txt 中的三个库。最后两条用于额外的独立数值复核。查看真实证据时请从 `benchmarks.csv` 的原始出处开始。拟合合成数字、两个求解器结果一致和物料守恒分别检验程序行为，均不能说明哪条机理在真实反应中占主导。

