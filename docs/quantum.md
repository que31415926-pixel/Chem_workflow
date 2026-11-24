# 从结构到可复核的量子化学结果

默认工作目录是仓库根目录。在 PowerShell 中激活相同环境或用 README 中的虚拟环境 Python 替换下列 `python`。每个输出路径使用一次，失败的外部任务会留下 `.partial-*` 日志。先运行 `python -m chemflow list co` 或 `list fe` 查看物种名、净电荷和允许的自旋多重度。

## 1　构象及坐标

```bash
python -m chemflow conformers co veratrole --quick --cores 2 --out runs/co-crest
python -m chemflow conformers fe pivalic_acid --quick --cores 2 --out runs/fe-crest
python -m chemflow prepare co veratrole --xyz runs/co-crest/crest_conformers.xyz --frame 1 --out runs/co-veratrole-c1
python -m chemflow prepare fe pivalic_acid --xyz projects/fe/examples/pivalic_acid_crest.xyz --frame 2 --out runs/fe-acid-c2
```

`conformers` 调用真实 CREST/GFN2，记录命令和日志。没有安装 CREST 时也可以用随包实跑的两个构象集合练习选帧。编号从 1 开始。金属自旋与配位几何需要独立构建，脚本对金属跳过 CREST 预设。铁项目的 CREST 用 DMSO 作极性筛选代理；ORCA 的默认介质以 ε=30、n=1.40 示范近似，两者都不代表真实 DMA/iPrOH 混合物。运行 `--solvation low` 或 `high` 比较 ε=22/38，研究时根据实测溶剂性质和显式溶剂/对离子模型复核。钴项目默认 THF。

结构由 RDKit ETKDGv3/UFF 和可复查的几何放置生成。母体钴的中心 N、两个卡宾 C 和外侧 N-iPr 定义于 `tools/build_seeds.py`。Fe–OEP 基于游离八乙基卟啉连接表，移去两个内侧 N–H 后加入铁，分别生成 FeII、FeIII–Cl 和铁烷基。Fe–Cl 是预催化剂；不应把 Cl 保留在所有步骤上。起始文件是坐标猜测。读图检查金属–供体距离、配体翻转、取代、轴配位和电荷，随后用量化优化。需要独立重建时在可写副本运行 `python tools/build_seeds.py`，会重写 `xyz/` 下同名种子。连接表来源见 [PubChem CID 102311](https://pubchem.ncbi.nlm.nih.gov/compound/102311)。

## 2　优化、频率和自旋

```bash
python -m chemflow prepare co co_me --mult 1 --out runs/co-me-s1
python -m chemflow prepare co co_me --mult 3 --out runs/co-me-s3
python -m chemflow prepare fe fe_sec_bn --mult 2 --out runs/fe-bn-d
python -m chemflow prepare fe fe_sec_bn --mult 4 --out runs/fe-bn-q
python -m chemflow prepare fe fe_sec_bn --mult 6 --out runs/fe-bn-s
python -m chemflow run-orca runs/fe-bn-d --exe /opt/orca/orca
python -m chemflow audit runs/fe-bn-d
```

请把演示路径替换为实际安装的 **ORCA 6.1 系列可执行文件绝对路径**。默认单核、每核 1500 MB，金属生产任务可据集群资源调整 `--cores` 与 `--memory`（MB/核，给操作系统与额外数组留余量）。`--method screen` 使用 r²SCAN-3c；`--method check` 用 PBE0-D4/def2-TZVP。后者既可执行独立优化/频率，也可用 `--stage sp --xyz .../job.xyz` 对**已经优化**的几何算单点。只做不同电子单点时，两套电子能均不应被直接当作活化自由能。如采用单点和低层热修正组合，应明确记录来源并单列。

`audit` 检查正常结束、优化收敛、虚频数、热力学量与输入文件哈希。配方/构型应先在输入端冻结。有 TS 输出时亲自查看唯一虚频的振动方向；有弱虚频时重新优化或加密积分。比较低自旋时还应检查波函数稳定性、自旋污染、金属/配体自旋密度，必要时用破缺对称初猜或多参考方法复核。默认命令不会自动搜索不同电子态。若不同模型改变能垒排序，保留两种解释。

## 3　过渡态及自由能

`co_bpin_ar`、`sh2_pri`、`sh2_sec` 给出完整元素和原子顺序，两个片段**初始彼此分离**。使用分子编辑器、约束扫描或 NEB 生成所需构型。例如先建立二级苄基碳/叔丁基碳接近、Fe–苄基键适当拉长的初猜并存为 `private/fe_ts_guess.xyz`。总自旋和电子耦合需要另行检验；不能把任一面上未找到鞍点解释成反应不存在。

```bash
python -m chemflow prepare fe sh2_sec --stage ts --mult 3 --xyz private/fe_ts_guess.xyz --out runs/fe-ts
python -m chemflow run-orca runs/fe-ts --exe /opt/orca/orca
python -m chemflow prepare fe sh2_sec --stage irc --mult 3 --xyz runs/fe-ts/job.xyz --out runs/fe-irc
python -m chemflow run-orca runs/fe-irc --exe /opt/orca/orca
```

IRC 两端分别确认连接预期反应物与产物；自旋交叉可能要另外分析。组合充分收敛的同方法频率结果时，复制 `projects/fe/examples/homolysis_recipe.json` 或 `projects/co/examples/binding_recipe.json` 到工作目录，将其中相对路径指向实际作业目录。运行：

```bash
python -m chemflow thermo private/my_recipe.json
```

程序逐项检查元素、电荷、温度、压力、溶剂设定、方法与计算结束状态，再报告 ΔG 和按粒子数变化的 1 atm → 1 M 理想标准态校正。合并多个构象与自旋布居需要另外处置。同物种 Fe–C 解离热力学与热均裂**动力学**不能等同；H₂ 的气相化学势、光激发和单电子转移参比不适用此全溶质计算。对证据充足的单一基元步骤，可在给定标准态和传输系数时按 Eyring 关系估算速率，再用实测约束。有效池模型中合并的若干步骤不应直接赋予某一个 DFT 势垒。各场景速率储存在 `network.json`；不自动覆盖。

## 4　用已发表坐标做钴体系参照

取得 [Shimozono 补充材料](https://doi.org/10.1021/jacs.5c05965.s001)，保存为 `private/Shimozono_SI.pdf`。另需安装 Poppler 的 `pdftotext`。导入器已用作者原始 PDF 校验，可提取两个芳烃的四个文献鞍点。示例首先在气相用文献方法重做优化与频率，再在优化后几何做苯单点：

```bash
python tools/import_co_reference.py private/Shimozono_SI.pdf runs/reference
python -m chemflow prepare co ref_m_ts1 --xyz runs/reference/m-TS1.xyz --stage ts --method literature --reference-solvent gas --temperature 298.15 --out runs/co-ref-m
python -m chemflow run-orca runs/co-ref-m --exe /opt/orca/orca
python -m chemflow prepare co ref_m_ts1 --xyz runs/co-ref-m/job.xyz --stage sp --method literature --reference-solvent benzene --out runs/co-ref-m-sp
```

温度 298.15 K 是示范复核设定，重现时核对原作者的热处理设置。按同一路线对 p-TS1、m-TS3、p-TS3 计算。组合能量用苯单点 (E) 加气相频率任务的 (G-E)，相同底物分别取 para 减 meta，与 `co/benchmarks.csv` 的数值比。参照物与游离藜芦醚、目标溶剂和 50 °C 条件不同。铁的文献参照可用 `fe_methyl`、`fe_ethyl`、`fe_neopentyl` 和 `--reference-solvent benzene` 发起，须与实测热均裂活化焓匹配。

软件文档：[CREST](https://crest-lab.github.io/crest-docs/) · [ORCA 过渡态与 IRC](https://www.faccts.de/docs/orca/6.1/tutorials/react/irc.html) · [热力学输出](https://www.faccts.de/docs/orca/6.1/tutorials/prop/thermo.html) · [连续介质设置](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html)。

