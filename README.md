# 钴催化硼化与铁卟啉偶联的计算–实验工作流

把有限的合作实验优先放在会改变研究决策的对照上。钴项目关注游离藜芦醚在 HBpin 单一外加硼源、无烯烃受体时能否周转。铁项目关注二级苄基盐和特戊酸的偶联能否发生，以及低效源于前体损失、铁端过程还是竞争反应。

| 任务 | 入口 | 结果 |
|---|---|---|
| 先走一遍，不需商业软件或集群 | `python -m chemflow demo --out runs/first` | 两套 ORCA 输入、24 个可导入 COPASI 的 SBML、模型时间曲线和两份首轮实验清单 |
| 计算具体结构 | [量子化学操作](docs/quantum.md) | 构象、不同自旋的输入、计算结果检查、配平反应的自由能 |
| 由实测选择下一次合作实验 | [数据接入与实验决策](docs/data.md) | 守恒检查、留出条件的模型比较、配对区间、下一批实验 CSV |
| 确定做哪条支路 | [钴项目](projects/co/README.md) · [铁项目](projects/fe/README.md) | 各自的阳性、条件、关键观测、前进或停止条件 |
| 核对哪些内容已运行 | [交付核验](docs/verification.md) | 实跑软件版本、验证范围及待办 |

## 第一次运行

在解压后的目录打开终端，目录中应能看到本文件和 `chemflow`。推荐 Python 3.12。只跑示例、数据分析和输入生成时，Windows PowerShell 可用下面四行：

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m chemflow doctor
.\.venv\Scripts\python.exe -m chemflow demo --out runs/first
```

Linux、macOS 或 WSL 可用：

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python -m chemflow demo --out runs/first
```

有 conda/mamba 的 Linux 或 WSL 用户也可以使用下列环境，一并安装 CREST、xTB 和重建几何所需的 RDKit。脚本随包附带起始 XYZ，常规使用无需重新生成。安装 ORCA 后再运行量子化学任务。COPASI 图形界面按需另装。

```bash
conda env create -f environment.yml
conda activate catalysis-workflows
python -m chemflow doctor
python -m chemflow demo --out runs/conda-first
```

无需安装即可查看 [钴项目首批清单](projects/co/first_batch.csv) 和 [铁项目首批清单](projects/fe/first_batch.csv)。运行后打开 `runs/first/co/wet/next_batch.csv` 与 `runs/first/fe/wet/next_batch.csv`。两份表默认只显示核实母体/前体、分析方法和阳性的第一步。打开任一项目的 `scenarios/diagnostics.csv`，可见同一干预对不同假说产生的预测差异。**起始速率和噪声是假设值。这些数只能演示怎样做决策，不能写成研究结果。** `projects/*/benchmarks.csv` 单独保存有出处的文献记录。两项目的目标收率均留空。

在 COPASI 用 **File → Import SBML** 打开 `scenarios/*.xml`，进入 **Tasks → Time Course**，把 Duration 设为钴 48、铁 24，时间单位为小时，然后运行。CSV 中的浓度为 mol/L。图形界面与 Python 使用同一模型定义。实际验证过的跨软件对照见交付核验。

## 如何上传 GitHub

以本文件所在目录作为仓库根目录，上传解压后的源码和项目文件即可。项目包含自建起始结构、两个实跑的有机底物构象集合、示例数据以及文献/软件链接。个人计算输出放 `runs/`，未公开的合作资料放 `private/`，两者已列入 `.gitignore`。压缩包不附 ORCA 二进制、商业许可证或期刊全文。代码许可由仓库所有者自行决定。

软件入口在 [CREST](https://crest-lab.github.io/crest-docs/)、[xTB](https://xtb-docs.readthedocs.io/)、[ORCA](https://www.faccts.de/docs/orca/6.1/manual/) 和 [COPASI](https://copasi.org/Support/User_Manual/Importing_and_Exporting/Importing_and_Exporting_SBML_Files/)。

