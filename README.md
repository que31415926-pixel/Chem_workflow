# ChemFlow

面向催化研究的可复用计算与实验工作流，连接分子建模、量子化学、反应动力学、数据分析与实验决策。

通过 JSON 配置定义研究课题、计算步骤及实验判断规则，无需为每个新反应重新编写整套流程。

## 快速开始

推荐 Python 3.12。在仓库根目录运行：

```bash
git clone https://github.com/que31415926-pixel/Chem_workflow.git
cd Chem_workflow
python -m pip install -r requirements.txt
python -m chemflow doctor
python -m chemflow demo --out runs/first
```

建议在虚拟环境中安装依赖。Windows 如使用 `py -3.12`，可将上述 `python` 替换为 `py -3.12`。

默认依赖为 NumPy、SciPy 和 python-libSBML。内置两个课题的演示无需 ORCA、CREST 或集群，生成 ORCA 输入、动力学轨迹、24 个 SBML 模型及实验清单。新增课题后，`demo` 会执行所有已发现课题的配置步骤，包括已启用的外部计算。

需要 RDKit、CREST 和 xTB 时，Linux / WSL 用户可使用：

```bash
conda env create -f environment.yml
conda activate catalysis-workflows
python -m chemflow doctor
```

实际 DFT 计算需要另外安装获授权的 ORCA。

## 研究项目

| 项目 | 研究问题 | 入口 |
| --- | --- | --- |
| Co | 钴催化游离藜芦醚脱氢硼化 | [项目说明](projects/co/README.md) · [首批实验](projects/co/first_batch.csv) |
| Fe | 铁催化二级苄基与三级烷基偶联 | [项目说明](projects/fe/README.md) · [首批实验](projects/fe/first_batch.csv) |

只运行指定项目：

```bash
python -m chemflow validate co
python -m chemflow run co --out runs/co
```

运行结果位于 `runs/` 下。默认演示的速率与观测值是示例假设或合成数据，并非实测研究结论；文献基准单独记录在各项目的 `benchmarks.csv`。

## 创建新课题

```bash
python -m chemflow init my_reaction
python -m chemflow plan my_reaction
python -m chemflow validate my_reaction
python -m chemflow run my_reaction --out runs/my_reaction
```

`init` 从独立的教学模板创建项目，也可使用 `init my_reaction --from co` 复用已有项目。修改项目目录中的配置、分子结构和数据即可。每次运行请指定尚不存在的输出目录。

详见 [新课题配置与工作流复用](docs/reuse.md)。

## 仓库结构

```text
chemflow/           通用计算与实验决策引擎
projects/
  co/               钴催化研究
  fe/               铁催化研究
templates/starter/  新课题模板
docs/               使用指南与验证记录
tests/              回归测试
tools/              结构重建、文献坐标导入与独立复核
requirements.txt    基础 Python 依赖
environment.yml     含构象搜索工具的 Conda 环境
```

计算输出放在 `runs/`，未公开合作资料放在 `private/`，均已加入 `.gitignore`。

## 使用指南

| 文档 | 内容 |
| --- | --- |
| [量子化学](docs/quantum.md) | 分子结构、构象、自旋、ORCA 计算及自由能 |
| [数据与实验决策](docs/data.md) | 数据格式、动力学拟合、物料衡算与实验交接 |
| [新课题配置](docs/reuse.md) | 项目配置、步骤组合与中间结果传递 |
| [验证记录](docs/verification.md) | 已执行测试、跨软件复核及适用边界 |

查看全部命令：

```bash
python -m chemflow --help
```

## 验证

```bash
python -m unittest discover -s tests -v
```

修正版交付记录包含 11 项回归测试、CREST 实际构象搜索和 SciPy/COPASI 数值交叉核验。ORCA 输入生成与结果审核逻辑已验证，但尚未在授权环境完成新的 DFT 计算。详见 [验证记录](docs/verification.md)。

外部软件：[CREST](https://crest-lab.github.io/crest-docs/) · [xTB](https://xtb-docs.readthedocs.io/) · [ORCA](https://www.faccts.de/docs/orca/6.1/manual/) · [COPASI](https://copasi.org/)
