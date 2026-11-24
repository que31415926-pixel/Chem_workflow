# 用同一套工作流接入新课题

安装方法见仓库首页。每个项目对应一个目录。命令按目录名发现项目。数值计算、软件调用、读数检查和阶段交接共用一套实现。现有项目可作为研究预设独立修改。

## 从空白课题开始

```bash
python -m chemflow init my_reaction
python -m chemflow projects
python -m chemflow plan my_reaction
python -m chemflow validate my_reaction
python -m chemflow run my_reaction --out runs/my_reaction
```

`init` 从独立的 E/Z-2-丁烯教学模板生成 `projects/my_reaction/`。默认执行输入制备、合成数据拟合、带入拟合参数的时间程模拟、物料衡算和首轮实验清单。首次运行结果应带有 `synthetic_fixture`。用自己的研究数据替换模板后再作科学判断。若已有相似项目，可用 `python -m chemflow init my_variant --from co` 复制其配置、示例与结构，另行修改。命令拒绝覆盖已有目录。

项目目录也可以在仓库外。把全局选项放在命令前面。例如 `python -m chemflow --projects-dir ../lab-projects init kinetics1`，后续的 `validate` 和 `run` 使用同一选项。也可设置 `CHEMFLOW_PROJECTS_DIR`。默认根目录为仓库内的 `projects/`。新项目 ID 支持字母、数字和下划线，首字符为字母或下划线；中文全名写入 `title`。

## 项目文件怎样分工

| 想改变的内容 | 修改位置 | 实例 |
|---|---|---|
| 用哪些步骤、先后顺序、数据交给哪一步 | `project.json → workflow` | 直接做实验交接可只保留 `["handoff"]`；新建目录只需这一项目清单、wet_plan.json 和 state.json |
| 温度、介质、自旋、量化方法 | `project.json → qm`、`species.json` | 每种分子定义净电荷和多重度；每个作业选底物、方法、自旋和介质。命名的连续介质或 ε/n 均可设置 |
| 换底物 | `species.json` 和 `xyz/` | 可导入已审核 XYZ；有机物也可配置 SMILES，执行 `build-xyz`。金属配位构型需提供合适的起始几何 |
| 反应步骤、条件、时间窗、基线、拟合参数 | `network.json` | `lhs/rhs` 写化学计量，`baseline` 指定参照条件；`fit` 可以是单一速率 ID 或多个速率 ID 的列表。只有一个假说也能积分并导出 SBML |
| 化学计量与衡算口径 | `material.json` | 每个输出由分子和分母的字段加权和、比例系数组成。空白保留为未测，可单独设置含两份片段的二聚物权重 |
| 下一批最值得做的实验 | `wet_plan.json`、`state.json` | 自定状态字段、可选值、阶段、实验及读数。`rules` 从上到下取第一个匹配；末项用 `when: {}` 作兜底；一个字段可匹配多个允许值 |
| 百分产率的配对判定 | `project.json → paired` | 自定 `practical`、`gap`、`min_blocks`；运行 `paired 项目名 数据.csv`。没有这种决策需求的项目可以省略 |

## 步骤怎样把结果传下去

以下是新建模板已能运行的一段配置。相对输入路径从项目目录解析。`@estimate/calibrated.json` 从本次运行的 `estimate` 步骤取文件，因此改动数据并重新运行后，模拟会读取新的拟合结果。

```json
"workflow": [
  "prepare",
  {"step": "fit", "id": "estimate", "data": "examples/observations.csv"},
  {"step": "simulate", "id": "forecast", "calibrated": "@estimate/calibrated.json"},
  {"step": "balance", "id": "closure", "data": "examples/material.csv"},
  "handoff"
]
```

若需要从新分子的 SMILES 开始，可在已安装 RDKit 和 CREST 的环境中把工作流改为下列三个步骤。模板中的 E_butene 可直接用于此示例。换成其他分子时一起更新 `species.json` 的 SMILES、电荷、多重度和任务物种名。研究需要区分构象时，分别传入不同帧。运行前可用 `plan` 查看引用关系，`validate` 检查配置和本地已有输入。依赖上一步生成的文件在真正执行时检查。

```json
"workflow": [
  {"step": "build-xyz", "id": "seeds", "species": ["E_butene"]},
  {"step": "conformers", "id": "sample", "species": "E_butene", "xyz": "@seeds/E_butene.xyz", "quick": true},
  {"step": "prepare", "jobs": [{"species": "E_butene", "xyz": "@sample/crest_conformers.xyz", "frame": 1}]}
]
```

可以在最后加 `{"step":"run-orca","id":"energy","job":"@orca_input","exe":"/实际路径/orca"}` 与 `{"step":"audit","id":"check","job":"@energy"}`。多个制备作业在 `@orca_input/作业ID` 下分别取用。`thermo` 接收 `recipe` 对配平作业求自由能。运行时 `--state`、`--calibrated` 可覆盖对应步骤的默认输入；新结果请指定新的输出目录。软件名、溶剂关键字和所选方法仍需受本机安装版本支持。

## 适用范围与边界

目前提供 10 类可组合步骤，包括有机结构生成、构象搜索、ORCA 输入/运行/审核、热力学汇总、动力学拟合/模拟、物料衡算和湿实验交接。`simulate` 实现以小时、mol/L 为单位的质量作用动力学，`paired` 实现配对百分产率的区间判定。新增底物、网络、条件或判断阈值可改配置；若要加入另一类求解器、任意非质量作用速率律或全新的统计设计，才需要扩展相应模块。拟合多个速率时应单独检验可辨识性；程序允许拟合不等于现有数据能区分它们。

钴、铁专属的金属起猜重建及文献坐标提取脚本用于复现这两个预设，扩展新课题无需修改。合作者提供的未公开输入可放在独立的项目目录。原有 `simulate co`、`handoff fe`、`paired-fe` 等命令继续可用。

