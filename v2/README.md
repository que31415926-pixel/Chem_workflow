# ChemFlow 2.0 dev3

笔记本编译研究协议，云端按依赖运行计算，结果可按任务取回并校验。客户端需要 Python 3.10+。从仓库根目录进入 v2 后运行以下命令。该示例只在本地生成输入和调度计划。

```sh
cd v2
python3 -m chemflow2 --version
python3 -m chemflow2 compile-study examples/water-calibration/study.json water-inputs
python3 -m chemflow2 validate water-inputs/plan.json
python3 -m chemflow2 batch plan examples/batch-profile.example.json water-inputs/plan.json
```

可在虚拟环境中用 `python3 -m pip install ./v2` 从仓库根目录安装命令行工具。安装后可从任意目录运行 `chemflow2 --help`。本地执行动力学积分时可安装 `./v2[kinetics]`。

## 运行自己的研究

1. 使用研究协议登记原件、经审核结构、方法矩阵、依赖及判据。research inventory 展开研究单元，research assemble 核对已提交和缺失的输入，research compile 生成 worker 任务。详见 [研究协议](docs/research-packs.md)。
2. 使用 cloud-setup 为自己的 Google Cloud 项目准备私有执行环境。需要 Google Cloud CLI、该项目的授权和符合私有云使用条件的 ORCA 安装包。先准备账户的实际额度和计价依据，然后执行代表性方法标定。详见 [私有云初始化](docs/cloud-setup.md) 和 [预算与对账](docs/budget.md)。
3. 使用 batch submit 提交。Workflows 在客户端退出后继续编排，计算任务向私有桶上传回执与产物。使用 batch status、fetch、cancel 和 cleanup 管理此次研究。详见 [任务生命周期](docs/batch-workflows.md)。

ORCA 路径支持经审核初猜、过渡态搜索、虚频、双向 IRC 和端点复核。离子对工作流调用 CREST 搜索显式构型。区间动力学按质量守恒、热力学约束及训练端点筛选可行域，并输出留出预测。方法与数据约定见 [反应路径与溶液计算](docs/pathway-and-solvation.md) 和 [区间动力学](docs/interval-kinetics.md)。

## 验证范围

从本目录运行 `python3 -m unittest discover -s tests`。SciPy 用于积分测试。真实 ORCA 生产验证需要用户授权的目标云账号、软件和科研输入。账户能否用免费额度完成一次完整研究，须以实际工作量、供应商计费记录和清理记录核对。

新引擎的软件源代码位于 chemflow2，云镜像通过私有构建安装科学程序。vendor/chemflow-v1/chemflow 是同一仓库第一代引擎的固定算子副本，其中包含兼容性修订。构建依赖由各软件提供方分发。公开目录不含研究原件、未发表的两篇论文协议或 ORCA 安装包。
