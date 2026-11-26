# 单账号私有计算环境

准备一台有 Python 3.10+ 的笔记本、Google Cloud CLI、一个独立的 Google Cloud 项目，以及获准在该项目的私有构建与计算环境使用的 ORCA 6 Linux x86-64 OpenMPI 4.1.6 官方安装包。项目应由账户持有人管理，已按提供商政策开通试用且具备所选区域配额。首次执行时复用 gcloud 当前账号，若尚未登录则启动浏览器授权。建议先用 `gcloud auth list` 检查当前账号。若使用 Google Cloud Shell，可在私有会话中运行相同命令。构建和 MPI 均在云端执行。

## 查看资源与预算

在引擎目录运行以下命令，将 `my-trial-project` 替换为自己的项目 ID。该命令离线打印区域、专用桶名、服务身份、Cloud Build 机型与最长构建时长。

```sh
python3 -m chemflow2.cloud_setup my-trial-project --region us-central1
```

将 examples/budget-policy.example.json 复制到自己的私有目录，在 budget-policy.json 填入上述项目、桶和区域以及经核实的试用状态、余额、有效期、保留时长、提供商价格和账单依据。构建的价格键为 `E2_HIGHCPU_8`，采用对应区域的 **Cloud Build** 小时上界。科学计算使用对应 Compute Engine 机型的上界。大小写有意义。存储、镜像、构建上传、失败与清理的费用应计入预算。申请一小时、100 GB 构建盘，账户和预算合格后才上传安装包。每次尝试使用独立预留，账单核对前保留失败尝试的预留。Cloud Billing 接口可核对项目有无账单关联，试用资格及剩余额度需要账户持有人提供实际依据。

## 私有构建与标定

```sh
python3 -m chemflow2.cloud_setup my-trial-project --region us-central1 \
  --profile "$HOME/chemflow-private/cloud.json" \
  --budget-policy "$HOME/chemflow-private/budget-policy.json" \
  --installer "$HOME/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418.run" \
  --calibration-plan "$HOME/chemflow-private/calibration/plan.json" \
  --byol-confirmed --apply
```

`--calibration-plan` 是可选参数，接受 compile-study 生成的工程标定计划。提供此项时，生成的 profile 只接受这个文件的 SHA-256。该计划请求的方法标签用于标定，须保留真实计算回执才能获得研究用资格。没有提供标定计划时，profile 的方法资格列表为空。生产 profile 使用已获相应方法资格的镜像摘要及节点配置。私有研究资料由协议输入传入。

安装器创建统一桶级访问且强制阻止公共读取的专用桶、私有镜像仓库和三个服务身份。构建身份可读取私有构建输入并写入项目专用仓库。作业身份可读取镜像、读写限定前缀下的新对象、写日志与报告 Batch 状态。编排身份可管理 Batch 作业、使用指定作业身份并读取作业区。项目需要授权创建上述资源的权限。启用 API 不会修改试用计划或自动升级付费。请在授权前核对专用项目的组织、文件夹和项目继承权限。

上传清单只包含引擎 Python 模块、上游 Python 模块、容器构建文件与指定的安装包。构建后删除暂存对象。构建前在新专用桶设置 private-build/ 前缀的一天删除策略，作为客户端意外关闭时的清理措施。云端安装标明的 OpenMPI、CREST、xTB、NumPy、SciPy、RDKit 等，在构建时检查架构、MPI 运行、ORCA 并行模块的动态链接和 Python 依赖。区域配额快照与最终私有镜像摘要写入 profile。实际方法、峰值资源和成本以研究计算与账单为准。

## 失败处理与拆除

在 profile 同目录保留以 `.resources.json` 结尾的资源清单。使用同一项目、区域和 profile 重试中断的初始化，已记录资源按云端标识核对后复用；若名称已被占用或其身份发生变化，命令会停止。在创建资源后、写清单前意外断电，可在项目活动日志中核对最后成功创建的资源，清理该单项后再重试。失败重建可以产生新的 Cloud Build 费用。每个研究完成后先执行 batch fetch，再用 batch cleanup 处理对应编排与作业对象。

确定已结束本项目中的 Cloud Build、Workflows 和 Batch 作业后，拆除专用镜像仓库及工作身份，保留桶和账单证据：

```sh
python3 -m chemflow2.cloud_setup my-trial-project --region us-central1 \
  --profile "$HOME/chemflow-private/cloud.json" --teardown
```

在下载并验证所有证据、预算账本、镜像与计费记录后，可在以上命令末尾附加 `--purge-bucket` 永久删除整个专用桶。资源清单记录实际移除的资源。保留的 Cloud Logging、账单导出及任何其他资源按提供商策略管理。安装包、账户记录、预算文件、运行结果和 profile 应置于私有目录。

[ORCA 官方安装参数](https://www.faccts.de/docs/orca/6.1/manual/contents/quickstartguide/installation.html) · [Google Workflows / Batch 示例](https://cloud.google.com/workflows/docs/tutorials/batch-and-workflows) · [Cloud Build 私有身份](https://cloud.google.com/build/docs/securing-builds/configure-user-specified-service-accounts)

## 从实测标定生成研究 profile

先取回标定结果，由研究审核人核对实际输入、原始 ORCA 日志、运行时间和资源记录。按方法、电子态与导数类别核对标定协议内预先分配的能力标签。独立从授权发行包核对 Linux ORCA 可执行文件本身的 SHA-256。在清理作业之前，执行：

```sh
python3 -m chemflow2.qualify "$HOME/chemflow-private/cloud.json" "$RUN" \
  "$HOME/chemflow-private/qualified.json" --reviewer "审核人姓名" \
  --authorized-orca-sha256 "$VERIFIED_BINARY_SHA256"
```

此入口通过登录的 Google 账号重新取回完整作业和冻结输入。它要求对应 Workflows 成功、全部 Batch 作业有提供商生成的 UID 且使用同一镜像与服务账号、每个 ORCA 任务来自同一授权可执行文件并通过该阶段的日志和产物核验。研究 profile 将每个标签映射到标定任务、回执与日志摘要并保留审核人。该记录仅支持已标定计算能力；科学闸门仍需各项原始基准数据。
