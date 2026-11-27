# 单账号额度准入与标定

在本地 profile 中设置 `budget_policy`，值为基于该 profile 目录解析的私有 JSON 路径。字段见 `examples/budget-policy.example.json`。将服务商账号状态、余额、试用条件、扣款限制、适用区域价目和用量上限依据保存在政策文件同目录，登记各自 SHA-256。保留原件，勿上传账号记录。确认账号处于供应商会阻止自费扣款的状态，且仅此项目使用额度。软件核对字段、时间、适用范围、文件完整性与数值；账号条款与原件内容由账号持有人核实。没有确认的价格和余额不能填零。

对每个 Batch 资源组，使用同时或依次可能产生费用的 VM 槽位总数、最大运行时长、启动和释放预留时长、整机小时价、盘容量与每 GB 小时价计算上界。每个资源组保留完整任务 ID。额外分别填写整次研究的证据存储、私有镜像、日志、网络、操作（含 Workflows、存储 API、云构建等未计价项目）、失败开销与清理费用上限。须覆盖证据下载、保留周期及最终删除；填零须有证据支持。币种为美元，账本逐次向上取整到分。单价、附加上限和额度不采用内置示意价格。

`observed_at` 是账号复核时刻，`valid_until` 取账号状态、授信及价格依据可确认的最早到期时刻。`tail_hours` 包含结束后的存储保留和清理窗口。引擎从有效期扣除运行上限、VM 启停预留及尾段时间，得到最晚派发时刻。队列等待和供应商实际费用存在不确定性。用户需维持已确认的账户级扣费限制；策略文件和软件预算不能替代该限制。

提交前把整张冻结任务图的上限预留到该账号唯一的私有桶 `chemflow-budget/v1/ledger.json`。读取精确对象版本，以 GCS 对象 generation 条件写入。发生竞争、权限故障或余额不足即停止。所有该账户的 ChemFlow 客户端须共用项目与桶。应禁止工作流及作业服务身份修改预算账本，并将清理限定为单次运行前缀。未提交成功、终止或执行失败的任务保留原预留，避免云账务延迟导致二次授信。若额度受阻，应先由账号持有人对照最终账单审计。不得通过改桶或删除账本重置额度。

## 实际资源与费用核对

选择真实研究矩阵中的高阶单点、频率和铁开壳层或 TS 作业进行代表性标定，分别保存失败及成功回执。Batch 作业的 ORCA/CREST worker 会记录资源采样。其他 Linux worker 可通过 `CHEMFLOW_MEASURE=1` 启用相同记录。回执包含引擎版本（ORCA 日志可识别时）、输入及产物摘要、墙钟与工作进程 CPU 时间。Linux cgroup v2 可用时记录内存组的生存期峰值和运行区间的采样高值，另记录私有目录已分配块使用量的采样高值。共享 cgroup 会包含同组其他进程，每 5 秒采样可能错过瞬时 scratch 峰值。缺少计量器的数据保留为空。用独占作业与供应商指标交叉校验后再确认峰值。

使用既有 `estimate` 将经审核可比的 pilot 回执映射到完整任务图。另可直接审查费用上限或核对实账：

```sh
python3 -m chemflow2.budget quote profile.json allocations.json
python3 -m chemflow2.budget reconcile private/billing-manifest.json
```

`allocations.json` 可直接使用 batch plan 的 JSON 输出，也接受数组。每行提供 `id, task_ids, machine_type, max_instances, max_hours, disk_gb`。对账清单包含 `project, run_id, billing_export, receipts`。其中 billing_export 为 `{"file":"billing-lines.csv","sha256":"..."}`，receipts 是同格式的回执引用数组。CSV 采用 `line_id,project,run_id,currency,gross_usd,credit_usd`；从完整提供商账单及资源标签提取归属于该次研究的费用，使用唯一行 ID，优惠额度填正数。结果给出观测任务、导出记录毛额、抵扣额及净额。导出的完整性、共享成本归属及账单最终性须独立复核。没有账单，无法计算一个账号最终覆盖的实验数。

## 结算已结束运行的预留

预留在运行结束后不会自动释放。账号持有人对照控制台账单（试用赠金余额、当月费用）确认实际花费之后，可以把**已经结束**的运行改按审计后的金额结算：

```sh
python3 -m chemflow2.budget settle profile.json build-xxxx=3.00 cf-xxxx=3.00 \
  --approved-by "账号持有人" --evidence "控制台：赠金剩余 USD 298.71，当月费用 0.00（2026-11-26 22:10）" \
  --backup-dir private/budget_audit
```

- 命令会查询各运行最近一次工作流执行的状态，**仍在运行的一律拒绝结算**；云构建 `build-*` 视为已结束，从未部署的运行也视为已结束。
- 结算金额必须为正，且不超过原预留。已经结算过的运行不会再结算第二次。
- 账本条目一律不删除，每条都保留原上限 `original_upper_usd` 和审计信息 `settlement_audit`。
- 写回账本使用 generation 条件，结算前的版本备份到 `--backup-dir`。
