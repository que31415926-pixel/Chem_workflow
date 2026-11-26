# ORCA 6.1.0 与工具链版本锁定（多设备统一）

所有为论文产出数字的设备（GCP 镜像、组里的 Mac、合作者的 Linux 工作站）都必须使用
[`toolchain.lock.json`](../toolchain.lock.json) 中的同一套版本。不同 ORCA 小版本或不同 xTB
版本得到的能量可能相差零点几 kcal/mol，混用会让不同设备算出的数据无法直接比较。

| 组件 | 锁定版本 |
| --- | --- |
| ORCA | 6.1.0（Linux 安装包 `orca_6_1_0_linux_x86-64_shared_openmpi418.run`） |
| Open MPI | 4.1.8 |
| CREST / xTB / tblite | 3.0.2 / 6.8.0 / 0.8.0 |
| Python / RDKit / NumPy / SciPy | 3.12.16 / 2026.09.2 / 2.5.5 / 1.18.1 |
| python-libsbml / google-cloud-storage | 5.20.5 / 3.17.2 |

这些版本来自 2026-11-26 GCP Cloud Build `b1230b1b` 实际解析出的环境。`v2/deploy/byol.Dockerfile`、
`environment.yml` 与锁文件逐项一致，`v2/tests/test_toolchain_lock.py` 会检查三者是否一致。

## 1. 下载 ORCA（只能手动）

ORCA 学术免费，但许可证不允许再分发，因此仓库里**没有也不能提交**安装包或二进制文件。

1. 在 <https://orcaforum.kofo.mpg.de> 注册并登录，进入 Downloads → ORCA 6.1.0。
2. 按平台下载：
   - **Linux x86-64（GCP 与集群）**：`orca_6_1_0_linux_x86-64_shared_openmpi418.run`，
     大小 583,759,473 字节，SHA-256 `3df2d25c36fe7e49f1c301b3ad6cfc9b837a370fd3cee31275a86ecd291ca944`。
   - **macOS 或其他平台**：下载对应平台的 **6.1.0** 版本。不同平台的二进制哈希本来就不同，
     所以只强制检查版本号和 MPI 能否运行。
3. 不要下载 6.1.1 或 6.0.x。不同设备之间只允许同一个版本号。

## 2. 校验与安装

```bash
# Linux x86-64：校验文件名、大小和 SHA-256，然后安装到 ~/orca_6_1_0
bash tools/toolchain/install_orca.sh ~/Downloads/orca_6_1_0_linux_x86-64_shared_openmpi418.run ~/orca_6_1_0

# macOS：同一个命令只校验 Linux 安装包（这是给 GCP 用的那个文件），不会安装；
# Mac 本机请自行安装 6.1.0 的 macOS 版本
```

在 Linux 上，ORCA 并行需要 Open MPI 4.1.8 在 `PATH` 和 `LD_LIBRARY_PATH` 上。
没有的话，可以按锁文件中的 `openmpi.source_url` 编译，方式与 Dockerfile 相同。
容器里没有 ssh 时，还需要设置 `OMPI_MCA_plm=isolated`（镜像中已经设好）。

## 3. 每台设备都要跑一次检查

```bash
python3 tools/toolchain/check_toolchain.py --orca /绝对路径/orca
```

这一步会跑两个水分子 HF 小任务（单进程和 2 进程 MPI），然后逐项比对 ORCA 版本、
Linux 二进制哈希、Open MPI、CREST、xTB 和 Python 包的版本，最后写出 `toolchain_receipt.json`。
任何 `FAIL` 都要先修好再产出结果。只做 DFT 的设备没装 CREST/xTB 时会显示 `WARN`，这是允许的。
把 receipt 和计算结果一起归档，就能追溯每个数字是在哪台设备、用什么版本算出来的。

## 4. GCP 镜像

`python -m chemflow2.cloud_setup` 只接受锁文件中的 Linux 安装包文件名（openmpi418）。
`chemflow2.qualify --authorized-orca-sha256` 使用的是**安装后 orca 二进制**的哈希
`903987dbd651df922ef6a7f91763abb25e7a79531f01018f5aa49fea104f2e91`，不是安装包本身的哈希。
详见 [`v2/docs/cloud-setup.md`](../v2/docs/cloud-setup.md)。

## 5. 升级版本的规则

只有在全组同意后才能改版本，并且要在同一个提交里同时修改 `toolchain.lock.json`、
`v2/deploy/byol.Dockerfile`（及打包副本）和 `environment.yml`，然后重新构建镜像、重新跑 qualify，
最后在所有设备上重跑 `check_toolchain.py`。已经算出的数据不能和新版本的数据混在同一张表里。
