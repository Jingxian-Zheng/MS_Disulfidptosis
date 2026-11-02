# GitHub 上传说明

建议仓库名：`ms-stromal-transcriptomics`。

## 需要上传的内容

`ms-stromal-transcriptomics_github_v1.0.0.zip` 是最终代码上传版。解压后，上传其中
`ms-stromal-transcriptomics` 文件夹里的文件和子目录，保持目录结构。不要只上传这个 ZIP
而让仓库里只有一个压缩文件。

当前目录也已提供同名的完整文件夹，可直接使用。仓库内有英文和中文 README、MIT 许可、
全部正式分析代码、历史科学分析、第三方许可、依赖记录、固定参数、40 个结果分表和验证程序。

## 推荐：GitHub Desktop

本机同名代码文件夹已初始化为本地 Git 仓库。在 GitHub Desktop 添加该本地仓库，检查待提交
文件，填写首次提交说明并提交，然后使用 Publish repository 发布到你的 GitHub 账户。
若使用另行解压的 ZIP，可先用 GitHub Desktop 为该文件夹创建本地仓库，再提交和发布。

账号、远程仓库和发布权限由你在 GitHub 中选择。本次交付未向任何远程仓库上传。
官方操作说明：[向 GitHub 添加本地代码](https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github)。

## 使用 GitHub 网页上传

创建空仓库后，使用 Add file → Upload files，分批上传解压后的文件，包含 `.github`、
`.gitignore` 和 `.gitattributes`。仓库文件超过 100 个，网页上传需要分批操作。
已有 README 和 LICENSE，创建远程仓库时无需再生成另一份。
官方限制与操作见：[上传文件](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository)。

## 配套数据与本地文稿工具

- `S1_Data.zip` 是已有修订分析的数据包，约 92 MB。它放在代码文件夹外，不会混入代码压缩包。
  可将它作为 Release 附件或保存到数据仓库，公开后在代码仓库记录真实下载地址。
  网页普通文件上传有 25 MiB 的单文件限制，见 [GitHub 大文件说明](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)。
- `manuscript_tools_local_only.zip` 保留 6 个文稿与审稿回复生成工具，供本地使用。
  这些工具含未发表审稿回复的文字载荷，未放入公共代码包，不影响分析复现。
- `release_audit/` 为本次整理核查记录；分析结果和原始代码未被覆盖。

## 他人如何运行

进入公开代码仓库目录后：

```text
python run.py check
python run.py list
python run.py prepare --archive 配套数据包的实际路径
Rscript environment/install_dependencies.R
python run.py run validate
python run.py run figures
```

完整说明在 `docs/REPRODUCIBILITY.md`；默认复现工作目录为 `work/reproduction/`，已被 Git 忽略。
已在本机实际完成数据解包校验、17 项关键统计复核和全部 10 张统计图重建。
全新电脑的依赖安装以及全部上游重训练不在本次已执行核查范围内。

## 与审稿回复衔接

实际公开后再填写真实仓库链接；GitHub 仓库并不会自动产生 DOI。若仍需要永久 DOI，须在
实际完成数据/代码归档后使用归档平台给出的编号。本地文件不应写成已经公开。
