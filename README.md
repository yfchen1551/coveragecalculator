# 探针质量检测工具

本工具根据目标序列与探针序列的 BLAST 比对结果，计算探针在目标序列上的覆盖率，并生成探针覆盖情况热图。当计算覆盖率时，需要用户输入表格格式（以 '.csv' 为扩展名）的 BLAST 比对结果文件；而生成探针覆盖热图时，需要用户输入 XML 格式（以 '.xml' 为扩展名）的 BLAST 比对结果文件。该工具可以适用于任何情况下对目标序列覆盖率及覆盖情况的计算和检测，只要提供相应的 BLAST 输出文件。

## 功能

* **计算序列覆盖率**：统计 BLAST 比对结果中，每条 subject sequence 被 query 覆盖的总覆盖率，用该序列被覆盖序列的碱基数/该序列总碱基数，最大值为 1 ，并将结果输出到指定 '.csv' 格式文件。
* **生成覆盖情况热图**：按照覆盖率大小将所有序列由高到底排列，根据指定的排名范围（一次最多显示 30 条序列），利用 Python 中的 Matplotlib 生成探针覆盖情况的 PNG 图像。

## 依赖

本工具依赖于以下 Python 版本及包：

- [Python](https://www.python.org) &gt;= 3.8
- [NumPy]
- [Pandas]
- [IntervalTree]
- [BioPython]
- [Seaborn]
- [Matplotlib]

可以使用 conda 或 pip 安装上述包，例如：
```bash
conda install pandas
```
或
```bash
pip install pandas
```
这里推荐使用 conda 安装，配合清华镜像源可以有效提高安装速度。

本工具在以上环境下开发并测试，在更旧的 Python 环境中或许也可以运行，但是没有经过测试。

## 使用方法

### BLAST 结果文件的获取（Linux系统）

* **安装最新版BLAST**：
1. 下载：
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.14.0+-x64-linux.tar.gz
```
2. 解压：
```bash
tar -zxvf ncbi-blast-2.14.0+-x64-linux.tar.gz
```
3. 移动至目标目录：
```bash
mv ncbi-blast-2.14.0+ 目标目录
cd 目标目录路径
mv ncbi-blast-2.14.0+ blast
```
4. 添加环境变量：
```bash
export PATH=目标目录路径: $PATH
source ~/.bashrc
```
5. 验证是否安装成功
```bash
blastn -version
```
若显示版本则安装成功，否则安装失败，检查安装的过程中有没有出现问题。

* **建库**：
```bash
makeblastdb -in <input.fasta> -dbtype nucl -parse_seqids -out <out_db>
```
* **比对**：

1. **'.csv'格式输出结果**：
```bash
blastn -db <out_db> -query <query.fasta> -word_size n -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen' -out <out.csv>
```
输出的csv文件中需要包括query序列和目标序列的id、长度，以及覆盖的起始区域和终止区域等信息，否则将无法计算覆盖率。

2. **'.xml'格式输出结果**:
```bash
blastn -db <out_db> -query <query.fasta> -word_size n -outfmt 5 -out <out.xml>
```

### 帮助及参数查看
如果需要知道该工具的具体使用方法和详细细节，请运行：
```bash
python3 blast_coverage_heatmap.py -h
```
工具会提示用法：
```bash
python3 blast_coverage_heatmap.py [-h] [-c output_file] [-m start_rank end_rank output_file] blast_output
```
其中，output_file为输出文件名称及路径； start_rank 和 end_rank 分别代表序列起始排名和终止排名，在热图中最多显示 30 条序列的覆盖情况； blast_output 表示 BLAST 输出文件，在计算覆盖率时要求为 csv 格式，生成热图时要求为 xml 格式。

以下是本工具的重要参数：
* `-c/--coverage`：计算覆盖率并将结果输出到指定文件。保证输入的 BLAST 输出文件为表格格式（其文件扩展名应为 '.csv'）。如果不是，则输出错误信息并退出程序。
* `-m/--heatmap`：生成覆盖率热图。保证输入的 BLAST 输出文件为 XML 格式（其文件扩展名为 '.xml'）。如果不是，则输出错误信息并退出程序。

### 覆盖率计算
计算覆盖率时，输入以下命令：
```bash
python3 blast_coverage_heatmap.py -c output_file blast_output.csv
```
将会得到 BLAST 结果中所有有比对结果的序列名称及其对应的覆盖率，覆盖率在 0-1 之间。

### 热图生成
生成热图时，输入以下命令：
```bash
python3 blast_coverage_heatmap.py -m start_rank end_rank output_rile blast_output.png
```
将会得到 BLAST 结果中，指定区间内的 subject sequences 覆盖情况热图，一次最多显示 30 条序列。

