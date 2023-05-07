#!/home/usr/bin/python

import sys
import argparse
import numpy as np
import pandas as pd
import csv
import os
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from Bio.Blast import NCBIXML
import seaborn as sns
import matplotlib.pyplot as plt

# 处理BLAST输出结果并计算覆盖率
def process_blast_output(blast_output_file):
    # 使用IntervalTree存储覆盖区间
    coverage = defaultdict(IntervalTree)
    # 存储参考序列长度
    lengths = defaultdict(int)

    # 打开文件，逐行处理
    with open(blast_output_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row_num, row in enumerate(reader, start=1):
            try:
                # 获取BLAST结果中每一行的数据，提取query序列、参考序列（ref）、起始位点（start）、结束位点（end）和参考序列的长度（s_len）
                query, ref, start, end, s_len = row[0], row[1], int(row[8]), int(row[9]), int(row[12])
                # 创建新区间new_interval，表示当前行中的比对区间，使用Interval对象存储区间的起始和结束位点
                new_interval = Interval(min(start, end), max(start, end))
                
                '''
                查找在`coverage[ref]`（一个`IntervalTree`对象）中与`new_interval`重叠的所有区间，以便在后续代码中合并这些区间并更新覆盖范围。
                `coverage`是一个字典，其中的key是参考序列（ref）的ID，value是一个`IntervalTree`对象。`IntervalTree`对象用于存储已覆盖的区间。
                `new_interval`是一个`Interval`对象，表示当前BLAST比对结果中的起始和结束位置。`new_interval.begin`和`new_interval.end`分别表示这个区间的起始和结束位置。
                `coverage[ref].overlap(new_interval.begin, new_interval.end)`函数调用在`IntervalTree`中查找与给定范围（`new_interval.begin`至`new_interval.end`）重叠的所有区间。返回的`overlapping_intervals`是一个包含重叠区间的集合。
                '''
                overlapping_intervals = coverage[ref].overlap(new_interval.begin, new_interval.end)

                # 处理重叠区间
                '''
                如果没有重叠区间，直接将新区间添加到'IntervalTree'中。如果有重叠区间，则将所有的重叠区间与新区间合并。
                首先将所有的重叠区间放入一个集合'merged_intervals'中，并将新区间加入。
                然后使用'difference_update'函数从'IntervalTree'中移除所有重叠的区间，这样在计算coverage时就不会计算重叠区间，使结果更准确。
                计算合并后区间的最小起始位点和最大结束位点，创建一个新的'Interval'对象，将它添加到'IntervalTree'中。这个新区间是原来重叠区间合并的结果，包括当前比对结果。
                如果处理比对结果时发生索引错误，例如输入文件格式不正确，则打印提示，要求用户检查输入文件格式，退出程序并返回错误代码1。
                '''
                if not overlapping_intervals:
                    coverage[ref].add(new_interval)
                else:
                    merged_intervals = set(overlapping_intervals)
                    merged_intervals.add(new_interval)
                    coverage[ref].difference_update(overlapping_intervals)

                    start = min(interval.begin for interval in merged_intervals)
                    end = max(interval.end for interval in merged_intervals)
                    coverage[ref].add(Interval(start, end))

                lengths[ref] = s_len
            except IndexError:
                print(f"Error: IndexError occurred while processing row {row_num}. Please make sure the BLAST output file has the correct format.")
                sys.exit(1)

    #计算覆盖率
    '''在遍历完BLAST比对结果后，计算每个参考序列的覆盖率。
    对于每个参考序列，使用'length'字典查找参考序列的长度，对于'IntervalTree'中的每个区间，使用'interval.length()'计算该区间的长度，
    然后将所有的区间的长度相加，得到覆盖的碱基数目。
    最后将覆盖的碱基数目除以参考序列长度，得到覆盖率。
    '''
    coverage_ratio = {}
    for ref, length in lengths.items():
        covered_bases = sum(interval.length() for interval in coverage[ref])
        coverage_ratio[ref] = covered_bases / length

    return coverage_ratio, lengths

# 将每条参考序列的覆盖率输出的文件
def output_coverage(coverage_ratio, output_file):
    with open(output_file, 'w') as outfile:
        outfile.write("Ref_ID\tCoverage\n")
        for ref, cov in coverage_ratio.items():
            outfile.write(f"{ref}\t{cov:.6f}\n")

# 解析BLAST输出结果并计算热图覆盖率
'''
使用BioPython中的NCBIXML模块解析输入的BLAST比对结果的XML格式文件，将结果储存在'blast_records'变量中。
创建一个名为 coverage 的 defaultdict 对象，其默认值为另一个 defaultdict，而这个内部的 defaultdict 的默认值为整数 0。这样的嵌套结构使得在使用 coverage 时，如果某个键不存在，就会自动创建一个具有默认值的新条目。
用于迭代处理 BLAST 记录、比对和高分对等片段（HSP）的四个嵌套 for 循环：
for blast_record in blast_records：遍历 blast_records 中的每个 BLAST 记录。
for alignment in blast_record.alignments：遍历每个 BLAST 记录中的比对结果。
for hsp in alignment.hsps：遍历每个比对结果中的高分对等片段（HSP）。
for pos in range(hsp.sbjct_start, hsp.sbjct_end)：遍历 HSP 的每个覆盖到的目标序列的位置。
对于每个覆盖到的目标序列位置（pos），在 coverage 字典中找到相应的目标序列（alignment.title），并将该位置的覆盖计数加1。
最后，将计算得到的覆盖率字典 coverage 返回给调用者。
'''
def parse_blast_output_heatmap(blast_output_file):
    blast_records = NCBIXML.parse(blast_output_file)
    coverage = defaultdict(lambda: defaultdict(int))

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                for pos in range(hsp.sbjct_start, hsp.sbjct_end):
                    coverage[alignment.title][pos] += 1

    return coverage

# 覆盖率热图的可视化参数
def visualize_coverage_heatmap(coverage, output_file, start_rank, end_rank):
    df = pd.DataFrame(coverage).fillna(0).T
    df.sort_index(axis=1, inplace=True)

    # 根据总覆盖率和指定排名范围过滤序列（因为热图显示的比对结果有限）
    df['total_coverage'] = df.sum(axis=1)
    df.sort_values(by='total_coverage', ascending=False, inplace=True)
    df = df.iloc[start_rank-1:end_rank]
    del df['total_coverage']
    
    # 生成热图并保存为PNG格式文件
    fig, ax = plt.subplots(figsize=(14, 6))
    sns.heatmap(df, cmap='viridis', ax=ax)

    ax.set_xlabel("Probe Position")
    ax.set_ylabel("Target Sequences")
    ax.set_title("Coverage Depth Heatmap")

    plt.yticks(rotation=0)
    plt.subplots_adjust(bottom=0.25, top=0.90)
    plt.savefig(output_file, dpi=300)
    plt.close()

# 主函数
'''
解析命令行参数并设置程序运行所需的参数。
首先，创建一个 ArgumentParser 对象，并提供程序的描述信息。
添加一个位置参数（positional argument），命名为 blast_output，并提供相应的帮助信息。该参数用于指定输入的 BLAST 输出文件。
添加一个可选参数（optional argument），短选项为 -c，长选项为 --coverage，并提供相应的帮助信息。该参数用于指定输出覆盖率结果的文件名，并要求输入的 BLAST 输出文件为表格格式。
添加另一个可选参数，短选项为 -m，长选项为 --heatmap，并提供相应的帮助信息。该参数接受三个值（由'nargs=3'指定），分别为开始排名、结束排名和热图输出文件名。该参数用于生成覆盖率热图，并要求输入的 BLAST 输出文件为 XML 格式。
最后解析命令行参数，并将结果存储在 args 对象中。该对象可以通过其属性访问解析到的参数值。例如，可以使用 args.blast_output 访问 BLAST 输出文件名。
'''
def main():
    try:
        parser = argparse.ArgumentParser(description="Calculate coverage and visualize heatmap from BLAST output")
        parser.add_argument("blast_output", help="BLAST output file")
        parser.add_argument("-c", "--coverage", metavar="output_file", help="Output coverage results to a file (requires tabular BLAST output)")
        parser.add_argument("-m", "--heatmap", nargs=3, metavar=("start_rank", "end_rank", "output_file"), help="Output heatmap to a PNG file (requires XML BLAST output)")
        args = parser.parse_args(args=[])

    # 检查输入文件是否存在，如果不存在，打印错误提示，要求用户输入BLAST输出结果，退出并返回错误代码1
        if not os.path.isfile(args.blast_output):
            print(f"Error: Input file {args.blast_output} does not exist.")
            sys.exit(1)
    
        # 根据参数执行相应功能
        '''
        1. 如果用户输入了 -c 或 --coverage 选项及其参数，表示需要计算覆盖率并将结果输出到文件。
        检查输入的 BLAST 输出文件是否为表格格式（假设其文件扩展名为 '.csv'）。如果不是，则输出错误信息并退出程序。
        调用 process_blast_output 函数处理输入的表格 BLAST 输出文件，并返回覆盖率（coverage_ratio）和参考序列长度（lengths）。
        将覆盖率输出到用户指定的文件。
        最后输出一条消息，告知用户覆盖率计算已完成，结果已保存到指定文件。

        2. 如果用户提供了 -m 或 --heatmap 选项及其参数，表示需要生成覆盖率热图。
        首先检查输入的 BLAST 输出文件是否为 XML 格式（假设其文件扩展名为 '.xml'）。如果不是，则输出错误信息并退出程序。
        从用户提供的参数中获取开始排名、结束排名和热图输出文件名。
        以只读模式打开输入的 XML BLAST 输出文件。
        调用 parse_blast_output_heatmap 函数处理输入的 XML BLAST 输出文件，并返回覆盖率。
        根据给定的覆盖率、输出文件名和排名范围，生成覆盖率热图。
        最后输出一条消息，告知用户热图可视化已完成，结果已保存到指定文件。
        '''
        
        if args.coverage:
            if not args.blast_output.endswith('.csv'):  # Assuming tabular BLAST output has a '.csv' extension
                print("Error: Coverage calculation requires a tabular BLAST output file (e.g., file.csv).")
                sys.exit(1)
            coverage_ratio, lengths = process_blast_output(args.blast_output)
            output_coverage(coverage_ratio, args.coverage)
            print(f"Coverage calculations completed. Results saved to {args.coverage}")

        if args.heatmap:
            if not args.blast_output.endswith('.xml'):  # Assuming XML BLAST output has a '.xml' extension
                print("Error: Heatmap visualization requires an XML BLAST output file (e.g., file.xml).")
                sys.exit(1)
            start_rank, end_rank, output_file = int(args.heatmap[0]), int(args.heatmap[1]), args.heatmap[2]
            with open(args.blast_output, "r") as blast_output_file:
                coverage = parse_blast_output_heatmap(blast_output_file)
            visualize_coverage_heatmap(coverage, output_file, start_rank, end_rank)
            print(f"Heatmap visualization completed. Results saved to {output_file}")

    except Exception as e:
        print(f"Error: {str(e)}")
        print("Usage:")
        print("  python script.py blast_output [-c output_file] [-m start_rank end_rank output_file]")
        print("\nFor more help, use the -h or --help option.")
        sys.exit(1)

if __name__ == "__main__":
    main()
