#!/bin/bash

# 设置输入数据文件和结果文件的路径
input_file="/home/jh/code/JHPC/dataset/Yeast/PPI/Collins.txt"
result_dir="/home/jh/code/JHPC/result"

# 设置参数范围
min_threshold=0.1
max_threshold=1.0
step=0.1

# 生成1000个参数组合
count=0
for p1 in $(seq $min_threshold $step $max_threshold); do
    for p2 in $(seq $min_threshold $step $max_threshold); do
        for p3 in $(seq $min_threshold $step $max_threshold); do
            result_file="${result_dir}/collins_${p1}_${p2}_${p3}.txt"
            ./cns "$input_file" "$result_file" "$p1" "$p2" "$p3"
            let count=count+1
            if [ $count -eq 1000 ]; then
                break
            fi
        done
        if [ $count -eq 1000 ]; then
            break
        fi
    done
    if [ $count -eq 1000 ]; then
        break
    fi
done