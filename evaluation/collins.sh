# python evaluation/match.py -n dataset/Yeast/PPI/Collins.txt ./dataset/Yeast/complex/yeast_complex.txt ./result/collins_0.3_0.1_0.7.txt

!/bin/bash

ppi_file="/home/jh/code/JHPC/dataset/Yeast/PPI/Collins.txt"
ref_file="/home/jh/code/JHPC/dataset/Yeast/complex/yeast_complex.txt"
pre_dir="/home/jh/code/JHPC/result"

# 设置输入数据文件和结果文件的路径
# input_file="/home/jh/code/JHPC/dataset/Yeast/PPI/Collins.txt"
# result_dir="/home/jh/code/JHPC/result"

# 设置参数范围
min_threshold=0.1
max_threshold=1.0
step=0.1

# 生成1000个参数组合
count=0
for p1 in $(seq $min_threshold $step $max_threshold); do
	for p2 in $(seq $min_threshold $step $max_threshold); do
		for p3 in $(seq $min_threshold $step $max_threshold); do
			result_file="${pre_dir}/collins_${p1}_${p2}_${p3}.txt"
			# ./cns "$input_file" "$result_file" "$p1" "$p2" "$p3"
			python match.py -n "$ppi_file" "$ref_file" "$result_file"
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

