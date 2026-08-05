# lefse-format_input.py脚本将输入文件lefse.txt转换成LEfSe分析所需的内部格式，并保存为input.in
# -c指定分组行；-s指定亚组行，若没有可以不指定；-u指定样本编号行。-o 1000000指定了归一化后范围
lefse_format_input.py lefse_example_standard.txt temp/input.in -c 1 -u 2 -o 1000000
# lefse_run.py脚本执行LEfSe分析，输入文件是input.in，输出结果保存为input.res。
# 如果样本有重复或分组信息有误，分析可能会报错。
# -l 指定LDA score阈值，默认为2，结果太多可调整为3或者4。
lefse_run.py temp/input.1m.in temp/input.1m.res
# lefse-plot_cladogram.py脚本根据分析结果temp/input.res绘制物种树，并以PDF格式保存到指定路径。
# 这个脚本有很多参数可供调整，比如字体、水平、颜色、显著性、图例等。
lefse_plot_cladogram.py \
  temp/input.1m.res \
  temp/lefse_cladogram.pdf \
  --format pdf \
  --dpi 300
# lefse-plot_res.py脚本绘制所有差异显著的特征（features）的柱状图，并以PDF格式保存。
lefse_plot_res.py \
  temp/input.1m.res \
  temp/lefse_LDA.pdf \
  --format pdf \
  --dpi 300 \
  --width 10 \
  --feature_font_size 8
# 查看显著差异物种，按丰度排序
grep -v '-' input.res | sort -k3,3n 
# 使用lefse-plot_features.py脚本绘制特定特征（例如Proteobacteria）的柱状图。
# -f one表示只绘制一个特征。
# --feature_name后面跟的是物种的名称。
mkdir -p temp/lefse_features

lefse_plot_features.py \
  -f diff \
  --archive none \
  --format pdf \
  --dpi 300 \
  temp/input.1m.in \
  temp/input.1m.res \
  temp/lefse_features/