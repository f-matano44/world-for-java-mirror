#! /usr/bin/env python3

# usage
# ja-world/$ python3 test/compareF0.py

import csv
import pandas as pd
import numpy as np


# Dio ----------------------------------------------------------
with open("build/f0Dio_java.csv", 'r') as f:
    reader = csv.reader(f)
    data_cpp = [row for row in reader]
with open("debug/f0Dio_cpp.csv", 'r') as f:
    reader = csv.reader(f)
    data_java = [row for row in reader]

Error = 0.0
maxError = 0.0
for i in range(1, len(data_cpp)):
    Error = float(data_cpp[i][1])-float(data_java[i][1])
    if maxError < Error:
        maxError = Error
print("dio diff: " + str(maxError))


# StoneMask ----------------------------------------------------
with open("build/f0Stone_java.csv", 'r') as f:
    reader = csv.reader(f)
    data_cpp = [row for row in reader]
with open("debug/f0Stone_cpp.csv", 'r') as f:
    reader = csv.reader(f)
    data_java = [row for row in reader]

Error = 0.0
maxError = 0.0
for i in range(1, len(data_cpp)):
    Error = float(data_cpp[i][1])-float(data_java[i][1])
    if maxError < Error:
        maxError = Error
print("Stone diff: " + str(maxError))


# CheapTrick ----------------------------------------------------
# CSVファイルからデータを読み込む
data1 = pd.read_csv('build/sp_java.csv', header=None).values
data2 = pd.read_csv('debug/sp_cpp.csv', header=None).values
# 要素ごとの差を計算
diff = np.abs(data1 - data2)
# 最大値を出力
print("sp diff: " + str(np.max(diff)))


# D4C ----------------------------------------------------------
# CSVファイルからデータを読み込む
data1 = pd.read_csv('build/ap_java.csv', header=None).values
data2 = pd.read_csv('debug/ap_cpp.csv', header=None).values
# 要素ごとの差を計算
diff = np.abs(data1 - data2)
# 最大値を出力
print("ap diff: " + str(np.max(diff)))
