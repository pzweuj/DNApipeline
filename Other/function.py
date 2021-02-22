#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.1"
__Author__ = "pzweuj"
__Date__ = "20210222"


import os
import time
import sys
import yaml

# 配置文件读入
def getRunningInfo(runInfo):
    with open(runInfo, "r") as stream:
        try:
            configDict = yaml.safe_load(stream)
            return configDict
        except yaml.YAMLError as exc:
            print(exc)

# 获得当前时间
def getTime():
    localTime = time.asctime(time.localtime(time.time()))
    return localTime

# 获得脚本路径
def getAbsPath():
    now = os.path.abspath(os.path.dirname(sys.argv[0]))
    return now

# 新建文件夹
def mkdir(folder_path):
    folder = os.path.exists(folder_path)
    if not folder:
        os.makedirs(folder_path)
        print("Create new folder " + folder_path + " done!")
    else:
        print("Folder " + folder_path + " already exists!")

# 获得当前路径
def getCurrentPath():
    return os.getcwd()