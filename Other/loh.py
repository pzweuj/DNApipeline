#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.01"
__Author__ = "pzweuj"
__Date__ = "20210402"

import os
import sys

FUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(FUN_DIR)
from Other.function import mkdir


class LOH(object):
	pass