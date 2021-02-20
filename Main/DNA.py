#! /usr/bin/python3
# -*- coding: utf-8 -*-

__Version__ = "0.1"
__Author__ = "pzweuj"
__Date__ = "20210220"

import os
import sys
import argparse
import yaml


def main(runinfo):
    config = open(runinfo, "r")
    configDict = yaml.load(config)
    
