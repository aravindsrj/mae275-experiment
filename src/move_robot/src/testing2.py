#! /usr/bin/env python

import os

p = [1,2,3,4]
os.system('python testing.py ' + str(p)[1:-1].replace(',',' '))