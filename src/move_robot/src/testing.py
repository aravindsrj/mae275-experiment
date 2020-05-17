#! /usr/bin/env python

import rospy
import sys
import os
from geometry_msgs.msg import Point

print("No. of args = %d", len(sys.argv))
c = []
for i in range(1,len(sys.argv)):
    c.insert(i-1,Point(1,2,0))
    a = int(sys.argv[i])
    print a

# c[0].x = 5
# c[0].y = 10
print c[0].x

if __name__ =='__main__':
    print("name = main")