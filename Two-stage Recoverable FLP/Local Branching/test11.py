# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 10:02:15 2019

@author: User
"""

# colours = [[[3],[2]],[[2],[3]],[[1],[3]]]
# c=[1]
#
# aa = [(i, colour.index(c))
#  for i, colour in enumerate(colours)
#  if c in colour]
#
# b = next(((i, colour.index(c))
#       for i, colour in enumerate(colours)
#       if c in colour),
#      None)

colours = [[[0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0], [330.93640888316133, 3]], [[1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [310.3468033268924, 1]]]
c= [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
aa = [(i, colour.index(c))
 for i, colour in enumerate(colours)
 if c in colour]
b = next(((i, colour.index(c))
      for i, colour in enumerate(colours)
      if c in colour),
     None)
