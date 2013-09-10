#!/usr/bin/env python
# coding: utf-8

file = open("detector_config.txt", "w")

for z in range(1, 61):
	group = 0
	for y in range(1, 11):
		for x in range(1, 11):
			a = x * 100.
			b = y * 100.
			c = z * 10.
			group += 1.
			if (((z >= 25) and (z <= 35)) or (z == 43) or (z == 49) or (z == 54)):
				perm = 0.2
			else:
				perm = 1.0
			file.write("{:.0f}\t{:.0f}\t{:.0f}\t{:.1f}\t{:.0f}".format(a, b, c, perm, group))
			file.write("\n")
file.close()
