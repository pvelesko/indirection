#!/usr/bin/env python

import os
import subprocess

NUM_THREADS=[1, 2, 4]
# 1MB, 2MB, 10MB, 100MB
SIZE = [25000]
#SIZE = [25000, 50000, 250000, 250000]
FILL = [i * 0.1 for i in range(0, 10)]

def parse_test(s):
  test = s.split("\n")
  for i in range(len(test)):
    if "Test0" in test[i]:
      t0 = float(test[i].split(" ")[0])
    if "Test1" in test[i]:
      t1 = float(test[i].split(" ")[0])
    if "Test2" in test[i]:
      t2 = float(test[i].split(" ")[0])
    if "Test3" in test[i]:
      t3 = float(test[i].split(" ")[0])
    if "Test4" in test[i]:
      t4 = float(test[i].split(" ")[0])
  return [t1, t2, t3, t4]

e = dict(os.environ)

test_run = []
for size in SIZE:
  for num_t in NUM_THREADS:
    for fill in FILL:
      it = 10
      cmd = ["./intel", str(size), str(it), str(fill)]
      e['OMP_NUM_THREADS'] = str(num_t)
      p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=e, universal_newlines=True)
      out, err = p.communicate()
      print((out))
      print((err))
      res = parse_test(out)
      dat = [num_t, size, it, fill] + res
      test_run.append(dat)
      print(dat)

print(test_run)
