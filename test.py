#!/usr/bin/env python

import os
import subprocess

NUM_THREADS=[1, 2, 4]
# 1MB, 2MB, 10MB, 100MB
#SIZE = [25000]
#SIZE = [25000, 50000, 250000, 250000]
SIZE_IT = [ (25000, 500), (50000, 100), (250000, 50), (2500000, 10) ]
FILL = [round(i * 0.1, 2) for i in range(0, 11)]

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
f = open("timing.txt", "w+")
f.write("")
f.close()
for i in range(len(SIZE_IT)):
  for num_t in NUM_THREADS:
    for fill in FILL:
      size, it = SIZE_IT[i]
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
      with open("timing.txt", "a") as f:
        f.writelines(str(dat).strip("[]") + "\n")

print(test_run)
