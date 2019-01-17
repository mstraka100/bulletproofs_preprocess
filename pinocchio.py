#!/usr/bin/python3
# Original code by Peter Wuille

import re
import sys

OUTPUT = re.compile("output ([0-9]+)")
CMD = re.compile("([a-z0-9-]+) in ([0-9]+) <([0-9 ]+)> out ([0-9]+) <([0-9 ]+)>")
MUL_CONST = re.compile("const-mul-([0-9a-f]+)")
MUL_NEG_CONST = re.compile("const-mul-neg-([0-9a-f]+)")
SKIP = re.compile("total | input")

assnfile = open("SHA256_Sample_Run1_optimized.in")

inputs = {}
for line in assnfile:
    vals = line.split(" ")
    inputs[int(vals[0])] = int(vals[1][:-1], 16)

for x in inputs.items():
    print("v%i = #%i" % (x[0], x[1]))

for line in sys.stdin:
    if "total" in line:
        continue
    if "input" in line:
        continue
    if "assert" in line:
        continue
    line = line.strip()
    sp = OUTPUT.fullmatch(line)
    if sp:
        print("debug v%i" % int(sp.group(1)))
        continue
    sp = CMD.fullmatch(line)
    if not sp:
        raise Exception("Unknown line: %s" % line)
    cmd = sp.group(1)
    num_in = int(sp.group(2))
    ins = [int(x) for x in sp.group(3).split(" ")]
    assert(len(ins) == num_in)
    num_out = int(sp.group(4))
    outs = [int(x) for x in sp.group(5).split(" ")]
    assert(len(outs) == num_out)
    if cmd == 'split':
        assert(len(ins) == 1)
        print("%s := v%i" % (",".join("v%i" % i for i in outs), ins[0]))
        continue
    if cmd == 'pack':
        assert(len(outs) == 1)
        print("v%i =: %s" % (outs[0], ",".join("v%i" % i for i in ins)))
        continue
    if cmd == 'xor':
        assert(len(ins) == 2)
        assert(len(outs) == 1)
        print("v%i = v%i ^ v%i" % (outs[0], ins[0], ins[1]))
        continue
    if cmd == 'add':
        assert(len(ins) >= 2)
        assert(len(outs) == 1)
        string = "v%i = v%i" % (outs[0], ins[0])
        for i in range(len(ins)-1):
            string += " + v%i" % ins[i+1]
        print(string)
        continue
    if cmd == 'mul':
        assert(len(ins) == 2)
        assert(len(outs) == 1)
        print("v%i = v%i * v%i" % (outs[0], ins[0], ins[1]))
        continue
    if cmd == 'zerop':
        assert(len(ins) == 1)
        assert(len(outs) == 2)
        print("v%i =? v%i" % (outs[1], ins[0]))
        continue 
    else:
        sp = MUL_CONST.fullmatch(cmd)
        if sp:
            assert(len(ins) == 1)
            assert(len(outs) == 1)
            print("v%i = v%i * %i" % (outs[0], ins[0], int(sp.group(1), 16)))
            continue
        sp = MUL_NEG_CONST.fullmatch(cmd)
        if sp:
            assert(len(ins) == 1)
            assert(len(outs) == 1)
            print("v%i = v%i * (-%i)" % (outs[0], ins[0], int(sp.group(1), 16)))
            continue
    raise Exception("Unknown command: %s" % cmd)