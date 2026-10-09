#!/usr/bin/env python3

# Arguments:
# check_results.py <err_tol|compare> <base_filename> <filename> <errInf tolerance> <errL2 tolerance> <solInf tolerance>
# This assumes displacement variables and attenuation variables have the same tolerance

import sys
import numbers

def run_checks(checks):
    fail = False
    for elem in checks:
        pass_check = (checks[elem][0] > checks[elem][1])
        print("%(check)s    tolerance: %(tol)f   calculated: %(val)f    pass: %(pass)s" %
                {"check": elem, "tol": checks[elem][0], "val": checks[elem][1], "pass": pass_check})
        if not pass_check: fail = True
    return fail

def main():
    if len(sys.argv) != 7: exit(1)

    if sys.argv[1] == "err_tol":
        # base_filename (argument 2) is ignored
        fp = open(sys.argv[3])
        fail = False
        tols = [float(sys.argv[i+4]) for i in range(3)]
        contents = fp.readlines()
        for i in range(len(contents)):
            line = contents[i]
            if line.find("Atten") == 0 or line.find("Disp") == 0:
                num_line = contents[i+1]
                vals = [float(x) for x in num_line.split()]
                checks = {"errInf": [tols[0], vals[0]], "errL2": [tols[1], vals[1]], "solInf": [tols[2], 1.0-vals[2]]}
                fail = fail or run_checks(checks)
        fp.close()
    elif sys.argv[1] == "compare":
        fail = False
        errInfTol = float(sys.argv[4])
        errL2Tol = float(sys.argv[5])
        base_file = open(sys.argv[2])
        test_file = open(sys.argv[3])
        base_data = base_file.readlines()
        test_data = test_file.readlines()
        base_file.close()
        test_file.close()
        if len(base_data) != len(test_data):
            print("ERROR: "+sys.argv[2]+" and "+sys.argv[3]+" are different lengths.")
            fail = True
        for i in range(min(len(test_data),len(base_data))):
            # some lines in the header are text only and must be ignored before converting to float
            base_line_ent = base_data[i].split()
            try:
                first_ent = float(base_line_ent[0]);
            except:
                #print("compare: line #", i, "base_line_ent[0]=' ", base_line_ent[0], "' could not be converted to a number, skipping this line")
                pass
                continue # skip this line

            base_line = [float(x) for x in base_data[i].split()]
            test_line = [float(x) for x in test_data[i].split()]
            if len(base_line) < 2:
                #print("base_line=' ", base_line, " ' has less than 2 numbers. Skipping this line")
                continue
            
            # TwilightErr blocks are "errInf errL2 solInf"; LambErr/PointSourceErr lines
            # are "t errInf errL2 solInf". Drop the time column so the error norms,
            # not t, are what gets compared.
            if len(base_line) == 4:
                base_line, test_line = base_line[1:], test_line[1:]
            if len(test_line) != len(base_line):
                print("ERROR: Line %d has %d values, reference has %d" % (i, len(test_line), len(base_line)))
                fail = True
                continue
            for name, tol, b, t in zip(("errInf", "errL2"), (errInfTol, errL2Tol), base_line, test_line):
                denom = max(abs(t), abs(b))
                rel = abs(b - t) / denom if denom > 0 else 0.0
                if rel > tol:
                    print("ERROR: Line %d %s tolerance: %g actual: %g (ref %g, got %g)" % (i, name, tol, rel, b, t))
                    fail = True
    else:
        print("Unknown test type: "+sys.argv[1])
        fail = True

    if fail: exit(1)
    exit(0)

if __name__ == "__main__":
    main()

