#! /usr/bin/env python

from IPython import embed
import argparse
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import re


def process_log(args):
    log = open(args.infile)
    # find start of transient simulation
    for line in log:
        if (str.find(line, 'qucsator') >= 0):
            break

    # collect data
    cpu = list()
    mem = list()
    times = list()
    icpu = list()
    imem = list()
    itimes = list()
    i = 0
    nit = [0]
    t = 0
    ld = list()
    f = list()
    I = list()
    h = list()
    tvec = list()
    itvec = list()
    cput = list()
    res = None
    for line in log:
        if str.find(line, 'CPU = ') >= 0:
            res = re.search(r"\([A-Za-z]+ (?P<date>.*), CPU = (?P<cpu>.*)s, Mem = (?P<mem>.*)Mb", line)
            # date format Wed Jun  3 01:15:42 2015
            itimes.append(datetime.strptime(res.group('date'),
                          '%b %d %H:%M:%S %Y'))
            icpu.append(float(res.group('cpu')))
            imem.append(float(res.group('mem')))
        if str.find(line, 'calcTR: ') >= 0 and str.find(line, '@') == -1:
            res = re.search(r"calcTR: (?P<tsim>.*); (?P<tit>.*); (?P<i>.*); "
                            r"(?P<t>.*); (?P<ld>.*); (?P<f>.*); (?P<I>.*)",
                            line)
            if (float(res.group('t'))-t > 0):
                i = i+1
                nit.append(0)
                h.append(float(res.group('t')) - t)
                t = float(res.group('t'))
            cput.append(float(res.group('tsim')))
            nit[i] = nit[i] + 1
            tvec.append(t)
            itvec.append(i)
            I.append(float(res.group('I')))
            f.append(float(res.group('f')))
            ld.append(float(res.group('ld')))
            times.append(itimes)
            cpu.append(icpu)
            mem.append(imem)
            itimes = list()
            icpu = list()
            imem = list()
#    realgetdptimes = [(t[-1]-t[0]).total_seconds() for t in times]
#    cpugetdptimes = [t[-1] for t in cpu]
#    realgetdpsolve = [(t[-4]-t[1]).total_seconds() for t in times]
#    cpugetdpsolve = [t[-4]-t[1] for t in cpu]
    if args.p:
        fig = plt.figure()
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
        ax1.plot(np.array(nit),
                 label='Newton Iterations/timestep')
        ax1.plot(np.array(itvec), -np.array(I),
                 label='Current')
        ax2.plot(np.array(h),
                 label='Stepsize')
        ax1.legend()
        ax2.legend()
        fig.show()
    if args.csv:
        np.savetxt(args.pre+'h.csv',
                   np.array([range(0, len(h)), h]).transpose(),
                   delimiter=',')
        np.savetxt(args.pre+'it.csv',
                   np.array([range(0, len(nit)), nit]).transpose(),
                   delimiter=',')
        np.savetxt(args.pre+'i.csv',
                   np.array([np.array(itvec), -np.array(I)]).transpose(),
                   delimiter=',')
    if args.i:
        embed()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Post Process ')
    parser.add_argument('infile', type=str, help='File to read log from')
    parser.add_argument('--pre', type=str, help='Prefix for output files')
    parser.add_argument('--csv', action='store_const', const=True,
                        default=False, help='write to csv files')
    parser.add_argument('-p', action='store_const', const=True,
                        default=False, help='enable plotting of results')
    parser.add_argument('-i', action='store_const', const=True,
                        default=False, help='enable interactive console')
    args = parser.parse_args()
    if args.infile:
        process_log(args)
