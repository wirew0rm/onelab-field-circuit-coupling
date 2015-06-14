#! /usr/bin/env python
# coding=utf-8

# Wrapper to invoke qucsator with the onelab socket
# set in the modified netlist file.
#
# Alexander Krimm <alex@stud.tu-darmstadt.de>
# TU-Darmstadt, 2015

import subprocess
import onelab
import os

c = onelab.client(__file__)

# Initialize Onelab Parameter Variables
# qucs solver Parameters
start = c.defineNumber('Circuit/Start Time [ms]', value=0)
stop = c.defineNumber('Circuit/Stop Time [ms]', value=0.6)
steps = c.defineNumber('Circuit/Number of Steps', value=50)
integrationMethod = c.defineString('Circuit/Integration Method',
                                   value='Trapezoidal',
                                   choices={'Euler',
                                            'Trapezoidal'})
order = c.defineNumber('Circuit/Integration Order', value=2,
                       choices={1, 2, 3, 4, 5, 6})
hinit = c.defineNumber('Circuit/Initial Step [ns]', value=1)
MinStep = c.defineNumber('Circuit/Minimum Step Size [s]', value=10e-16)
MaxIter = c.defineNumber('Circuit/Maximum Number of Iterations', value=20)
initialDC = c.defineString('Circuit/initialize with DC solution', value='no',
                           choices={'yes', 'no'})
# File names
circuitfile = c.defineString('Circuit/Circuit Name', value='rl',
                             choices={'rlc', 'rl', 'psu'})
getdpcheckcmd= c.defineString('Circuit/Command to check getdp',
                value='getdp inductor_simple/inductor.pro -pre -msh inductor_simple/inductor.msh',
                choices={'getdp inductor_simple/inductor.pro -pre -msh inductor_simple/inductor.msh',
                         'getdp inductor/inductor.pro -pre -msh inductor/inductor.msh',
                         'getdp transformer/transformer.pro -pre -msh transformer/transformer.msh'})

# Parameters for onelabgetdp.so module
verbose = c.defineString('Circuit/verbose', value='yes',
                         choices={'yes', 'no'})
gedtpcmd = c.defineString('Circuit/Command to launch GetDP',
                          value='getdp inductor_simple/inductor.pro -pre '
                          '-cal -pos -msh inductor_simple/inductor.msh '
                          '-gmshread inductor_simple/res/a.pos',
                          choices={'./licurve.py',
                                   'getdp inductor/inductor.pro -pre '
                                   '-cal -pos -msh inductor/inductor.msh '
                                   '-gmshread inductor/res/a.pos',
                                   'getdp transformer/transformer.pro -pre '
                                   '-cal -pos -msh transformer/transformer.msh '
                                   '-gmshread transformer/res/a.pos',
                                   'getdp inductor/inductor_simple.pro -pre '
                                   '-cal -pos -msh inductor_simple/inductor.msh '
                                   '-gmshread inductor_simple/res/a.pos', })
Iparam = c.defineString('Circuit/Parameter to export current to',
                        value='Input/4Coil Parameters/0Current (rms) [A]',
                        choices={'LICurve/I',
                                 'Input/4Coil Parameters/0Current (rms) [A]'})
FluxParam = c.defineString('Circuit/Parameter for Flux',
                           value='Output/40Flux [Wb]',
                           choices={'LICurve/Flux',
                                    'Output/40Flux [Wb]',
                                    'Output/54Inductance from Flux[H]'})
LdParam = c.defineString('Circuit/Parameter for differential Inductance',
                         value='Output/55Inductance from Flux Diff[H]',
                         choices={'Output/55Inductance from Flux Diff[H]',
                                  'diffquot',
                                  'LICurve/Ld',
                                  'Lchord'})
# let getdp announce its parameters
if c.action == 'check':
    c.runSubClient('getdp', getdpcheckcmd)

if c.action == 'compute':
    infile = os.getcwd() + '/' + circuitfile + '/' + circuitfile + '.sch'
    outfile = os.getcwd() + '/' + circuitfile + '/' + circuitfile + '.dat'

    # define output files, that can be archived along with the
    # onelab database after each run
    c.outputFiles([circuitfile + '/' + circuitfile + '.dat',
                   circuitfile + '/' + circuitfile + '.dpl'])

    # convert schematic file to netlist
    subprocess.call(["qucs", "-n", "-i", infile, "-o", infile + ".netlist"])

    # manipulate netlist
    ol_device = False
    solver = False
    with open(infile + ".netlist", 'r') as inf, \
            open(infile + ".netlist.tmp", 'w') as outf:
        for line in inf:
            cont = line.split()
            # Modify lines containing Inductors named "onelab"
            if (len(cont) and cont[0] == "L:onelab"):
                c.sendInfo("substituting inductance with Field device")
                args = ["OL:onelab"]    # change component id to qucsonelab
                args.append(cont[1])    # Node 1 keep from original
                args.append(cont[2])    # Node 2 keep from original
                # args.append("I=\"1 A\"")# Set initial Current
                args.append("\n")
                line = ' '.join(args)
                ol_device = True
            # Modify lines containing Mutal inductances named "onelab"
            if (len(cont) and cont[0] == "MUT:onelab"):
                c.sendInfo("substituting Mutal Inductance with Field device")
                args = ["OL4:onelab"]   # change component id to qucsonelab
                args.append(cont[1])   # Node 1 keep from original
                args.append(cont[2])   # Node 2 keep from original
                args.append(cont[3])   # Node 3 keep from original
                args.append(cont[4])   # Node 4 keep from original
                args.append("\n")
                line = ' '.join(args)
                ol_device = True
            # Modify lines containing transient solver named solve
            elif (len(cont) and cont[0] == '.TR:solve' and not solver):
                c.sendInfo("Setting circuit simulation parameters")
                args = ['.TR:solve', 'Type="lin"', 'reltol="0.001"',
                        'abstol="1 pA"', 'vntol="1 uV"', 'Temp="26.85"',
                        'LTEreltol="1e-3"', 'LTEabstol="1e-6"',
                        'LTEfactor="1"', 'Solver="CroutLU"',
                        'relaxTSR="no"', 'MaxStep="0"']
                args.append('Start="' + str(start) + ' ms"')
                args.append('Stop="' + str(stop) + ' ms"')
                args.append('Points="' + str(steps) + '"')
                args.append('IntegrationMethod="' + integrationMethod + '"')
                args.append('Order="' + str(order) + '"')
                args.append('InitialStep="' + str(hinit) + ' ns"')
                args.append('MinStep="' + str(MinStep) + '"')
                args.append('MaxIter="' + str(MaxIter) + '"')
                args.append('initialDC="' + initialDC + '"')
                args.append("\n")
                line = ' '.join(args)
                solver = True
            outf.write(line)

    if not solver:
        c.sendError("Netlist doesn't contain transient simulation named solve")
    elif not ol_device:
        c.sendError("Netlist doesn't contain inductivity named onelab")
    else:
        c.runSubClient('qucsator', "./qucs-wrapper.py -i " + infile +
                       " -o " + outfile)
    # remove temporary files
    os.remove(infile + ".netlist")
    os.remove(infile + ".netlist.tmp")
