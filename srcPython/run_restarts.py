#!/usr/bin/env python3

import argparse
import datetime as dt
import json
import os

# ----------------------------------------------------------------------
# Function to parse input arguments
# ----------------------------------------------------------------------

def get_args_restart():

    parser = argparse.ArgumentParser(description =
                                     'Run series of restarts for Aether')
    parser.add_argument('-input',
                        help = 'aether json file for whole run',
                        default = 'aether.json.whole')

    parser.add_argument('-mpi',
                        help = 'mpirun command',
                        default = 'mpirun')
    
    parser.add_argument('-post',
                        help = 'post processing command',
                        default = '~/bin/postAether.py')
    
    parser.add_argument('-aether',
                        help = 'aether command',
                        default = './aether')
    
    parser.add_argument('-rundir',
                        help = 'path to the default run directory',
                        default = './share/run')
    
    parser.add_argument('-restarts',
                        help = 'number of restarts',
                        default = 1, type = int)
    
    parser.add_argument('-ensembles',
                        help = 'number of ensemble members',
                        default = 3, type = int)
    
    parser.add_argument('-blocks',
                        help = 'number of blocks in each ensemble member',
                        default = 1, type = int)

    parser.add_argument('-totaltime',
                        help = 'total time (in seconds) of run (-1 = ignore)',
                        default = -1, type = int)

    parser.add_argument('-test', \
                        help='just output the files and dont run', \
                        action="store_true")
    
    parser.add_argument('-dowhole', \
                        help='include the whole run for comparison', \
                        action="store_true")
    
    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------
# Get the start time, end time, and delta-t from the whole run dict:
# ----------------------------------------------------------------------

def get_times(wholeDict):

    if ('StartTime' in wholeDict):
        st = wholeDict['StartTime']
        startTime = dt.datetime(st[0], st[1], st[2], st[3], st[4], st[5])
    else:
        startTime = None
    if ('EndTime' in wholeDict):
        et = wholeDict['EndTime']
        endTime = dt.datetime(et[0], et[1], et[2], et[3], et[4], et[5])
    else:
        endTime = None
    if (startTime and endTime):
        deltaTime = (endTime - startTime).total_seconds()
    else:
        deltaTime = -1.0
        
    return startTime, endTime, deltaTime

# ----------------------------------------------------------------------
# run os command
# ----------------------------------------------------------------------

def run(command, isTest, isVerbose):
    if (isVerbose):
        print(' -> Running : ', command)
    if (not isTest):
        os.system(command)
    return

# ----------------------------------------------------------------------
# Main code:
# ----------------------------------------------------------------------

if __name__ == '__main__':

    args = get_args_restart()

    isVerbose = True
    
    inFile = args.input
    nRestarts = args.restarts
    nMembers = args.ensembles
    nBlocks = args.blocks

    nProcs = nBlocks * nMembers

    aetherCommand = args.aether

    runDir = args.rundir
    notInRunDir = False
    
    if (not os.path.exists(aetherCommand)):
        print('-> Can not find aether command : ', aetherCommand)
        notInRunDir = True
        
    if (not os.path.exists('./UA')):
        print('-> Can not find UA directory!')
        notInRunDir = True

    if (notInRunDir):
        if (os.path.exists(runDir)):
            print('Found rundir... copying')
            dirPre = './run.restarts'
            cdPre = 'cd ' + dirPre + ' ; '
            cdPost = ' ; cd ..'
            command = 'rm -rf ' + dirPre
            command = 'cp -r ' + runDir + ' ' + dirPre
            run(command, args.test, isVerbose)
        else:
            print('-> Can not find rundir to copy!')
            print('use -rundir to tell where to find the default run directory')
            exit()
    else:
        cdPre = ''
        cdPost = ''
        dirPre = '.'
            
    runCommand = cdPre + args.mpi + ' -np %d ' % nProcs
    runCommand = runCommand + aetherCommand + cdPost

    postCommand = args.post + ' -rm'

    print('Number of Ensemble Members : ', nMembers)
    print('Number of blocks in each member : ', nBlocks)
    print(' --> nProcessors : ', nProcs)
    
    with open(inFile, 'r') as filePointer:
        # Reading from json file
        wholeDict = json.load(filePointer)
        filePointer.close()

    sTime, eTime, dTime = get_times(wholeDict)

    if (args.totaltime > 0):
        dTime = args.totaltime
    
    if (dTime < 0.0):
        print('Some sort of error determining time interval for run')
        print('please check StartTime and EndTime in whole run file')
        exit()

    if ('Perturb' in wholeDict):
        perturb = wholeDict['Perturb']
    else:
        perturb = {'f107': {'Mean' : 1.0,
                            'Std' : 0.0,
                            'Add': False,
                            'Constant': True}}
        wholeDict['Perturb'] = perturb

    if ('Ensembles' in wholeDict):
        wholeDict['Ensembles']['nMembers'] = nMembers
    else:
        wholeDict['Ensembles']['nMembers'] = nMembers
        
    interval = dTime / (nRestarts + 1)

    if ('Restart' in wholeDict):
        wholeDict['Restart']['do'] = False
        wholeDict['Restart']['dt'] = interval
    else:
        wholeDict['Restart'] = {'do' : False,
                                'dt' : interval}
        
    subDict = {'Ensembles' : {'nMembers': nMembers},
               'Restart': {'do' : True, 'dt': interval},
               'Perturb': perturb}

    print(' --> StartTime set to : ', sTime)
    
    for iRun in range(nRestarts + 2):

        isTest = args.test
        if (iRun == 0):
            if (not args.dowhole):
                isTest = True

        cRun = '%04d' % iRun
        print('Iteration: ', iRun)
        if (iRun > 0):
            endtime = sTime + dt.timedelta(seconds = interval * iRun)
        else:
            endtime = sTime + dt.timedelta(seconds = dTime)
        endtimeAsArray = [endtime.year,
                          endtime.month,
                          endtime.day,
                          endtime.hour,
                          endtime.minute,
                          endtime.second]
        print(' --> EndTime set to : ', endtime)

        if (args.test):
            outFile = dirPre + '/aether_' + cRun + '.json'
        else:
            outFile = dirPre + '/aether.json'

        if (iRun <= 1):
            wholeDict['EndTime'] = endtimeAsArray
            jsonObject = json.dumps(wholeDict)
        else:
            subDict['EndTime'] = endtimeAsArray
            jsonObject = json.dumps(subDict)

        print(' --> Writing File : ', outFile)
        with open(outFile, "w") as outfile:
            outfile.write(jsonObject)        
    
        # ------------------------------------------------------------------
        # run Aether
        command = runCommand + ' 2>&1 | tee output_' + cRun + '.log'
        run(command, isTest, isVerbose)

        # ----------------------
        # remove target directories
        command = \
            'cd ' + dirPre + '/UA ; ' + \
            'rm -rf output_' + cRun + \
            ' ; cd ..' + cdPost
        run(command, isTest, isVerbose)

        command = \
            'cd ' + dirPre + '/UA ; ' + \
            'rm -rf restartOut_' + cRun + \
            ' ; cd ..' + cdPost
        run(command, isTest, isVerbose)
            
        # ----------------------
        # Post Processing
        command = \
            'cd ' + dirPre + '/UA/output ; ' + \
            postCommand + \
            ' ; cd .. ; ' + \
            'mv output output_' + cRun + \
            ' ; mkdir output ' + \
            ' ; cd ..' + cdPost
        run(command, isTest, isVerbose)

        # ----------------------
        # Set up restart directories
        command = \
            'cd ' + dirPre + '/UA ; ' + \
            'mv restartOut restartOut_' + cRun + \
            ' ; rm -f ./restartIn ' + \
            ' ; ln -s restartOut_' + cRun + ' ./restartIn ' + \
            ' ; mkdir restartOut ' + \
            ' ; cd ..' + cdPost
        run(command, isTest, isVerbose)
        
