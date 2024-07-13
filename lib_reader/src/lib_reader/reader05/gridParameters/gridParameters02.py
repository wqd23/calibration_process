# -*- coding:utf-8 -*- 

"""
GRID02 Internal Parameters
----------

Auxiliary module to store internal parameters for GRID02 and GRID02-based datapack\n
`v1.0.0` for GRID02 calibration result analysis, also quick view of GRID02 in-orbit data, by ghz and ydx\n
"""

#****************************************************************************************************************************************************************
#********************************************************************Internal parameters********************************************************************
#****************************************************************************************************************************************************************

#--------------------------------------------------------------------------Config files--------------------------------------------------------------------------
#Config file for file import
importConfig = './config_files/02/import_config.json'

#EC config file for x-ray data readout
xrayConfig = './config_files/02/ec_xray_config.csv'

#EC config file for source data readout
sourceConfig = './config_files/02/ec_src_config.csv'

#EC config file for HPGe data readout
HPGeConfig = './config_files/02/hpge_config.csv'

#EC config file for x-ray data processing
xrayProcessConfig = './config_files/02/ec_xray_config.json'

#EC config file for source data processing
sourceProcessConfig = './config_files/02/ec_src_config.json'

#Import config file for temp-bias processing
importTempbiasConfig = './config_files/02/tempbias_import_config.json'

#Import config file for EC processing
importECConfig = './config_files/02/ec_import_config.json'

#------------------------------------------------------------------General readout options------------------------------------------------------------------
#Deduplicate the data instead of readout, internal flag
#Use this for LARGE files, for the sake of RAM usage
#DO NOT ADD ANY OPTIONS OTHER THAN '-hn' TO AVOID ERROR
deduplicateFile = False

#In-pack uscount reset point finding, for special cases. Usually not needed, but in cases of usocunt overflow, switch this flag to True
inPackReset = False

#Reference of files to ignore CRC check
ignoreCrcRef = ['191229212255_Na22_1h_20mv_COM6-Data.txt', '191230145530_Co60_2h_20mv_COM6-Data.txt', \
    '191229193211_Am241_5m_20mv_COM6-Data.txt', '191229141203_Am241_5m_COM6-Data.txt', '191229172209_Am241_5m_20mv_COM6-Data.txt', \
    '191231120244_0_10m_20mv_COM6-Data.txt']

#Pattern reference for raw data readout
#bin: binary files; hex: heprint text files; -New: new programme(6th ver.) format
pattrenRef = {
        'binSci':                     rb'\xAA\xBB\xCC.{' + bytes('{}'.format(504), encoding = 'utf-8') + rb'}\xDD\xEE\xFF', 
        'binTel':                     rb'\x12\x34\x56.{' + bytes('{}'.format(490), encoding = 'utf-8') + rb'}\x78\x9A\xBC', 
        'hexSciNew':              r'AA BB CC ([0-9A-F]{2} ){504}DD EE FF', 
        'hexTelNew':              r'12 34 56 ([0-9A-F]{2} ){490}78 9A BC', 
        'hexSci':                    r'AA BB CC ([0-9A-F]{2} ){496}DD EE FF', 
        'hexTel':                    r'01 23 45 ([0-9A-F]{2} ){490}67 89 01', 
}

#Frequency of internal crystal oscillator, in Hz
#For data before 20200620, internal frequency was 24.05MHz, and 20MHz afterwards
internalFreq = 20.0e6

#Resistance of internal voltage and current monitor circuit, in kOhm
internalResistance = 1.1

#Launch time in the format 'YYYYmmdd', for in-orbit data qiuck-view
launchTime = '20201106'

#-----------------------------------------------------------------Advanced readout options-----------------------------------------------------------------
#---------------------------------------------------------------------------Time cut----------------------------------------------------------------------------
#Reference of time cut for some specified files
cutFileRef = {
        'qh-Am241-30度-0.dat':         650, 
        'qh-Am241-35度-2.dat':         220, 
        'qh-Am241-40度-1.dat':         110, 
        'Am241-29degree-0.dat':       2960, 
        'Am241-29degree-1.dat':       60, 
        'Am241-35degree-0.dat':       1220, 
        'Am241-35degree-1.dat':       70, 
        'Am241-40degree-0.dat':       400, 
        'Am241-40degree-1.dat':       80, 
        'Am241-40degree-2.dat':       90, 
        'Am241-0度-1.dat':               180, 
        'Am241-10度-1.dat':             1710, 
        'Am241-20度-1.dat':             1430, 
        'Am241-25度-1.dat':             1350, 
        'Am241-30度-1.dat':             1300, 
        'Am241CI-0度-0.dat':            1100, 
        'Am241CI-10度-0.dat':          660, 
        'Am241CI-20度-0.dat':          1480, 
        'Am241CI-25度-0.dat':          1250, 
        'Am241CI-30度-0.dat':          1300, 
        'Am241CI-35度-0.dat':          1120, 
    }

#-----------------------------------------------------------------Advanced readout options-----------------------------------------------------------------
#------------------------------------------------------------------Multiple runs (.dat) files------------------------------------------------------------------
#Reference of required scan for some specified files
reqFileRef = {
        'qh-Am241-30度-0.dat':       [2, 2], 
        'qh-Am241-35度-2.dat':       [2, 1], 
        'qh-Am241-40度-1.dat':       [2, 1], 
        'Am241-40degree-0.dat':       [1, 0], 
        'grid-unpack.dat':                  [7, 1], 
    }

#Reference of multiple bias runs .dat files, for multiple temp-bias runs files ONLY, in volts
biasRef = {
        'bias_set':         [27.0, 27.5, 28.0, 28.3, 28.5, 28.7, 29.0],
        'bound':            0.05,
    }

#Reference for temperature runs, for multiple temp-bias runs files ONLY
tempRef = {
        'set':          [0, 10, 20, 25, 30],
        'bound':        1,
    }

#Reference for unused runs, for multiple temp-bias runs files ONLY
unusedRef = {
        'DUMMY':         [],
        # 'Am241CI-20度-0.dat':        [27.5, 28.0, 28.3, 28.5, 28.7, 28.9],
        # 'Am241CI-30度-0.dat':        [],
    }

#Reference for unused science runs, for multiple runs files ONLY
unusedSciRef = {
    'DUMMY':        [], 
    '20201207_1607_unpack-20201207-1607281201.dat':     [0], 
    '20201208_1408_unpack-20201208-1607407201.dat':     [0, 5, 20, 23, 24], 
    '20210121_2219_unpack-20210121-1611238339.dat':     [0], 
    '20210125_1033_unpack-202101231033.dat':                [14, 15], 
}

#Reference for unused telemetry runs, for multiple runs files ONLY
unusedTelRef = {
    'DUMMY':        [], 
    '20201206_1424_unpack.dat':         [31], 
    '20201207_1607_unpack-20201207-1607281201.dat':     [1, 27], 
    '20201208_0000_unpack-20201207-1607356644.dat':     [1, 27], 
    '20201208_1408_unpack-20201208-1607407201.dat':     [4, 19, 22, 23], 
}

#---------------------------------------------------------------Live time fit and plot options---------------------------------------------------------------
#Cut-off value for live time fit, in seconds
liveTimeCut = 0.1

#Internal MCU dead time
deadTimeMCU = 50e-6
