# -*- coding:utf-8 -*- 

"""
GRID03 Internal Parameters
----------

Auxiliary module to store internal parameters for GRID03 and GRID03-based datapack\n
`v1.0.0` for GRID03 calibration result analysis, also quick view of GRID03 in-orbit data, by ghz and cjr\n
"""

#****************************************************************************************************************************************************************
#********************************************************************Internal parameters********************************************************************
#****************************************************************************************************************************************************************

#--------------------------------------------------------------------------Config files--------------------------------------------------------------------------
#Config file for file import
importConfig = './config_files/03/import_config.json'

#EC config file for x-ray data readout
xrayConfig = './config_files/03/ec_xray_config.csv'
xray_json_config = './config_files/03/ec_xray_json_config.json'
#EC config file for source data readout
sourceConfig = './config_files/03/ec_src_config.csv'
source_json_config='./config_files/03/ec_src_json_config.json'
#EC config file for HPGe data readout
HPGeConfig = './config_files/03/hpge_config.csv'

#EC config file for x-ray data processing
xrayProcessConfig = './config_files/03/ec_xray_config.json'

#EC config file for source data processing
sourceProcessConfig = './config_files/03/ec_src_config.json'

#Import config file for temp-bias processing
importTempbiasConfig = './config_files/03/tempbias_import_config.json'

#Import config file for EC processing
importECConfig = './config_files/03/ec_import_config.json'

#------------------------------------------------------------------General readout options------------------------------------------------------------------
#Reference of files to ignore CRC check
ignoreCrcRef = []

#Number of events contained in feature mode data packs
featureEventNum = 20

#Length of UDP packages in rundata, in bytes
udpPackLen = 16384

#Length of HK data packages, in bytes
#hkDataLen = 52
hkDataLen = 95

#Length of IV scan data, in bytes
IVPackLen = 600

#Length of Vbr scan data, in bytes
VbrPackLen = 600

#Length of timeline data packages, in bytes
timelineDataLen = 16

#Pattern reference for raw data readout
pattrenRef = {
        'udp':                         rb'\x1a\xcf\xfc\x1d.{2}\x11\x19.{2}\x88.{' + bytes('{}'.format(udpPackLen), encoding = 'utf-8') + rb'}.{1}\x2e\xe9\xc8\xfd', 
        'HK':                          rb'\x1a\xcf\xfc\x1d.{2}\x11\x19.{2}\x89.{1}.{' + bytes('{}'.format(hkDataLen), encoding = 'utf-8') + rb'}.{1}\x2e\xe9\xc8\xfd',
        'HK_new':                      rb'\x1a\xcf\xfc\x1d\x00\x00\x00\x00.{4}\x00\x00\x00\x01',
        'HK_03b': b'\\x1a\\xcf\\xfc\\x1d.{2}\\x11\\x19.{2}\\x89.{1}.{52}.{1}\\x2e\\xe9\\xc8\\xfd',
        'timeline':                  rb'\x1a\xcf\xfc\x1d.{2}\x11\x19.{2}\x90.{' + bytes('{}'.format(timelineDataLen), encoding = 'utf-8') + rb'}.{1}\x2e\xe9\xc8\xfd', 
        'time_new':                    rb'\x1a\xcf\xfc\x1d\x00\x00\x00\x00.{4}\x00\x00\x00\x07', 
        'sciNewHead':            rb'\x5a\x5a\x5a\x5a\x5a\x5a\x5a\x5a.{', 
        'sciNewTail':               rb'}\xaa\xaa\xaa\xaa\xaa\xaa\xaa\xaa', 
        'sciFeatureNew':         rb'\x5a\x5a\x99\x66\x99\x66\x5a\x5a.{' + bytes('{}'.format(featureEventNum * 24 + 8), encoding = 'utf-8') + \
                rb'}\xaa\xaa\x99\x66\x99\x66\xaa\xaa', 
        'sciHead':                   rb'\x5a\x5a\x5a\x5a\x5a\x5a\x5a\x5a.{', 
        'sciTail':                     rb'}\x51\x52\x53\x54\x5a\x5b\x5c\x5d', 
        'IV':                           rb'\x1a\xcf\xfc\x1d.{2}\x11\x19.{2}\xAA.{' + bytes('{}'.format(IVPackLen), encoding = 'utf-8') + rb'}.{1}\x2e\xe9\xc8\xfd', 
        'Vbr':                         rb'\x1a\xcf\xfc\x1d.{2}\x11\x19.{2}\xAB.{' + bytes('{}'.format(VbrPackLen), encoding = 'utf-8') + rb'}.{1}\x2e\xe9\xc8\xfd', 
}

#Frequency of internal crystal oscillator, in Hz
internalFreq = 100.e6

#Number of points in the head of pulse array to use for baseline fit
baseLen = 1

#Use measured baseline value to subtract baseline from data
useFitBaseline = True

#Use baseline value in science data pack for each event to subtract baseline from data
useEventBaseline = True

#Number of UDP packages to process in single run when in multiple readout mode, or -1 for single readout mode
maxUdpReadout = 10000

#Launch time in the format 'YYYYmmdd', for in-orbit data qiuck-view
launchTime = '20201222'

#Check pulse shape for pulse mode data
checkPulseShape = False

#Plot figures in data readout
plotFigReadout = True

#Threshold for splitting science and telemetry runs, in seconds
splitRunTime = 30

#-----------------------------------------------------------------Advanced readout options-----------------------------------------------------------------
#---------------------------------------------------------------------------Time cut----------------------------------------------------------------------------
#Reference of time cut for some specified files
cutFileRef = {
        '13.4_ch2_100s_rundata2020-09-05-20-24-52.dat':         12020, 
        '15.6_ch0_60s_rundata2020-09-05-20-12-16.dat':           11260, 
        '20.9_ch0_30s_rundata2020-09-05-19-50-11.dat':           9900, 
        '24.2_ch2_30s_rundata2020-09-05-19-40-52.dat':           9372, 
        '34.5_ch3_30s_rundata2020-09-05-19-13-14.dat':           7724, 
        '38.3_ch1_40s_rundata2020-09-05-19-03-25.dat':           7120, 
        '42.5_ch1_60s_rundata2020-09-05-18-49-11.dat':           6280, 
        '49.2_ch3_30s_rundata2020-09-16-13-49-00.dat':           4541, 
        '51.4_ch0_30s_rundata2020-09-16-14-18-21.dat':           6300, 
        '92.0_ch0_60s_rundata2020-09-17-13-46-38.dat':           12200, 
        '92.0_ch2_60s_rundata2020-09-17-13-48-41.dat':           12674, 
        '92.0_ch3_60s_rundata2020-09-17-13-51-08.dat':           12571, 
        '101.8_ch2_60s_rundata2020-09-17-14-05-56.dat':          13595, 
        '113.7_ch0_60s_rundata2020-09-17-14-21-03.dat':          14400, 
        '10C_28.0V_4m_rundata2020-09-11-18-46-06.dat':          2500, 
        '15C_29.0V_4m_rundata2020-09-11-21-41-17.dat':          6940, 
        '20C_26.5V_4m_rundata2020-09-11-23-22-50.dat':          13000, 
        '25C_26.5V_4m_rundata2020-09-12-02-06-13.dat':          22800, 
        '25C_28.3V_4m_rundata2020-09-12-02-46-06.dat':          25200, 
        '30C_27.5V_4m_rundata2020-09-12-05-08-21.dat':          5050, 
        '35C_27.5V_4m_rundata2020-09-12-08-09-29.dat':          3030, 
        '40C_28.5V_4m_rundata2020-09-12-10-55-33.dat':          3570, 
        #--------------------------------------------------GRID03b files--------------------------------------------------
        'jly_18p0_ch0_30s_rundata2021-04-29-15-21-50.dat':          5025, 
        'jly_24p2_ch0_30s_rundata2021-04-29-15-10-52.dat':          4380, 
        'jly_31p2_ch3_30s_rundata2021-04-29-14-57-31.dat':          3580, 
        'jly_38p2_ch0_30s_rundata2021-04-29-14-38-47.dat':          2460, 
        'jly_45p9_ch0_30s_rundata2021-04-29-17-41-58.dat':          6050, 
        'jly_45p9_ch3_30s_rundata2021-04-29-17-38-24.dat':          5840, 
        'jly_77p0_ch3_30s_rundata2021-04-29-19-00-01.dat':          10740, 
    }

#-----------------------------------------------------------------Advanced readout options-----------------------------------------------------------------
#------------------------------------------------------------------Multiple runs (.dat) files------------------------------------------------------------------
#Reference of required scan for some specified files
reqFileRef = {
        'grid-unpack.dat':                  [7, 1], 
    }

#Reference of multiple bias runs .dat files, for multiple temp-bias runs files ONLY, in volts
biasRef = {
        'bias_set':         [26.5, 27.0, 27.5, 28.0, 28.3, 28.5, 28.7, 29.0],
        'bound':            0.05, 
    }

#Reference for temperature runs, for multiple temp-bias runs files ONLY
tempRef = {
        'set':          [10, 15, 20, 25, 30, 35, 40], 
        'bound':        1, 
    }

#Reference for unused runs, for multiple temp-bias runs files ONLY
unusedRef = {
        'DUMMY':         [],
    }

#Reference for unused science runs, for multiple runs files ONLY
unusedSciRef = {
    'DUMMY':        [], 
}

#Reference for unused telemetry runs, for multiple runs files ONLY
unusedTelRef = {
    'DUMMY':        [], 
}

#--------------------------------------------------------------------General plot options--------------------------------------------------------------------
#Plot and fit range for baseline fits
baselineRange = [500, 2000]

#Use specified bin edges to generate spectrum for baseline fit (NOT RECOMMENDED)
useBaselineEdges = False

#Bin width for baseline fits
baselineBinWidth = 5

#Pulse number to check
checkPulseNum = 0
