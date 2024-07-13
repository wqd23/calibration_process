# -*- coding:utf-8 -*- 

"""
GRID common Internal Parameters
----------

Auxiliary module to store common internal parameters for all datapack versions\n
`v1.0.0` for GRID02 calibration result analysis, also quick view of GRID02 in-orbit data, by ydx and ghz\n
"""

#****************************************************************************************************************************************************************
#********************************************************************Internal parameters********************************************************************
#****************************************************************************************************************************************************************

#----------------------------------------------------------------------------App info----------------------------------------------------------------------------
#Current GridDataProcessor version
currentVersion = '1.0.0'

#--------------------------------------------------------------------------Config files--------------------------------------------------------------------------
#Config file for file output
outputConfig = './config_files/02/output_config.json'

#Config file for file import
importConfig = './config_files/02/import_config.json'

#Temperature-bias correction reference file
tempBiasConfig = './config_files/02/tempbias_coef_202005.json'

#EC correction reference files
ecConfig = ['./config_files/02/ec_coef_ch0.json', 
                 './config_files/02/ec_coef_ch1.json', 
                 './config_files/02/ec_coef_ch2.json', 
                 './config_files/02/ec_coef_ch3.json', ]

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

#Import config file for EC processing
importECConfig = './config_files/02/ec_import_config.json'

#------------------------------------------------------------------General readout options------------------------------------------------------------------
#Deduplicate the data instead of readout, internal flag
#Use this for LARGE files, for the sake of RAM usage
#DO NOT ADD ANY OPTIONS OTHER THAN '-hn' TO AVOID ERROR
deduplicateFile = False

#Report step of deduplicating files, any non-positive integer for no output
deduplicateReportStep = 10000

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
internalFreq = 24.05e6

#Resistance of internal voltage and current monitor circuit, in kOhm
internalResistance = 1.1

#Check temperature data for error
tempCheck = False

#UTC check to remove data with UTC error
utcCheck = False

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

#--------------------------------------------------------------------General plot options--------------------------------------------------------------------
#Plot settings for saving figure directly or creating a new plot window
saveFigNoPlot = False

#Path to save figures
saveFigPath = './figure'

#Savefig settings, whether the data file uses its own savefig path(create new sub-directories for each new data file)
saveFigOwnPath = False

#Plot spectrum with log scale, to make some peaks more visible
plotSpecLogScale = False

#Starting position when plotting spectrum with log scale
logScaleBegin = 100.

#Display UTC time in all raw spectrum and light curve figures, for in-orbit data quick-view
displayUTC = False

#Maximum ADC value
adcMax = 65536.

#Upper energy limit for plots
energyPlotUpper = 2000.

#------------------------------------------------------------------Plot light curve options------------------------------------------------------------------
#Light curve save root path
lightSaveRootPath = './'

#Light curve save directory name
lightSaveDir = '_lightcurve_all/'

#Light curve save directory name, for seperate-channel plots
lightSaveChDir = '_lightcurve_channel/'

#Plot option for light curve
lightAll = True

#Split runs when plotting light curve
splitRunsLight = False

#Split length when plotting run-split light curve, in seconds
splitRunsLen = -1

#Time step for count rate variation plot, in seconds
timeStep = 1.0

#---------------------------------------------------------------Live time fit and plot options---------------------------------------------------------------
#Cut-off value for live time fit, in seconds
liveTimeCut = 0.1

#Internal MCU dead time
deadTimeMCU = 50e-6

#----------------------------------------------------------General spectrum fit and plot options----------------------------------------------------------
#Available sources
sourceAvailable = ['Am241', 'Ba133', 'Cs137', 'Na22', 'Th228', 'Co60', 'x']

#Reference of nbins for available sources
nbinsRef = {
                'Am241':        2048, 
                'Ba133':        2048, 
                'Th228':        1024, 
                'Cs137':        512, 
                'Na22':         512, 
                'Co60':         512, 
    }

#Reference of bin width for available sources
binWidthRef = {
                'Am241':        32, 
                'Ba133':         32, 
                'Th228':         64, 
                'Cs137':         128, 
                'Na22':          128, 
                'Co60':          128, 
    }

#Reference of fit ranges for some sources
fitRangeRef = {
        'Am241':    [[1100, 2800], [1100, 2300], [1100, 2600], [1200, 2900]], 
        'Ba133':     [[2500, 4500], [2800, 5500], [2800, 5500], [3200, 5500]], 
        'Na22':      [[14000, 19000], [14000, 18000], [15000, 20000], [15000, 21000]], 
        'Cs137':     [[19000, 24000], [18000, 23000], [19000, 26000], [20000, 27000]], 
        'Th228':    [[6000, 9000], [6000, 9000], [6000, 10000], [7000, 10000]], 
        'Co60':      [[41000, 47000], [39500, 45000], [43500, 51000], [46500, 53000]], 
        # 'Co60_low':      [[33000, 41000], [32000, 39000], [35000, 43000], [37000, 44500]], #Currently unused
    }

#Temp-bias correction of fit range
fitRangeCorr = False

#Form of internal background used in spectrum fit
bkgForm = 'lin'

#Path to save fit results
saveFitPath = './fit_result'

#Save fit results settings, whether the data file uses its own save fit results path(create new sub-directories for each new data file)
saveFitOwnPath = False

#Option to import fit result from file
importFitResult = False

#---------------------------------------------------------Advanced spectrum fit and plot options---------------------------------------------------------
#----------------------------------------------------------------------For individual files----------------------------------------------------------------------
#Reference of nbins for individual files
nbinsFileRef = {
                '0-265-5.txt':            8192, 
                '0-267-5.txt':            8192, 
                '0-269-9.txt':            4096, 
                '0-271-9.txt':            4096, 
                '0-273-9.txt':            4096, 
                '0-275-15.txt':            4096, 
                '0-277-15.txt':            4096, 
                '5-265-5.txt':            8192, 
                '5-267-5.txt':            8192, 
                '5-269-9.txt':            4096, 
                '5-271-9.txt':            4096, 
                '5-273-9.txt':            4096, 
                '5-275-15.txt':            4096, 
                '5-277-15.txt':            4096, 
                '10-265-5.txt':            8192, 
                '10-267-5.txt':            8192, 
                '10-269-9.txt':            8192, 
                '10-271-9.txt':            4096, 
                '10-273-9.txt':            4096, 
                '10-275-15.txt':            4096, 
                '10-277-15.txt':            4096, 
                '10-279-15.txt':            4096, 
                '15-265-5.txt':            8192, 
                '15-267-5.txt':            8192, 
                '15-269-9.txt':            8192, 
                '15-271-9.txt':            4096, 
                '15-273-9.txt':            4096, 
                '15-275-15.txt':            4096, 
                '15-277-15.txt':            4096, 
                '15-279-15.txt':            4096, 
                '15-281-15.txt':            4096, 
                '20-265-5.txt':            8192, 
                '20-267-5.txt':            8192, 
                '20-269-9.txt':            8192, 
                '20-271-9.txt':            8192, 
                '20-273-9.txt':            4096, 
                '20-275-15.txt':            4096, 
                '20-277-15.txt':            4096, 
                '20-279-15.txt':            4096, 
                '20-281-15.txt':            4096, 
                '20-283-15.txt':            4096, 
                '25-265-5.txt':            16384, 
                '25-267-6.txt':            8192, 
                '25-269-9.txt':            8192, 
                '25-271-9.txt':            8192, 
                '25-273-9.txt':            4096, 
                '25-275-15.txt':            4096, 
                '25-277-15.txt':            4096, 
                '25-279-15.txt':            4096, 
                '25-281-15.txt':            4096, 
                '25-283-15.txt':            4096, 
                '30-265-5.txt':            16384, 
                '30-267-5.txt':            8192, 
                '30-269-9.txt':            8192, 
                '30-271-9.txt':            8192, 
                '30-273-9.txt':            4096, 
                '30-275-15.txt':            4096, 
                '30-277-15.txt':            4096, 
                '30-279-15.txt':            4096, 
                '30-281-15.txt':            4096, 
                '30-283-15.txt':            4096, 
    }

#Reference of bin width for individual files
binWidthFileRef = {
                'DUMMY':                                 adcMax, 
    }

#Reference of fit range for some specified files
fitFileRef = {
        '0-265-5.txt':          [[400, 1000], [350, 900], [400, 1000], [400, 1000]], 
        '0-267-5.txt':          [[450, 1100], [400, 1000], [400, 1200], [400, 1200]], 
        '0-269-9.txt':          [[500, 1250], [500, 1200], [500, 1250], [500, 1400]], 
        '0-271-9.txt':          [[600, 1400], [600, 1400], [600, 1400], [600, 1600]], 
        '0-273-9.txt':          [[700, 1600], [700, 1600], [700, 1600], [700, 1600]], 
        '0-275-15.txt':          [[600, 1400], [600, 1400], [600, 1400], [600, 1400]], 
        '0-277-15.txt':          [[800, 2000], [900, 2000], [800, 2000], [900, 2200]], 
        '0-279-15.txt':          [[1000, 2200], [1000, 2200], [900, 2200], [1000, 2500]], 
        '0-287-15.txt':          [[1500, 3300], [1500, 3300], [1500, 3300], [1500, 3300]], 
        '0-289-15.txt':         [[1500, 3500], [1500, 3500], [1500, 3500], [1750, 4000]], 
        '5-265-5.txt':          [[350, 900], [350, 900], [300, 900], [400, 900]], 
        '5-267-5.txt':          [[400, 1000], [400, 1000], [400, 1000], [400, 1100]], 
        '5-269-9.txt':          [[500, 1200], [500, 1100], [500, 1200], [500, 1200]], 
        '5-271-9.txt':          [[600, 1300], [500, 1200], [500, 1400], [600, 1400]], 
        '5-273-9.txt':          [[600, 1500], [600, 1400], [600, 1500], [600, 1600]], 
        '5-275-15.txt':          [[800, 1600], [700, 1600], [700, 1700], [700, 1800]], 
        '5-281-15.txt':          [[1000, 2300], [1000, 2200], [1000, 2200], [1000, 2500]], 
        '5-283-15.txt':          [[1000, 2500], [1200, 2500], [1200, 2500], [1200, 2700]], 
        '5-285-15.txt':          [[1300, 2700], [1300, 2700], [1400, 2800], [1500, 3000]], 
        '5-289-15.txt':          [[1500, 3500], [1500, 3500], [1500, 3500], [1500, 3500]], 
        '10-265-5.txt':          [[300, 800], [300, 800], [300, 900], [300, 800]], 
        '10-267-5.txt':          [[400, 900], [400, 900], [300, 1000], [400, 1000]], 
        '10-269-9.txt':          [[400, 1100], [400, 1100], [400, 1200], [400, 1200]], 
        '10-275-15.txt':          [[600, 1400], [600, 1400], [600, 1400], [600, 1400]], 
        '10-281-15.txt':          [[1000, 2200], [1000, 2200], [1000, 2300], [1000, 2500]], 
        '10-289-15.txt':          [[1500, 3300], [1500, 3200], [1500, 3300], [1500, 3500]], 
        '15-265-5.txt':          [[300, 700], [300, 700], [300, 700], [300, 800]], 
        '15-267-5.txt':          [[350, 900], [350, 850], [350, 900], [300, 900]], 
        '15-269-9.txt':          [[400, 1000], [400, 1000], [400, 1000], [400, 1000]], 
        '15-271-9.txt':          [[500, 1200], [500, 1100], [500, 1200], [500, 1300]], 
        '15-273-9.txt':          [[600, 1400], [600, 1300], [600, 1400], [600, 1400]], 
        '15-275-9.txt':          [[600, 1500], [600, 1500], [700, 1600], [600, 1600]], 
        '15-283-15.txt':          [[1000, 2500], [1000, 2300], [1000, 2500], [1000, 2500]], 
        '15-285-15.txt':          [[1200, 2500], [1200, 2500], [1200, 2600], [1300, 2700]], 
        '20-265-5.txt':          [[250, 650], [250, 650], [250, 700], [250, 700]], 
        '20-267-5.txt':          [[300, 800], [300, 800], [300, 800], [300, 850]], 
        '20-269-9.txt':          [[400, 1000], [400, 900], [400, 1000], [400, 1000]], 
        '20-271-9.txt':          [[500, 1100], [500, 1100], [500, 1100], [500, 1200]], 
        '20-279-15.txt':          [[750, 1800], [800, 1800], [800, 1800], [800, 2000]], 
        '20-281-15.txt':          [[900, 2000], [900, 2000], [900, 2000], [1000, 2200]], 
        '20-289-15.txt':          [[1400, 2900], [1500, 3000], [1500, 3000], [1500, 3500]], 
        '25-265-5.txt':          [[200, 600], [200, 600], [200, 650], [250, 650]], 
        '25-267-6.txt':          [[250, 700], [250, 700], [300, 750], [300, 800]], 
        '25-269-9.txt':          [[400, 800], [400, 800], [400, 800], [400, 900]], 
        '25-271-9.txt':          [[400, 1000], [400, 1000], [400, 1000], [400, 1100]], 
        '25-273-9.txt':          [[500, 1200], [500, 1200], [500, 1200], [500, 1200]], 
        '25-275-12.txt':          [[600, 1400], [600, 1400], [600, 1400], [600, 1400]], 
        '25-283-15.txt':          [[900, 2100], [900, 2100], [900, 2100], [1000, 2200]], 
        '25-285-15.txt':          [[1000, 2300], [1000, 2200], [1000, 2400], [1000, 2500]], 
        '30-265-5.txt':          [[200, 550], [200, 500], [200, 600], [200, 600]], 
        '30-267-5.txt':          [[250, 650], [250, 650], [250, 700], [300, 750]], 
        '30-273-9.txt':          [[500, 1000], [500, 1000], [500, 1100], [500, 1200]], 
        '30-275-9.txt':          [[500, 1200], [500, 1200], [500, 1200], [500, 1300]], 
        '30-281-15.txt':          [[800, 1800], [800, 1700], [800, 1800], [800, 1900]], 
        '30-283-15.txt':          [[800, 2000], [800, 2000], [900, 2000], [1000, 2200]], 
        '30-285-15.txt':          [[900, 2300], [1000, 2200], [900, 2300], [1000, 2500]], 
        'Am241-0度-1.dat':       [[[400, 1200], [400, 1200], [400, 1200], [400, 1200]], [[600, 1600], [600, 1600], [600, 1600], [600, 1600]], [], \
            [[], [], [], [1000, 3000]], [[1000, 3000], [1000, 3000], [1000, 3000], [1000, 3500]], [], []],
        'Am241-10度-1.dat':      [[[], [], [], [400, 1200]], [], [[], [], [], [800, 2500]], [], [], [[], [], [], [1000, 3500]], [[], [], [], [1500, 4000]]],
        'Am241-20度-1.dat':      [[[300, 900], [300, 900], [300, 1000], [400, 1000]], [[], [], [], [500, 1400]], [[500, 1800], [600, 1800], [500, 2000], [700, 2000]], \
            [], [], [], [[1000, 3000], [1000, 3000], [1000, 3000], [1200, 3500]]],
        'Am241-25度-1.dat':      [[[300, 800], [300, 800], [300, 900], [300, 1000]], [[400, 1200], [400, 1200], [400, 1300], [500, 1400]], \
            [[600, 1600], [600, 1600], [600, 1700], [600, 1800]], \
            [[], [700, 2000], [], []], [[700, 2200], [700, 2200], [800, 2300], [1000, 2500]], [[800, 2100], [800, 2100], [800, 2200], [1000, 2500]], [[], [], [], []]],
        'Am241-30度-1.dat':      [[[200, 800], [200, 800], [300, 800], [300, 900]], [[], [400, 1100], [], []], [[500, 1500], [600, 1500], [500, 1600], [600, 1800]], \
            [[600, 1800], [600, 1800], [700, 2000], [700, 2200]], [], [], []],
        'Am241CI-0度-0.dat':      [[[500, 1100], [500, 1100], [], [500, 1300]], [[], [], [], [750, 1750]], [], [], [], [], []],
        'Am241CI-10度-0.dat':      [[[400, 1000], [400, 1000], [400, 1100], [500, 1200]], [[], [], [], [500, 1500]], [], [], [], [], []],
        'Am241CI-20度-0.dat':      [[[300, 900], [300, 900], [300, 1000], [400, 1000]], [[500, 1200], [], [], [700, 1400]], [], [], [], [], []],
        'Am241CI-25度-0.dat':      [[[300, 800], [300, 800], [300, 900], [300, 1000]], [[400, 1200], [500, 1200], [500, 1300], [700, 1400]], [], [], [], [], []],
        'Am241CI-30度-0.dat':      [[[200, 800], [200, 800], [300, 900], [300, 1000]], [[400, 1100], [400, 1100], [400, 1200], [400, 1300]], [], [], [], [], []],
        'Am241CI-35度-0.dat':      [[[200, 700], [200, 700], [200, 800], [200, 800]], [[300, 1000], [300, 1000], [300, 1100], [400, 1200]], [], [], [], [], []], 
    }

#-----------------------------------------------------------Temperature-bias correction options-----------------------------------------------------------
#Standard temperature for temperature-bias correction, in degrees celcius
tempStandard = 25.0

#Standard bias for temperature-bias correction, in volts
biasStandard = 28.5

#----------------------------------------------------------Plot telemetry data variation options----------------------------------------------------------
#Telemetry plot settings. group by scans or channels
telemetryGroupByScans = False

#----------------------------------------------------------Temperatue-bias fit and plot options----------------------------------------------------------
#Temperature-bias fit settings for contour plot
tempBiasCont = False

#Path to save temp-bias fit results and figures
saveTempBiasPath = './tempbias_result'

#Plot 3D surface or temperature/bias split curves
tempBiasPlot3D = True

#Split data by set bias or temperature values when doing split curves plot
tempBiasSplitBias = True

#Systematic error for temperature measurement, in degrees celcius
tempErrSys = 0.5

#Set bias values, for bias-split plot
biasSetVal = [26.5, 26.7, 26.9, 27.1, 27.3, 27.5, 27.7, 27.9, 28.1, 28.3, 28.5, 28.7, 28.9]

#Range for choosing data points for each set bias value
biasSetStd = 0.05

#Set temperature values, for temperature-split plot
tempSetVal = [0, 5, 10, 15, 20, 25, 30]

#Range for choosing data points for each set temperature value
tempSetStd = 1.

#-------------------------------------------------------------Leak current fit and plot options-------------------------------------------------------------
#Path to save leak current fit results and figures
saveCurrentPath = './current_result'

#---------------------------------------------------------------Angular responce plot options---------------------------------------------------------------
#Angular response settings, scan interval in degrees
angleScanIntv = 15

#Path to save angular response results and figures
saveAnglePath = './angle_result'

#-----------------------------------------------------------Energy-channel fit and plot options-----------------------------------------------------------
#Reference of sources' energies, in keV
sourceEnergyRef = {
        'Am241':    59.5,
        'Na22':      511.,
        'Cs137':     661.7,
        'Th228':    239.,
        'Co60':      1332.2,
}

#Path to save EC plot and fit results and figures
saveECPath = './ec_result'

#Import processed EC data from external file
importECData = False

#Fit EC settings
fitEC = True

#Fit EC settings, fit EC curve as energy-channel(False) or channel-energy(True)
fitECAsEnergy = False

#EC and resolution split energy point, in keV
energySplitEC = [49., 55.]

#Energy of Gd K edge, in keV
energyEdge = 50.2

#EC gain correction
doGaiCorrEC = True

#EC correction factor for gain drop of source data
gainCorrEC = [1.1116, 1.0971, 1.1042, 1.1352]

#Error of EC correction factor for gain
gainCorrECErr = [2.7756e-3, 2.1635e-3, 3.1506e-3, 2.2561e-3]

#Bounds for finding roots as upper and lower bounds of medium energy range
channelBoundRef = [1000, 3000]

#-----------------------------------------------------------------HPGe data process options-----------------------------------------------------------------
#Path to save figures for HPGe data
saveHPGeFigPath = './figure_HPGe'

#------------------------------------------------------------Efficiency process and fit options-------------------------------------------------------------
#HPGe simulation data
HPGeSimuFile = './config_files/02/HPGe_eff_simu.npy'

#GRID simulation data, energy
simuFile = './config_files/02/simu_eff.npy'

#Path to save figures for efficiency data
saveEffFigPath = './eff_result'

#------------------------------------------------------------------Plot stage figure options------------------------------------------------------------------
#Path to save figures for GRID and HPGe stage plot
saveStagePath = './stage_fig'
