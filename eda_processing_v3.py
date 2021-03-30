###############################################################################
# IMPORTAR MODULOS

import sys
import matplotlib.pyplot as plt
import math



###############################################################################
# FUNCOES

def FreqToTime(samples,freq):

    period = float(1/freq)
    time = [i*period for i in samples]
    return time


def EDAConversion(value):

    n = 10 # in BITalino the first four channels are sampled using 10-bit resolution

    # USE THIS BLOCK FOR NORMAL BITALINO
    # source: https://bitalino.com/datasheets/EDA_Sensor_Datasheet.pdf
    #R = 1-(value/(2**n)) #MOhm
    #EDA = float(1/R) #uS

    # USE THIS BLOCK FOR BITALINO REVOLUTION
    # source: http://bitalino.com/datasheets/REVOLUTION_EDA_Sensor_Datasheet.pdf  
    vcc = 3.3 #value in V
    EDA = float(((value/(2**n))*vcc)/0.132) # value in uS

    return EDA


def ReadData(filepath):

    print("\n\tExtracting data...")

    file = open(filepath,mode='r',newline='\r\n') # opens file
    data = file.readlines() # extracts the lines of the file for a list ex: ['2\t4.8\t0.1\n','1\t3.2\t7\n']
    file.close() # closes file

    dataframe = [line.split() for line in data] # separates columns for each line ex: ['2', '4.8', '0.1']
    eda = []
    samples = []
    for i in dataframe:
        samples.append(int(i[0]))
        eda.append(float(i[1]))
    dim = len(samples)
    print("\tFound " + str(dim) + " lines of data")

    return samples,eda,dim


def MovingAverage(y,dim,k,n):
    # k is the window size - base is 0.01
    # n is the step

    # criar uma lista com o EDA dentro da amplitude criada, para nao dar erros de dominio no array
    add_before = [y[0]]*k
    add_after = [y[-1]]*k
    y_new = add_before + y + add_after # concatena as 3 listas

    mov_avg = []
    for i in range(k,dim+k,n):
        y_sum = 0
        for j in range(-k,k+1):
            y_sum += y_new[i+j]
        y_sum = y_sum/((2*k)+1)
        mov_avg.append(round(y_sum,3))

    return mov_avg


def ExtractionSCR(eda,dim):

    window_pc = float(input("\tWindow size percentage >>> "))
    window = int(dim*window_pc)
    print('\tWindow size:',window)
    n = 1

    print('\tApplying moving average method for SCL calculation...')
    test = [i for i in eda]
    scl = MovingAverage(test,dim,window,n)

    scr = [round((eda[i]-scl[i]),5) for i in range(len(scl))]

    return scl,scr


def FindPeaks(x,y,k):

    k = k-k%2+1
    kHalf = int(math.trunc(k/2))

    new_y = y[kHalf+1:]
    new_y = new_y + [-math.inf]*int(kHalf+1)
    lastComparedValue = [-math.inf] * kHalf
    leftHalf = y[0:int(kHalf+1)]

    Max = []
    Min = []
    TMax = []
    TMin = []
    peaks_info = []

    max_num = 0
    min_num = 0
    for i in range(len(new_y)-1):

        curVal = leftHalf.pop(0)
        validation = False

        if curVal > max(lastComparedValue) and curVal >= max(leftHalf) and min_num > max_num:
            # ensures that the first value is a minimum
            max_num += 1
            Max.append(curVal)
            TMax.append(x[i])
            peaks_info.append(["max",x[i],curVal])
        if curVal < min(lastComparedValue) and curVal <= min(leftHalf) and min_num == max_num:
            min_num += 1
            Min.append(curVal)
            TMin.append(x[i])
            peaks_info.append(["min",x[i],curVal])

        lastComparedValue.pop(0)
        lastComparedValue.append(curVal)
        leftHalf.append(new_y[i])

    # ensures that the last value is a minimum
    if peaks_info[-1][0] == 'max':
        peaks_info.pop(-1)
        Max.pop(-1)
        TMax.pop(-1)

    return peaks_info, Max, TMax, Min, TMin


def AnalysePeaks(Max, TMax, Min, TMin,AmpVal,data):

    #Filtrar por amplitude maior ou igual a AmpVal
    Max2 = []
    Min2 = []
    TMax2 = []
    TMin2 = []
    Amplitude = []

    # save points with intendent amplitudes
    for i in range(len(Max)):
        #print(i)
        Amp = float(Max[i]-Min[i])
        #print(Amp)
        if abs(Amp) >= AmpVal:
            TMax2.append(TMax[i])
            Max2.append(Max[i])
            TMin2.append(TMin[i])
            Min2.append(Min[i])
            Amplitude.append([(TMax[i],Max[i]),round(Amp,4)])
    if len(TMin)>len(TMax):
        TMin2.append(TMin[-1])
        Min2.append(Min[-1])

    # calculate rise time and recovery time for those points
    RiseTime = []
    RecoveryTime = []
    for i in range(len(Max2)):
        Amp = float(Max[i]-Min[i])
        # calculate rise time
        rise = TMax2[i]-TMin2[i]
        RiseTime.append([(TMax2[i],Max2[i]),round(rise,4)])
        # calculate recovery time for 63pc decreasing
        halfamp = Amp*(1-0.63)
        difval = [0,math.inf]
        i0 = data[1].index(TMax2[i])
        i1 = data[1].index(TMin2[i+1])
        for j in range(i0,i1+1):
            dif = abs(data[2][j])-halfamp
            if dif < difval[1]:
                difval[0] = data[1][j]
                difval[1] = dif
        halfampTime = difval[0]
        recovery = halfampTime - TMax2[i]
        RecoveryTime.append([(TMax2[i],Max2[i]),round(recovery,4)])

    return Max2, TMax2, Min2, TMin2, Amplitude, RiseTime, RecoveryTime


def Downsampling(x,y):

    new_y = []
    new_x = []
    factor = 170
    for i in range(0,len(y),factor):
        new_y.append(y[i])
        new_x.append(x[i])

    return new_x, new_y


def GraphicFilter(time, eda_raw, eda_movavg):
    fig = plt.figure()
    plt.title('Electrodermal Activity (EDA)', fontsize = 40)
    plt.xlabel('Time (s)', fontsize = 30)
    plt.ylabel("EDA (uS)", fontsize = 30)
    plt.plot(time, eda_raw, '-', color='#0E2A47', label = 'Raw')
    plt.plot(time, eda_movavg, '-', color='#00ffffff', label = 'Filtered')
    plt.grid()
    plt.legend(fontsize = 30)
    plt.show(block=True)


def GraphicSCL(time, eda_movavg, scl):
    fig = plt.figure()
    plt.title('Skin Conductance Level (SCL)', fontsize = 40)
    plt.xlabel('Time (s)', fontsize = 30)
    plt.ylabel("SCL (uS)", fontsize = 30)
    plt.plot(time, eda_movavg, '-', color='#0E2A47', label = 'EDA')
    plt.plot(time, scl, '-', color='#e613d5', label = 'SCL')
    plt.grid()
    plt.legend(fontsize = 30)
    plt.show(block=True)


def GraphicSCR(time, eda_movavg, scr):
    fig = plt.figure()
    plt.title('Skin Conductance Response (SCR)', fontsize = 40)
    plt.xlabel('Time (s)', fontsize = 30)
    plt.ylabel("SCR (uS)", fontsize = 30)
    plt.plot(time, scr, '-', color='#0E2A47')
    plt.grid()
    plt.show(block=True)


def GraphicPeaks(time, Max, Min, TMax, TMin):
    fig = plt.figure()
    plt.title('Skin Conductance Response (SCR)', fontsize = 40)
    plt.xlabel('Time (s)', fontsize = 30)
    plt.ylabel("SCR (uS)", fontsize = 30)
    plt.plot(time, scr, '-', color='#0E2A47')
    plt.plot(TMax, Max, 'o', color='#e613d5')
    plt.plot(TMin, Min, 'o', color='#e613d5')
    plt.grid()
    plt.show(block=True)


def Display():
    display = input("\n\tDisplay graphics? (y - yes; n - no) >>> ")
    while display not in ('y','n'):
        print("\tInvalid character!")
        display = input("\n\tDisplay graphics? (y - yes; n - no) >>> ")
    return display


def Menu():

    option = str(input("\nWhat do you want to do? (a) read data, (b) filter EDA, (c) extract SCR, (d) analyse peaks, (q)uit >>> ")).lower()
    while option not in ('a','b','c','d','q'):
        print("Invalid character!")
        option = str(input("\nWhat do you want to do? (a) read data, (b) filter EDA, (c) extract SCR, (d) analyse peaks, (q)uit >>> ")).lower()

    return option





###############################################################################
# MAIN

print("\nWelcome to the Hei-lab EDA processing tool!")

option = Menu() # executar menu inicial

while option != 'q':


    ######################################################
    # Reading raw data
    if option == 'a':

        filepath = './eda_raw.csv'
        freq = int(input('\n\tSignal frequency (Hz) >>> '))

        # executar funcao 'ReadData' para extracao de dados - mudar localizacao e nome do ficheiro de acordo com o desejado
        samples,eda_raw,dim = ReadData(filepath)
        time = FreqToTime(samples,freq)
        print('\n\tData size (time and signal dimension):',len(time),',',dim)

        # DISPLAY GRAPHICS
        if Display() == 'y':
            fig = plt.figure()
            plt.title('EDA raw', fontsize = 40)
            plt.xlabel('Samples', fontsize = 30)
            plt.ylabel("EDA (uS)", fontsize = 30)
            plt.plot(samples, eda_raw, '-', color='#0E2A47', linewidth = 8)
            plt.grid()
            plt.tick_params(labelsize = 25)
            plt.show(block=True)

        option = Menu() # asks for new action


    ######################################################
    # FILTERING RAW DATA (MOVING AVERAGE)
    if option == 'b':

        print('\n\tSmoothing raw data with Moving Average...')

        window_pc = float(input("\tWindow size percentage >>> "))
        window = int(dim*window_pc)
        print('\tWindow size:',window)
        #n = int(dim*0.001)
        n = 1
        eda_movavg = MovingAverage(eda_raw,dim,window,n)
        print("\tFiltered data size:",len(eda_movavg))

        # SAVE DATA
        if str(input("\n\tSave data? (y - yes; n - no) >>> ")) == 'y':
            eda_filtered = open('./eda_movavg.csv','w')
            for i in range(dim):
                eda_filtered.write(str(time[i]) + '\t' + str(eda_movavg[i]) + '\n')
            eda_filtered.close()
            print ("\tData saved")

        # DISPLAY GRAPHICS
        if Display() == 'y':
            GraphicFilter(time,eda_raw,eda_movavg)
            print ("\tGraphic created.")

        option = Menu() # asks for new action


    ######################################################
    # EXTRACTING SCL AND SCR FROM MOVING AVERAGE
    if option == 'c':

        print('\n\tExtracting SCL and SCR...')
        scl,scr = ExtractionSCR(eda_movavg,dim)
        print('\tSCL size:',len(scl))
        print('\tSCR size:',len(scr))

        # SAVE DATA
        if str(input("\n\tSave data? (y - yes; n - no) >>> ")) == 'y':
            scl_data = open('./scl.csv','w')
            scr_data = open('./scr.csv','w')
            for i in range(dim):
                scl_data.write(str(time[i]) + '\t' + str(scl[i]) + '\n')
                scr_data.write(str(time[i]) + '\t' + str(scr[i]) + '\n')
            scl_data.close()
            scr_data.close()
            print ("\tData saved")

        # DISPLAY GRAPHICS
        if Display() == 'y':
            GraphicSCL(time, eda_movavg, scl)
            GraphicSCR(time, eda_movavg, scr)
            print ("\tGraphic created.")

        option = Menu() # asks for new action


    ######################################################
    # ANALYSING PEAKS
    if option == 'd':

        print('\n\tAnalysing peaks...')
        window_pc = float(input("\tWindow percentage >>> "))
        window = int(dim*window_pc)
        print('\n\tWindow size:',window)
        peaks_info, Max, TMax, Min, TMin = FindPeaks(time,scr,window)
        amplitude = float(input("\n\tMinimum peak amplitude for SCR >>> "))
        data = [[i for i in range(len(time))],[i for i in time], [i for i in scr]]
        Max2, TMax2, Min2, TMin2, Amplitude, RiseTime, RecoveryTime = AnalysePeaks(Max, TMax, Min, TMin,amplitude,data)

        print('\n\t> Peak anaysis results')
        for i in range(len(Amplitude)):
            print('\tPeak #' + str(i+1) + '\n\tMax: ' + str(Amplitude[i][0][1]) + ', Amplitude: ' + \
                  str(Amplitude[i][1]) + ', Rise Time: ' + str(RiseTime[i][1]) + ', Recovery Time: ' + \
                  str(RecoveryTime[i][1]))

        # SAVE DATA
        if str(input("\n\tSave data? (y - yes; n - no) >>> ")) == 'y':
            print("\tSaving peaks info...")
            peaks_info_data = open('./peaks_info.csv','w')
            for i in peaks_info:
                peaks_info_data.write(str(i) + '\n')
            peaks_info_data.close()
            print("\tPeaks info saved.")

        # DISPLAY GRAPHICS
        if Display() == 'y':
            GraphicPeaks(time, Max, Min, TMax, TMin)
            GraphicPeaks(time, Max2, Min2, TMax2, TMin2)
            print ("\tGraphic created.")

        option = Menu() # asks for new action


print("Program finished.")
sys.exit() #simulador termina


