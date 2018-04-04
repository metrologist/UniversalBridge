from __future__ import print_function
from __future__ import division

"""
ubridge.py creates a UNIVERSALBRIDGE class that includes inbuilt uncertainties due to the bridge design and takes in
updated calibration data for the bridge, including calibration uncertainties.
bridge_value() takes dial readings from the Universal Impedance Bridge and returns the corrected impedance value.
All calibration work with the UB should use this module.
"""
import GTC as gtc  # likely to return values as ureal or ucomplex
import math
import csv


class UNIVERSALBRIDGE(object):
    """
    All methods related to the universal bridge are held in this class.
    """

    def __init__(self, path_name, temp):
        self.path_name = path_name
        self.temp = temp
        self.ranges = ['1', '2', '3', '4', '5', '6', '7']  # UB range switch
        self.config = ['Z', 'Y']  # either impedance or admittance configuration
        self.standard_names = ['C1', 'R4A', 'R4B', 'R4C', 'G1', 'G2']  # internal standards
        self.standard_nominal = {'C1': 1e-9, 'R4A': 1e4, 'R4B': 1e5, 'R4C': 1e6, 'G1': 1e-5, 'G2': 1e-5}
        self.caldata, self.acdc, self.tempcoeffs, self.stability, self.ivd = self.get_csv_data(self.path_name)
        #  set up gtc values of the internal components
        self.R4A = self.standard_nominal['R4A'] * (1 + self.caldata['R4A']) * (1 + self.acdc['R4A']) * \
                   (1 + gtc.function.mul2(self.tempcoeffs['R4A'], (self.temp - self.caldata['cal_temp']))) * (
                       1 + self.stability['R4A'])
        self.R4B = self.standard_nominal['R4B'] * (1 + self.caldata['R4B']) * (1 + self.acdc['R4B']) * \
                   (1 + gtc.function.mul2(self.tempcoeffs['R4B'], (self.temp - self.caldata['cal_temp']))) * (
                       1 + self.stability['R4B'])
        self.R4C = self.standard_nominal['R4C'] * (1 + self.caldata['R4C']) * (1 + self.acdc['R4C']) * \
                   (1 + gtc.function.mul2(self.tempcoeffs['R4C'], (self.temp - self.caldata['cal_temp']))) * (
                       1 + self.stability['R4C'])
        self.C1 = self.standard_nominal['C1'] * (1 + self.caldata['C1']) * (1 + self.stability['C1']) * \
                  (1 + gtc.function.mul2(self.tempcoeffs['C1'], (self.temp - self.caldata['cal_temp'])))
        self.G1 = self.standard_nominal['G1'] * (1 + self.caldata['G1']) * (1 + self.acdc['G1']) * \
                  (1 + gtc.function.mul2(self.tempcoeffs['G1'], (self.temp - self.caldata['cal_temp']))) * (
                      1 + self.stability['G1'])
        self.G2 = self.standard_nominal['G2'] * (1 + self.caldata['G2']) * (1 + self.acdc['G2']) * \
                  (1 + gtc.function.mul2(self.tempcoeffs['G2'], (self.temp - self.caldata['cal_temp']))) * (
                      1 + self.stability['G2'])
        self.amp = (1 + self.caldata['Amp'])  # * (1 + 1j * self.caldata['Amp_phase']) # Amp_phase is frequency scaled
        self.res_aux = 1e-5  # low Z ranges have 10 ppm uncertainty due to difficulty of managing auxiliary balances

        # the unit_scales dictionaries hold the dial multiplier for each range
        self.x_unit_scales = {'1Z': 1e-10, '2Z': 1e-9, '3Z': 1e-8, '4Z': 1e-7, '5Z': 1e-6, '6Z': 1e-5, '7Z': 1e-4,
                              '1Y': 1e-10, '2Y': 1e-11, '3Y': 1e-12, '4Y': 1e-13, '5Y': 1e-14, '6Y': 1e-15, '7Y': 1e-16}
        self.r_unit_scales = {'1Z': 1e-7, '2Z': 1e-6, '3Z': 1e-5, '4Z': 1e-4, '5Z': 1e-3, '6Z': 1e-2, '7Z': 1e-1,
                              '1Y': 1e-7, '2Y': 1e-8, '3Y': 1e-9, '4Y': 1e-10, '5Y': 1e-11, '6Y': 1e-12, '7Y': 1e-13}
        self.scale_multiplier = {'1': 1e4, '2': 1e3, '3': 1e2, '4': 1e1, '5': 1e0, '6': 1e0, '7': 1e0}  # tran ratios
        self.xdialmax = 1000000  # LC dial is 6 dials, each reading -1 to 10
        self.rdialmax = 10000000  # RG dial is 7 dials, each reading -1 to 10
        self.Adial_lin = self.rdialmax * self.ivd['Adial_lin']
        self.Bdial_lin = self.xdialmax * self.ivd['Bdial_lin']

    def get_csv_data(self, path):
        """
        :param path: the csv file formatted to hold all the calibration dictionaries
        :return: four dictionaries caldata, acdc, tempcoeffs, stability, ivd
        The caldata dictionary contains the calibrated component values, acdc contains frequency dependence uncertainty,
        tempcoeffs has estimates of termperature coefficients, stability has uncertainties for time dependence and ivd
        has linearity estimates for the two main inductive voltage dividers.
        """
        cal_file = open(path, mode='rb')  # Open in 'binary' mode
        reader = csv.reader(cal_file)
        reader.next()  # read in unused row
        header = reader.next()  # read in unused row
        table = gtc.number_strings.sequence_parser(reader)
        cal_file.close()
        Nr = len(table)
        # Nc = len(header)
        caldata = {}
        acdc = {}
        tempcoeffs = {}
        stability = {}
        ivd = {}
        for r in range(Nr):
            dictionary = table[r][0]
            assert dictionary in ['caldata', 'acdc', 'tempcoeffs', 'stability',
                                  'ivd'], 'csv dictionary name not in list'
            key = table[r][1]
            value = gtc.ureal(table[r][2], table[r][3], table[r][4], table[r][5])
            if dictionary == 'caldata':
                caldata[key] = value
            elif dictionary == 'acdc':
                acdc[key] = value
            elif dictionary == 'tempcoeffs':
                tempcoeffs[key] = value
            elif dictionary == 'stability':
                stability[key] = value
            elif dictionary == 'ivd':
                ivd[key] = value
        return caldata, acdc, tempcoeffs, stability, ivd

    def unit_scale(self, ubrange):
        """
        range is entered as '4Z' and a dial factor is returned in a tuple with the unit
        for internal use with available_ranges
        """
        assert len(ubrange) == 2, "invalid range identifier, must be exactly two characters long"
        assert ubrange[0] in self.ranges, "first range identifier character must be in 1...7"
        assert ubrange[1] in self.config, "second range identifier character must by Y or Z"
        return self.x_unit_scales[ubrange], self.r_unit_scales[ubrange]

    def available_ranges(self):
        """
        Prints out the ranges on the Universal Impedance Bridge. For information only.
        Useful to cross-check that main calculation is correct at full-scale dial settings.
        """
        range_listZ = []
        range_listY = []
        for x in self.ranges:
            range_listZ.append(str(x) + 'Z')
            range_listY.append(str(x) + 'Y')
        for x in range_listZ:
            a = self.unit_scale(x)
            maxL = a[0] * self.xdialmax
            maxR = a[1] * self.rdialmax
            print(x, maxL, 'H', '\t', maxR, 'ohm')
        for x in range_listY:
            a = self.unit_scale(x)
            maxC = a[0] * self.xdialmax
            maxG = a[1] * self.rdialmax
            print(x, maxC, 'F', '\t', maxG, 'siemen')

    def bridge_value(self, ubrange, dialA, dialB, freq, res):
        """
        ubrange as usual 5Z or 2Y format
        dialA is the resistance dial, dialB is the L/C dial
        returns either impedance or admittance in ohm (F, H, S); i.e. no assumed multiplier
        note that the wC1 can have w dropped to let results be directly in L or F, but could have 2nd order effect
        with real/imag interaction
        """
        # add linearity uncertainty to dial settings, but create as new ureal so subtracting zero reading does not cancel the linearity error
        adial_lin = gtc.ureal(self.Adial_lin.x, self.Adial_lin.u, self.Adial_lin.df, 'Adial_lin' + repr(dialA))
        dialA = dialA + adial_lin
        bdial_lin = gtc.ureal(self.Bdial_lin.x, self.Bdial_lin.u, self.Bdial_lin.df, 'Bdial_lin' + repr(dialB))
        dialB = dialB + bdial_lin
        assert res in [0, 1], "res must be 0 if off or 1 if used"
        if res == 1:  # use resolution term unless actual scatter of dial readings available
            resolution = (2e3 / freq) ** (1.0 / 1.4) * 1e-6  # 1 ppm at 2 kHz, increasing to 6 ppm at 50 Hz
        else:
            resolution = 0  # the assumption is that you will then take a mean and sd of several measured values
        Adial_res = gtc.ureal(0, self.rdialmax * resolution, 10,
                              'Adialres' + repr(dialA))  # label identifies dial setting
        assert ubrange[0] in self.ranges, "first range identifier character must be in 1...7"
        assert ubrange[1] in self.config, "second range identifier character must by Y or Z"
        if ubrange[0] < 4 and ubrange[1] == 'Z':  # i.e. 1Z, 2Z or 3Z gets an aux resolution term
            aux_resolution = self.res_aux  # aux_resolution is a 10 ppm term to cover the difficulty of balancing I for R
        else:
            aux_resolution = 0.0  # only Z ranges 1, 2 and 3 present a challenge with repeatable auxiliary balances be cause of lead/switch resistance drift
        Adial_res = Adial_res + gtc.ureal(0, self.rdialmax * aux_resolution, 10, 'Adialauxbal' + repr(dialA))
        Bdial_res = gtc.ureal(0, self.xdialmax * resolution, 10, 'Bdialres' + repr(dialB))  # in a gtc budget
        dialA = dialA + Adial_res
        dialB = dialB + Bdial_res
        # select relevant R4 based on range
        range_number = ubrange[0]
        range_type = ubrange[1]
        if int(range_number) == 7:
            R4 = self.R4C
        elif int(range_number) == 6:
            R4 = self.R4B
        else:
            R4 = self.R4A
        # just to keep the formula short, the calibration errors are made local
        C1, G1, G2, amp = self.C1, self.G1, self.G2, self.amp
        # scale dials to give max of 1.1 for a and b
        a = dialA / self.rdialmax
        b = dialB / self.xdialmax
        phase_error = (1 + 1j * self.caldata['Amp_phase'] * freq / 1600)  # Amp_phase measured at 1.6 kHz
        # calculate impedance/admittance
        if range_type == 'Z':
            zz = R4 / G2 * (a * G1 + 1j * 2 * math.pi * freq * b * C1) * amp * phase_error / self.scale_multiplier[
                range_number]
            return zz
        if range_type == 'Y':
            G3 = 1 / R4
            yy = G3 / G2 * (a * G1 + 1j * 2 * math.pi * freq * b * C1) / amp * phase_error * self.scale_multiplier[
                range_number]  # note division by amp, the RHS bridge error
            return yy

    def cmc_uncert(self, ubrange, dialA, dialB, frequency):
        """
        The minimum uncertainty that can be quoted in a report under the CIPM MRA is determined by the officially
        recognised CMCs. These CMCs can only be an approximation to what is correctly calculated by bridge_value.
        :param ubrange: 3Y, 7Z style two character names for the rangee
        :param dialA: the resistance 7 dials
        :param dialB:  the reactance 6 dials
        :param frequency: in Hz
        :return: a tuple of real and reactive uncertainties, (g_uncert, c_uncert) or (ronly_uncert, r_uncert, l_uncert),
        in base units. The ronly_uncert is for a resistor, while the r_uncert is for the resistive part of an inductor.
        """
        immittance = self.bridge_value(ubrange, dialA, dialB, frequency, 1)  # note these methods could be combined
        range_type = ubrange[1]
        assert range_type in self.config, "second range identifier character must by Y or Z"
        if range_type == 'Y':
            capacitance = immittance.imag.x / (2 * math.pi * frequency)  # capacitance in farad
            capmicro = capacitance * 1e6  # capacitance in microfarad
            c_uncert = (0.2 / frequency + 22 * capmicro) * 1e-12  # convert from pF to F
            # conductance = immittance.real.x # not needed in CMC calculation
            cappf = capacitance * 1e12  # capacitance in picofarad
            if cappf == 0.0:
                cappf = 1e-6  # use 1 aF instead of zero to calculate tand
            tand = (0.000027 + 0.00027 / cappf)
            g_uncert = 2 * math.pi * frequency * capacitance * tand  # G =wC*tan(d)
            return g_uncert, c_uncert
        else:
            range_number = int(ubrange[0])  # will be 1...7
            inductance = immittance.imag.x / (2 * math.pi * frequency)  # inductance in henry
            inductmicro = inductance * 1e6  # inductance in microhenry
            resistance = immittance.real.x
            ronly_uncert = (2000 / frequency + 19 * resistance) * 1e-6
            factor1 = 0.2 * 10 ** (range_number - 1)  # (range_number-1) multiplies by 10 for each step up in range
            l_uncert = math.sqrt((0.000014 * inductmicro) ** 2 + (0.001 * resistance) ** 2 +
                                 (factor1 / frequency) ** 2) * 1e-6  # in microhenry
            factor2 = 0.002 * 10 ** (range_number - 1)  # (range_number-1) multiplies by 10 for each step up in range
            r_uncert = math.sqrt((0.000014 * resistance) ** 2 + (0.0000001 * inductmicro) ** 2 +
                                 (factor2 / frequency) ** 2)  # in ohm
            return ronly_uncert, r_uncert, l_uncert


if __name__ == "__main__":
    temperature = gtc.ureal(20, 0.5, 10, 'temperature')  # temperature at which bridge is being used
    calfile = r'ubdict_nov_2017.csv'  # the file that contains practically all calibration and uncertainty information
    # create the bridge object
    ubridge = UNIVERSALBRIDGE(calfile, temperature)  # temperature put here as it might be in common with the UUT

    print('Usage example, returning a value and full uncertainty budget')
    ##     answer = ubridge.bridge_value('1Z', 101190, 118, 1.59135e3, 1)
    answer = ubridge.bridge_value('7Y', 46, 10004, 1.59135e3, 1)
    print(repr(answer.real))
    print(repr(answer.imag / (2 * math.pi * 1.59135e3)))
    for l, u in gtc.reporting.budget(answer.real, trim=0): print(l, u)
    print()

    print('Calculate values for the dials set to full scale')
    zranges = ['1Z', '2Z', '3Z', '4Z', '5Z', '6Z', '7Z']
    yranges = ['1Y', '2Y', '3Y', '4Y', '5Y', '6Y', '7Y']
    freq = 1.6e3
    for x in zranges:
        z = ubridge.bridge_value(x, 10000000, 1000000, freq, 0)
        print(x, z.imag / (2 * math.pi * freq), 'H', z.real, 'ohm')
    for x in yranges:
        y = ubridge.bridge_value(x, 10000000, 1000000, freq, 0)
        print(x, y.imag / (2 * math.pi * freq), 'F', y.real, 'siemen')
    print()

    print('Compare with the UB\'s nameplate range values as, embedded in available_ranges')
    ubridge.available_ranges()
    print()

    print('Check standards have been correctly constructed, including uncertainties')
    all_standards = [ubridge.R4A, ubridge.R4B, ubridge.R4C, ubridge.G1, ubridge.G2, ubridge.C1, ubridge.amp]
    for x in all_standards:
        print(repr(x))
        for l, u in gtc.reporting.budget(x, trim=0): print(l, u)
        print()

    print('Compare calculated uncertainty and CMCs for different dial settings')
    rdials = [0, 1000, 10000, 10000000]
    xdials = [0, 100, 1000, 1000000]
    check_range = '3Z'  # choose range to explore
    frequency = 53  # choose frequency
    for r in rdials:
        for x in xdials:
            result = ubridge.bridge_value(check_range, r, x, frequency, 1)
            real_part = result.real
            reactive_part = result.imag / (2 * math.pi * frequency)
            cmc = ubridge.cmc_uncert(check_range, r, x, frequency)
            # print(r, x, real_part, reactive_part, cmc )
            if check_range[1] == 'Y':
                print(r, x, cmc[0] / real_part.u, cmc[1] / reactive_part.u)
            else:
                print(r, x, cmc[0] / real_part.u, cmc[1] / real_part.u, cmc[2] / reactive_part.u)

    dict_caldata, dict_acdc, dict_tempcoeffs, dict_stability, dict_ivd = ubridge.get_csv_data(calfile)
    component_list = ubridge.standard_names
    print('caldata', component_list)
    for x in component_list:
        print(repr(dict_caldata[x]))
    print(repr(dict_caldata['Amp']))
    print(repr(dict_caldata['Amp_phase']))
    print(repr(dict_caldata['cal_temp']))
    print('acdc', component_list)
    for x in component_list:
        print(repr(dict_acdc[x]))
    print('tempcoeffs', component_list)
    for x in component_list:
        print(repr(dict_tempcoeffs[x]))
    print('stability', component_list)
    for x in component_list:
        print(repr(dict_stability[x]))
    other_list = ['Adial_lin', 'Bdial_lin']
    print(other_list)
    for x in other_list:
        print(repr(dict_ivd[x]))

    print(repr(ubridge.C1))
"""
C:\Python27\python.exe Y:/Staff/KJ/PycharmProjects/UBridge/ubridge.py
Usage example, returning a value and full uncertainty budget
ureal(0.0101197133986692406477,1.34859840260857095715e-06,20.2)
ureal(1.180233177708003485054596e-08,1.351174443958782163395534e-10,20.4)
Bdial_lin118 7.07157888972e-07
Adial_lin101190 7.07156632855e-07
Bdialresureal(118,1,10) 6.30789197227e-07
Adialresureal(101190,10,10) 6.30788076763e-07
Amp_phase_err 1.20997776339e-07
Amperr 1.00185760757e-07
R4Aacdc 1.78905112144e-08
G2acdc 1.78905112144e-08
G1acdc 1.78892949197e-08
G2alpha 7.15620448576e-09
G2stab 7.15620448576e-09
R4Aalpha 7.15620448576e-09
R4Astab 7.15620448576e-09
G1alpha 7.15571796786e-09
G1stab 7.15571796786e-09
temperature 5.06020071944e-09
cal_temp 5.06020071944e-09
G2err 3.57807540732e-09
R4Aerr 3.57790867802e-09
G1err 3.57781283015e-09
C1stab 4.17223154494e-10
C1alpha 8.34446308987e-11
C1err 4.17164751428e-11

Calculate values for the dials set to full scale
1Z ?0.0001000198? H ?1.000071? ohm
2Z ?0.001000198? H ?10.00071? ohm
3Z ?0.01000198? H ?100.0071? ohm
4Z ?0.1000198? H ?1000.071? ohm
5Z ?1.000198? H ?10000.71? ohm
6Z ?10.00134? H ?100000.6? ohm
7Z ?100.0211? H ?1000084? ohm
1Y ?0.0001000089? F ?0.999962? siemen
2Y ?0.00001000089? F ?0.0999962? siemen
3Y ?0.000001000089? F ?0.00999962? siemen
4Y ?0.0000001000089? F ?0.000999962? siemen
5Y ?0.00000001000089? F ?0.0000999962? siemen
6Y ?0.000000001000153? F ?0.00001000026? siemen
7Y ?0.0000000001000076? F ?0.000000999949? siemen

Compare with the UB's nameplate range values as, embedded in available_ranges
1Z 0.0001 H 	 1.0 ohm
2Z 0.001 H 	 10.0 ohm
3Z 0.01 H 	 100.0 ohm
4Z 0.1 H 	 1000.0 ohm
5Z 1.0 H 	 10000.0 ohm
6Z 10.0 H 	 100000.0 ohm
7Z 100.0 H 	 1000000.0 ohm
1Y 0.0001 F 	 1.0 siemen
2Y 1e-05 F 	 0.1 siemen
3Y 1e-06 F 	 0.01 siemen
4Y 1e-07 F 	 0.001 siemen
5Y 1e-08 F 	 0.0001 siemen
6Y 1e-09 F 	 1e-05 siemen
7Y 1e-10 F 	 1e-06 siemen

Check standards have been correctly constructed, including uncertainties
ureal(10000.541000000001,0.030823693609093962,21.7)
R4Aacdc 0.0250013525
R4Aalpha 0.010000541
R4Astab 0.010000541
temperature 0.00707145035663
cal_temp 0.00707145035663
R4Aerr 0.005

ureal(99999.00999999999,0.3082177290635078,21.7)
R4Bacdc 0.249997525
R4Balpha 0.09999901
R4Bstab 0.09999901
temperature 0.0707099780829
cal_temp 0.0707099780829
R4Berr 0.05

ureal(1000067.2,3.08240867533355,21.7)
R4Cacdc 2.500168
R4Calpha 1.0000672
R4Cstab 1.0000672
temperature 0.707154298762
cal_temp 0.707154298762
R4Cerr 0.5

ureal(1.0000129000000000355501302e-05,3.0822457156332133711412704e-11,21.7)
G1acdc 2.50003225e-11
G1alpha 1.0000129e-11
G1stab 1.0000129e-11
temperature 7.07115902864e-12
cal_temp 7.07115902864e-12
G1err 5e-12

ureal(1.0000075000000000982483886e-05,3.0822295097088918821350794e-11,21.7)
G2acdc 2.50001875e-11
G2alpha 1.0000075e-11
G2stab 1.0000075e-11
temperature 7.07112084487e-12
cal_temp 7.07112084487e-12
G2err 5e-12

ureal(1.00014000000000002250705931985e-09,5.2208773715918670982318844249e-15,11.9)
C1stab 5.0007e-15
C1alpha 1.00014e-15
temperature 7.07205776136e-16
cal_temp 7.07205776136e-16
C1err 5e-16

ureal(1.0000109999999999832,1.399999999999999979e-05,10)
Amperr 1.4e-05

Compare calculated uncertainty and CMCs for different dial settings
0 0 0.372429464121 37.2429464121 37.2382134962
0 100 0.37242946412 37.2429464251 37.2381699451
0 1000 0.372429464057 37.2429477139 37.2338591447
0 1000000 0.372365670836 38.5223535901 2.59076353297
1000 0 0.374304373019 37.2429060754 37.2381616548
1000 100 0.374304373019 37.2429060884 37.2381181039
1000 1000 0.374304372955 37.2429073773 37.2338073215
1000 1000000 0.374240258721 38.5223118803 2.59076351842
10000 0 0.391140221797 37.2389133948 37.2330304224
10000 100 0.391140221796 37.2389134078 37.2329868895
10000 1000 0.39114022173 37.2389146964 37.2286778897
10000 1000000 0.391073238069 38.5181832821 2.59076207805
10000000 0 1.29608562515 2.69196585731 2.30428235204
10000000 100 1.29608562515 2.69196585814 2.30428234388
10000000 1000 1.29608562515 2.69196594043 2.30428153552
10000000 1000000 1.29608460535 2.77383574728 1.806010243
caldata ['C1', 'R4A', 'R4B', 'R4C', 'G1', 'G2']
ureal(0.0001399999999999999877355,4.999999999999999773741e-07,10,label='C1err')
ureal(5.410000000000000041061e-05,4.999999999999999773741e-07,10,label='R4Aerr')
ureal(-9.900000000000000081402e-06,4.999999999999999773741e-07,10,label='R4Berr')
ureal(6.719999999999999357094e-05,4.999999999999999773741e-07,10,label='R4Cerr')
ureal(1.29000000000000001574e-05,4.999999999999999773741e-07,10,label='G1err')
ureal(7.500000000000000190006e-06,4.999999999999999773741e-07,10,label='G2err')
ureal(1.0999999999999999714e-05,1.399999999999999979e-05,10,label='Amperr')
ureal(0,1.6999999999999999866e-05,10,label='Amp_phase_err')
ureal(20,0.5,10,label='cal_temp')
acdc ['C1', 'R4A', 'R4B', 'R4C', 'G1', 'G2']
ureal(0,0,inf,label='C1acdc')
ureal(0,2.50000000000000062802e-06,10,label='R4Aacdc')
ureal(0,2.50000000000000062802e-06,10,label='R4Bacdc')
ureal(0,2.50000000000000062802e-06,10,label='R4Cacdc')
ureal(0,2.50000000000000062802e-06,10,label='G1acdc')
ureal(0,2.50000000000000062802e-06,10,label='G2acdc')
tempcoeffs ['C1', 'R4A', 'R4B', 'R4C', 'G1', 'G2']
ureal(0,1.9999999999999999095e-06,10,label='C1alpha')
ureal(0,1.9999999999999999095e-06,10,label='R4Aalpha')
ureal(0,1.9999999999999999095e-06,10,label='R4Balpha')
ureal(0,1.9999999999999999095e-06,10,label='R4Calpha')
ureal(0,1.9999999999999999095e-06,10,label='G1alpha')
ureal(0,1.9999999999999999095e-06,10,label='G2alpha')
stability ['C1', 'R4A', 'R4B', 'R4C', 'G1', 'G2']
ureal(0,5.00000000000000040902e-06,10,label='C1stab')
ureal(0,9.99999999999999954748e-07,10,label='R4Astab')
ureal(0,9.99999999999999954748e-07,10,label='R4Bstab')
ureal(0,9.99999999999999954748e-07,10,label='R4Cstab')
ureal(0,9.99999999999999954748e-07,10,label='G1stab')
ureal(0,9.99999999999999954748e-07,10,label='G2stab')
['Adial_lin', 'Bdial_lin']
ureal(0,9.99999999999999954748e-07,10,label='Adial_lin')
ureal(0,9.99999999999999954748e-07,10,label='Bdial_lin')

Process finished with exit code 0
"""
