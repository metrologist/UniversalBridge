from __future__ import print_function
from __future__ import division

"""
calculate_imp.py calculates immittances of components measured with the Universal Impedance Bridge
"""
import ubridge as ub
import GTC as gtc
import math
import excel


class UUT(object):
    """
    A set of immittances that will be calibrated using the Universal Bridge. Bridge readings, environmental conditions
    and sensitivity coefficients are gathered in a spreadsheet. Calculated calibration results will be returned also in
    a spreadsheet for preparing a full report.
    """

    def __init__(self, bridge, input_workbook, target_in_sheet, block_descriptor, output_workbook, target_out_sheet):
        self.bridge = bridge
        self.input_workbook = input_workbook
        self.block_descriptor = block_descriptor  # this simply has to correctly match the spreadsheet.[9, 35, 1, 15]
        self.output_workbook = output_workbook
        self.target_in_sheet = target_in_sheet
        self.target_out_sheet = target_out_sheet
        self.data_block, self.datdict, self.excel_obj = self.get_measurement_data()

    def get_measurement_data(self):
        #  get component measurement data from a spread sheet
        this_uut = excel.CALCULATOR(self.input_workbook, self.output_workbook)
        my_copy_data = this_uut.getdata_block(self.target_in_sheet, self.block_descriptor)
        thisdict = []  # need a list of dictionaries
        for x in my_copy_data:
            #  put each 'column' value into a dictionary for clarity
            thisdict.append(
                {'item': x[0], 'nom_freq': x[1], 'ubrange': x[2], 'xdial': x[3], 'rdial': x[4], 'frequency': x[5],
                 'temperature': x[6], 'tempu': x[7], 'tempdf': x[8], 'uuttempcox': x[9], 'uuttempcoxu': x[10],
                 'uuttempcoxdf': x[11], 'uuttempcor': x[12], 'uuttempcoru': x[13], 'uuttempcordf': x[14]})
        return my_copy_data, thisdict, this_uut

    def calculate_values(self, room_temperature):
        y = []
        yzero = []
        for x in self.datdict:
            answer = ubridge.bridge_value(x['ubrange'], x['rdial'], x['xdial'], x['frequency'], 1)
            answer._set_label(x['item'] + str(x['nom_freq']) + x['ubrange'])  # str so e.g. 1k and 1000 both are strings
            #  next multiply the ubrdige values by factors relating to the uut
            if x['item'] != 'coax zero':  # no uncertainty in the coax zero uut
                uuttemp = gtc.ureal(x['temperature'], x['tempu'], x['tempdf'],
                                    'uuttemp' + x['item'])
                rtempco = gtc.ureal(x['uuttempcor'], x['uuttempcoru'], x['uuttempcordf'],
                                    'rtempco' + x['item'])
                r_value = answer.real * (1 + gtc.function.mul2(rtempco, uuttemp - room_temperature))
                xtempco = gtc.ureal(x['uuttempcox'], x['uuttempcoxu'], x['uuttempcoxdf'],
                                    'xtempco' + x['item'])
                x_value = answer.imag / (2 * math.pi * x['frequency']) * (
                    1 + gtc.function.mul2(xtempco, uuttemp - room_temperature))
            else:
                assert x['item'] == 'coax zero', "label not 'coax zero' when it should be"
                r_value = answer.real
                x_value = answer.imag / (2 * math.pi * x['frequency'])
            # create lists of the ubridge values
            y.append(
                (x['item'], x['ubrange'], str(x['nom_freq']), r_value,
                 x_value))  # all tuples of R-L or G-C, including coax zeros
            if x[
                'item'] == 'coax zero':  # a list of just the coax zeros is also created, duplicating the entry in the full list
                yzero.append((x['ubrange'], str(x['nom_freq']), r_value, x_value))
        return y, yzero

    def subtract_zeros(self, all_items, zeros):
        y_corrected = []  # zero corrected values
        for x in all_items:  # all_items is full list of answers
            for zero in zeros:  # zeros is the list of coax zeros
                if zero[1] == x[2]:  # find a zero of the same nominal frequency
                    if zero[0] == x[1]:  # find if the zero is also the correct range
                        zero_corrected = (x[3] - zero[2], x[4] - zero[3])  # correct r and x values
                        if x[
                            0] == 'coax zero':  # check that the zeros cancel, just in case there is a duplicate coax zero
                            assert zero_corrected[0].x == 0, 'zero is not cancelling when subtracted from itself'
                            assert zero_corrected[1].x == 0, 'zero is not cancelling when subtracted from itself'
                        y_corrected.append(zero_corrected)
        # if no matching coax zero was found then y_corrected will be missing the item
        assert len(y_corrected) == len(all_items), "value missing a matching coax zero"
        return y_corrected

    def cmc_check(self):
        cmcs = []  # for the list of cmcs
        for x in self.datdict:
            cmc = ubridge.cmc_uncert(x['ubrange'], x['rdial'], x['xdial'], x['frequency'])
            cmcs.append(cmc)
        return cmcs

    def cmc_replace(self, cmcs, calculated):
        pass

    def create_output(self, y_corrected, cmcs):
        # now put all info in spreadsheet
        for i in range(len(y_corrected)):
            rpart = y_corrected[i][0]
            xpart = y_corrected[i][1]
            # first calculate coverage factors
            if math.isnan(rpart.df):  # avoiding the coax zeros with no uncertainty
                kr = 0
            else:
                kr = gtc.reporting.k_factor(rpart.df)
            if math.isnan(xpart.df):
                kx = 0
            else:
                kx = gtc.reporting.k_factor(xpart.df)
            # then expanded unertainties
            xU = xpart.u * kx
            rU = rpart.u * kr
            self.data_block[i].append(xpart.x)
            self.data_block[i].append(xU)  # expanded uncertainty
            self.data_block[i].append(kx)
            self.data_block[i].append(rpart.x)
            self.data_block[i].append(rU)  # expanded uncertainty
            self.data_block[i].append(kr)
            # then compare with CMCs
            # first decide on which r CMC to use
            # note this test only applies to Z ranges
            self.data_block[i].append(cmcs[i][1])  # for now, add the CMCs in a column to test in excel
            self.data_block[i].append(cmcs[i][0])
            # put the block, with final answers added, into the output spreadsheet
        self.excel_obj.makeworkbook(self.data_block, self.target_out_sheet)
        return

    def inspect_budget(self, line_no, results):
        # Can also look at the budget for individual lines
        print('Budget for', self.data_block[line_no][0], self.data_block[line_no][1])
        item = results[line_no]
        print(repr(item[0]))
        for l, u in gtc.reporting.budget(item[0], trim=0):
            print(l, '\t\t\t', u)
        print()
        print(repr(item[1]))
        for l, u in gtc.reporting.budget(item[1], trim=0):
            print(l, '\t\t\t', u)
        print()
        return

    def muirhead_zeros(self, y_corrected):
        """
        Returns list of muirhead capacitor values with jig zero subtracted
        Includes a 0.05 pF standard uncertainty in jig definition
        y_corrected: capacitance values that have had coaxial zeros subtracted
        """
        jig_zero = gtc.ureal(0, 0.05e-12, 10, 'jig_zero')  # 0.05 pF 1 sigma
        names = ['M1', 'M0.1', 'M0.01', 'M0.001', 'M0.0005']
        muirhead_caps = []
        for i in range(len(y_corrected)):
            if self.data_block[i][0] in names:
                for j in range(len(y_corrected)):
                    if self.data_block[j][0] == self.data_block[i][0] + 'zero':  # e.g. 'M1zero, 'M0.0005zero'
                        print(repr(y_corrected[i][0]))
                        thing = jig_zero + y_corrected[i][0]
                        jig_corrected = (
                        y_corrected[i][0] - y_corrected[j][0], jig_zero + y_corrected[i][1] - y_corrected[j][1])
                        muirhead_caps.append(jig_corrected)
        return muirhead_caps

    def muirhead_tand(self, caps):
        """
        assumes 10,000 rad/s and returns G/wC in mrad
        also assumes a 0.005 ohm standar uncertainty in series connection
        caps:  is the list of jig_corrected muirhead capacitor values
        """
        w = 1e4
        seriesR = 0.005
        tan_deltas = []
        for i in range(len(caps)):
            seriesG = gtc.ureal(0, (w * caps[i][1].x) ** 2 * seriesR, 10, 'seriesG')
            calc_tand = gtc.atan((caps[i][0] + seriesG) / (caps[i][1] * w))
            tan_deltas.append(calc_tand * 1e3)
        return tan_deltas

    def muirhead_output(self, caps, tandeltas):
        """
        Creates spreadsheet with final report values, in F and mrad
        caps: is jig corrected capacitances from muirheadzeros
        tandeltas:  is tan delta values from muirhead_tand
        """
        # now put all info in spreadsheet
        data_block = [[], [], [], [], []]
        for i in range(len(caps)):
            rpart = caps[i][0]
            xpart = caps[i][1]
            # first calculate coverage factors
            kr = gtc.reporting.k_factor(rpart.df)
            kx = gtc.reporting.k_factor(xpart.df)
            # then expanded unertainties
            xU = xpart.u * kx
            rU = rpart.u * kr
            data_block[i].append(xpart.x)
            data_block[i].append(xU)  # expanded uncertainty
            data_block[i].append(kx)
            data_block[i].append(rpart.x)
            data_block[i].append(rU)  # expanded uncertainty
            data_block[i].append(kr)
            # now for loss angle
            td = tandeltas[i].x
            utd = tandeltas[i].u
            ktd = gtc.reporting.k_factor(tandeltas[i].df)
            Utd = ktd * utd
            data_block[i].append(td)
            data_block[i].append(Utd)
            data_block[i].append(ktd)
        self.excel_obj.makeworkbook(data_block, 'muirhead')


if __name__ == "__main__":
    calfile = r'ubdict_nov_2017.csv'
    room_temperature = gtc.ureal(20, 0.5, 10,
                                 'temperature')  # this should be the ambient temperature given in conditions
    # create the bridge object
    ubridge = ub.UNIVERSALBRIDGE(calfile,
                                 room_temperature)  # temperature put here as it might be in common with the UUT
    block_descriptor = [9, 24, 1, 15]  # this simply has to correctly match the spreadsheet.[9, 35, 1, 15]
    cap_set = UUT(ubridge, 'S22012.xlsx', 'pyUBreadings', block_descriptor, 'capResults.xlsx', 'pyUBresults')
    # note that a different temperature could be used for the UUT
    cap_answers, cap_zeros = cap_set.calculate_values(room_temperature)
    cap_zero_corrected = cap_set.subtract_zeros(cap_answers, cap_zeros)
    cmc_list = cap_set.cmc_check()
    ##     cap_set.create_output(cap_zero_corrected, cmc_list)  #creates capResults.xlsx with just C, G results
    a = cap_set.muirhead_zeros(cap_zero_corrected)
    for i in range(len(a)):
        print(a[i])
    b = cap_set.muirhead_tand(a)
    for i in range(len(b)):
        print(repr(b[i]))
    cap_set.muirhead_output(a, b)  # creates capResults.xlsx with tan delta done in GTC

    # Can also print the budget for individual lines
##     cap_set.inspect_budget(4, cap_zero_corrected)
##     cap_set.inspect_budget(10, cap_zero_corrected)

