from astropy.io.ascii import basic, core
from astropy.table import Table, MaskedColumn
from astropy import units as u, constants as c
import numpy as np
import dateutil.parser as dparser
from scipy.ndimage import median_filter


class MaglabHeader(basic.CsvHeader):
    comment = r'\s*;'
    write_comment = ';'
    start_line = 1
        
    def get_cols(self, lines):
        lines = self.process_lines(lines)
        start_line = self.start_line
        for i, line in enumerate(lines):
            if i == start_line:
                break
        else:  # No header line matching
            raise ValueError('No header line found in table')
            
        self.names = [x.strip() for x in next(self.splitter([line]))]
        self.units = next(self.splitter([next(lines)]))
        self._set_cols_from_names()
        for c, u in zip(self.cols, self.units):
            setattr(c, 'unit', u)


class MaglabData(core.BaseData):
    comment = r'\s*;'
    write_comment = ';'
    start_line = 5


class Maglab(basic.Csv):
    """Reads a Oxford Instruments Maglab data file."""
    _format_name = 'maglab'
    _io_registry_can_write = False
    _description = 'Oxford Instruments Maglab data file reader'
    header_class = MaglabHeader
    data_class = MaglabData


def normalize(table):
	data = []
	for col in table.columns.values():
		if isinstance(col, MaskedColumn) and col.mask.all():
			# drop completely empty columns
			continue
		data.append(col)
	return Table(data)


class PPMSHeader(basic.CsvHeader):
    UNITS = {
        # Heat Capacity Option units
        'Seconds': 'second',
        'seconds': 'second',
        'Oersted': '0.0001 * T',
        'Kelvin': 'K',
        'µJ/K': 'uJ/K',
        'µJ/K/K': 'uJ/K/K',
        'Seconds': 'second',
        # ACMS option units
        'sec': 'second',
        'emu': 'erg/gauss',
        'Oe': '0.0001 * T',
        'code': None,
    }
    comment = r'\s*;'
    write_comment = ';'

    def start_line(self, lines):
        return list(lines).index('[Data]') + 1
    
    def _set_cols_from_names(self):
        names, units = [], []
        for header in self.names:
            if '(' in header:
                h, u = [x.strip() for x in header.replace(')', '').split('(')]
            else:
                h, u = header.strip(), None
            names.append(h)
            
            if u in self.UNITS:
                units.append(self.UNITS[u])
            else:
                units.append(u)
                
        self.names = names
        super(PPMSHeader, self)._set_cols_from_names()
        for col, unit in zip(self.cols, units):
            if unit:
                col.unit = unit


class PPMSData(basic.CsvData):
    comment = r'\s*;'
    write_comment = ';'
    def start_line(self, lines):
        return list(lines).index('[Data]') + 2


class PPMSOutputter(core.TableOutputter):
    def __call__(self, cols, meta):
        cols = [c for c in cols if any(c.str_vals)]
        return normalize(super(PPMSOutputter, self).__call__(cols, meta))


class PPMS(basic.Csv):
    """Reads a Quantum Design PPMS data file."""
    _format_name = 'ppms'
    _io_registry_can_write = False
    _description = 'Quantum Design PPMS data file reader'
    header_class = PPMSHeader
    data_class = PPMSData
    outputter_class = PPMSOutputter
    
    
fu = u.def_unit(['f.u.', 'formula unit'], u.dimensionless_unscaled)

def acms_legacy(path, volume=None, formula_units=None, mode='acdc'):
    """Reduce and preprocess acms dataset.
    
    ..note::
        The colnames were changed. This function still uses the old
        schema.
    
    :param volume: The sample volume.
    :param formula_units: The numer of formula units of the sample.
    :param mode: Data modes, either 'acdc', 'ac' or 'dc'.
    
    """
    if volume:
        if not isinstance(volume, u.Quantity):
            raise ValueError('Missing type of volume parameter.')
    
    source = Table.read(path, format='ascii.ppms')
    # Boolean mask, True for DC magnetisation measurements
    dc_mask = source['Measure Type'] == 0
    ac_mask = source['Measure Type'] == 5
    
    if mode == 'acdc' and (np.sum(dc_mask) != np.sum(~ac_mask)):
        raise ValueError('Nonequal number of ac ({}) and dc ({}) measurements'.format(np.sum(ac_mask), np.sum(dc_mask)) )
    
    data = Table(masked=False)
    if mode == 'ac':
        data['B'] = source[ac_mask]['Magnetic Field'].to(u.T).round(4)
        data['T'] = source[ac_mask]['Temperature']
    else:
        data['B'] = source[dc_mask]['Magnetic Field'].to(u.T).round(4)
        data['T'] = source[dc_mask]['Temperature']
    
    if (mode == 'ac') or (mode == 'acdc'):
        data['B_ac'] = source[ac_mask]["Amplitude"].to(u.T)
        data['m_ac'] = source[ac_mask]["M'"]
        data["m'_ac"] = source[ac_mask]["M''"]
        if volume:
            H = data['B_ac'].quantity / c.mu0
            M = data['m_ac'].quantity / volume
            M_imag = data["m'_ac"].quantity / volume
            data["χ"] = (M / H).si
            data["χ'"] = (M_imag / H).si

    if (mode == 'dc') or (mode == 'acdc'):
        data['m'] = source[dc_mask]['M-DC']

        if formula_units:
            # calculate magnetic moment per formula unit
            data['m_fu'] = data['m'].to(c.muB) / formula_units
        if volume:
            # calculate magnetisation.
            data['M'] = (data['m'].quantity / volume).si
    
    data.meta['temperature'] = np.round(data['T'].mean(), 1)
    if volume:
        data.meta['volume'] = volume
    data.meta['z'] = source['Sample Center'].quantity[0].value
    if mode == 'ac' or mode == 'acdc':
        data.meta['frequency'] = source[ac_mask]['Frequency'][0]
    data.meta['path'] = path
    try:
        # Try to extract date information from filepath
        data.meta['date'] = dparser.parse(path,fuzzy=True)
    except ValueError:
        pass
    return data


def acms(path, volume=None, formula_units=None, demag=None, mode='acdc', scan=None, masked=False):
    """Reduce and preprocess acms dataset.
    
    :param volume: The sample volume.
    :param formula_units: The numer of formula units of the sample.
    :param demag: The demagnetizing factor. To calculate the demagnetizing correction,
        the magnetisation M is needed. Therefore the volume is mandatory and it only
        works with modes 'dc' or 'acdc'.
    :param mode: Data modes, either 'acdc', 'ac' or 'dc'.
    :param scan: The scan variable. if scan is 'B' then dM/dH can be calculated.
    
    """
    if mode not in {'ac', 'dc', 'acdc'}:
        raise ValueError("invalid mode. Must be one of 'ac', 'dc' or 'acdc'")
    if volume:
        if not isinstance(volume, u.Quantity):
            raise ValueError('Missing type of volume parameter.')
    if demag:
        if volume is None:
            raise ValueError(
                'volume parameter is neccessary to calculate the'
                'magnetisation used for demagnetizing correction.')
        if mode == 'ac':
            raise ValueError(
                "Can't calculate demagnetizing correction with mode"
                "'ac'. Magnetisation is neccessary.")
    source = Table.read(path, format='ascii.ppms')
    if masked:
        data = Table()
        dc_mask = ac_mask = np.ones(len(source), dtype=bool)
    else:
        data = Table(masked=False)
        # Boolean mask, True for DC magnetisation measurements
        dc_mask = source['Measure Type'] == 0
        ac_mask = source['Measure Type'] == 5

        if mode == 'acdc' and (np.sum(dc_mask) != np.sum(ac_mask)):
            raise ValueError('Nonequal number of ac ({}) and dc ({}) measurements'.format(np.sum(ac_mask), np.sum(dc_mask)) )
    
    if mode == 'ac':
        data['B'] = source[ac_mask]['Magnetic Field'].to(u.T).round(4)
        data['T'] = source[ac_mask]['Temperature']
    else:
        data['B'] = source[dc_mask]['Magnetic Field'].to(u.T).round(4)
        data['T'] = source[dc_mask]['Temperature']
    data['H'] = data['B'] / c.mu0
    
    if (mode == 'ac') or (mode == 'acdc'):
        data['t_ac'] = source[ac_mask]['Time Stamp']
        data['f'] = source[ac_mask]['Frequency']
        data['B_ac'] = source[ac_mask]["Amplitude"].to(u.T)
        data["m'_ac"] = source[ac_mask]["M'"]
        data["m''_ac"] = source[ac_mask]["M''"]
        if volume:
            H = data['B_ac'].quantity / c.mu0
            M = data["m'_ac"].quantity / volume
            M_imag = data["m''_ac"].quantity / volume
            data["χ'"] = (M / H).si
            data["χ''"] = (M_imag / H).si
        
        # Handle higher harmonic susceptibilities
        harmonics_real = [x for x in source[ac_mask].columns if x.startswith("M '[")]
        for colname in harmonics_real:
            i = int(colname[4:-1])
            data["m'_ac[{}]".format(i)] = source[ac_mask][colname]
            if volume:
                M_i = data["m'_ac[{}]".format(i)].quantity / volume
                data["χ'[{}]".format(i)] = (M_i / H).si

        harmonics_imag = [x for x in source[ac_mask].columns if x.startswith("M ''[")]
        for colname in harmonics_imag:
            i = int(colname[5:-1])
            data["m''_ac[{}]".format(i)] = source[ac_mask][colname]
            if volume:
                M_imag_i = data["m''_ac[{}]".format(i)].quantity / volume
                data["χ''[{}]".format(i)] = (M_imag_i / H).si        

    if (mode == 'dc') or (mode == 'acdc'):
        data['t_dc'] = source[dc_mask]['Time Stamp']
        data['m'] = source[dc_mask]['M-DC']

        if formula_units:
            # calculate magnetic moment per formula unit
            data['m_fu'] = data['m'].to(c.muB) / formula_units
        if volume:
            # calculate magnetisation.
            data['M'] = M = (data['m'].quantity / volume).si
            if scan == 'B':
                data['dM/dH'] = np.gradient(M) / np.gradient(data['H'])
    
    if demag:
        demagnetizing_correction(data, demag=demag)
    
    data.meta['temperature'] = np.round(data['T'].mean(), 1)
    if volume:
        data.meta['volume'] = volume
    data.meta['z'] = source['Sample Center'].quantity[0].value
    data.meta['path'] = path
    try:
        # Try to extract date information from filepath
        data.meta['date'] = dparser.parse(path, fuzzy=True)
    except ValueError:
        pass
    return data


def demagnetizing_correction(data, demag):
    """Calculates the demagnetizing correction.
    
    The ac susceptibility is corrected according to [1]
    
    [1]: Youssif, M. I., Bahgat, A. A. & Ali, I. A. AC magnetic
         susceptibility technique for the characterization of high
         temperature superconductors. Egyptian Journal of Solids 23,
         231–250 (2000).

    """
    Hext, M = data['H'], data['M']
    data['H_int'] = Hint = Hext - demag * M
    data['B_int'] = c.mu0 * Hint
    
    #scale = Hext / Hint
    #idx = Hext == 0
    #scale[idx] = median_filter(scale, size=5)[idx]
    
    for name, col in data.columns.items():
        if name == 'dM/dH':
            data['dM/dH_int'] = np.gradient(M) / np.gradient(data['H_int'])
        elif name == "χ'" or name.startswith("χ'["):
            chi_r = col
            chi_i = data[name.replace("'", "''")]
            
            data[name + '_int'] = (chi_r - demag * (chi_r**2 + chi_i**2)) / (demag**2 * (chi_r**2 + chi_i**2) - 2 * demag * chi_r + 1)

        elif name == "χ''" or name.startswith("χ''["):
            chi_i = col
            chi_r = data[name.replace("''", "'")]
            
            data[name + '_int'] = chi_i / (demag**2 * (chi_r**2 + chi_i**2) - 2 * demag * chi_r + 1)
        #data[name + '_int'] = col * scale


def magnetic_moment_in_fu(m, formula_units):
    """Converts the magnetic moment from si units to Bohr magneton per formula
    units.
    
    :param m: Magnetic moment.
    :param formula_units: The number of formula units.
    
    """
    return m.to(c.muB) / formula_units


def heatcapacity(path):
    # The HC option sometimes creates comment lines without commas.
    with open(path, 'r', encoding='cp1252', newline='') as f:
        buffer = ''.join([l for l in f.readlines() if 'Error' not in l])
    source = Table.read(buffer, format='ascii.ppms')
    data = Table(masked=False)
    data['B'] = source['Field'].to(u.T).round(4)
    data['T'] = source['System Temp']
    data['Tsample'] = source['Sample Temp']
    data['Tpuck'] = source['Puck Temp']
    data['C'] = source['Total HC']
    data['Csample'] = source['Samp HC']
    data['Caddenda'] = source['Addenda HC']
    data['Coffset'] = source['Addenda Offset HC']
    return data