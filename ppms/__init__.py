from astropy.io.ascii import basic, core
from astropy.table import Table
from astropy import units as u, constants as c
import numpy as np
import dateutil.parser as dparser


def normalize(table):
    data, names, mask = [], [], []
    for col in table.columns.values():
        if not col.mask.all():
        #if col.name == 'Comment':
        #    print(col)
        #    err = ValueError()
        #    err.col = col
        #    raise err
        #if col.any():
            data.append(col.quantity.si if col.unit else col)
            names.append(col.name)
            mask.append(col.mask)
    table = Table(data=data, masked=True, names=names, meta=table.meta)
    table.mask = np.array(mask)
    return table


class PPMSHeader(basic.CsvHeader):
    UNITS = {
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

def acms(path, volume=None, formula_units=None, mode='acdc'):
    """Reduce and preprocess acms dataset.
    
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
        data['B'] = source[ac_mask]['Magnetic Field'].round(4)
        data['T'] = source[ac_mask]['Temperature']
    else:
        data['B'] = source[dc_mask]['Magnetic Field'].round(4)
        data['T'] = source[dc_mask]['Temperature']
    
    if (mode == 'ac') or (mode == 'acdc'):
        data['B_ac'] = source[ac_mask]["Amplitude"].quantity
        data['m_ac'] = source[ac_mask]["M'"].quantity
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