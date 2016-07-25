from astropy.io.ascii import basic, core
from astropy.table import Table
import numpy as np


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