import requests

import pandas as pd
from bs4 import BeautifulSoup

def _remove_multiple_values(line, value=""):
    while (value in line):
        line.remove(value)
    return line


def sst_columns():
    return ['year', 'month', 'nino12', 'anom12', 'nino3', 'anom3', 'nino4', 'anom4', 'nino34', 'anom34']


def oni_columns():
    return ["seas", "year", "total", "anom"]


def _oni_values_to_row(values):
    seas, year, total, anom = values
    d = dict(
        seas=seas,
        year=int(year),
        total=float(total),
        anom=float(anom)
    )
    return d


def _sst_values_to_row(values):
    year, month, nino12, anom12, nino3, anom3, nino4, anom4, nino34, anom34 = values
    d = dict(
        year=int(year),
        month=int(month),
        nino12=float(nino12),
        anom12=float(anom12),
        nino3=float(nino3),
        anom3=float(anom3),
        nino4=float(nino4),
        anom4=float(anom4),
        nino34=float(nino34),
        anom34=float(anom34)
    )
    return d


smn_columns = dict(
    sst=sst_columns(),
    oni=oni_columns()
)

smn_separated = dict(
    sst="  ",
    oni=" ",
)

smn_values_to_row = dict(
    sst=_sst_values_to_row,
    oni=_oni_values_to_row,
)


def _rawdata_to_df(raw, var=None):
    if not var:
        print("")
        return None
    lines = raw.split("\n")
    df = pd.DataFrame(columns=smn_columns[var])
    for i in range(1, len(lines) - 1):
        values = _remove_multiple_values(lines[i].split(smn_separated[var]))
        df = df._append(smn_values_to_row[var](values), ignore_index=True).copy()
    return df 


def get_smn_data(url, var=None):
    response = requests.get(url)
    soup = BeautifulSoup(response.content, "html.parser")
    raw = soup.text
    if var:
        return _rawdata_to_df(raw, var=var)
    else:
        return raw
