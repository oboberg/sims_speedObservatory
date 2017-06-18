__all__ = ["unix2mjd", "mjd2unix", "mjd2djd"]

def unix2mjd(timestamp):
    return timestamp / 86400 + 40587

def mjd2unix(mjd):
    return (mjd - 40587) * 86400

def mjd2djd(mjd):
    """
    Convert Modified Julian Date to Dublin Julian Date (what pyephem uses).
    (this function adapted from Peter Yoachim's code)
    """
    doff = 15019.5 # this equals ephem.Date(0)-ephem.Date('1858/11/17')
    return mjd - doff
