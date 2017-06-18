__all__ = ["unix2mjd", "mjd2unix", "mjd2djd"]

def unix2mjd(timestamp):
    """Convert from unix timestamp to MJD

    Parameters
    ----------
    timestamp : float
        The unix timestamp to be converted.

    Returns
    -------
    The modified Julian date corresponding to `timestamp`

    Notes
    -----
    This does the crudest possible conversion, and does not take into
    account leap seconds (or potentially other discrepancies I don't know about)
    """
    return timestamp / 86400 + 40587

def mjd2unix(mjd):
    """Convert from MJD to unix timestamp

    Parameters
    ----------
    mjd : float
        The modified Julian date to be converted.

    Returns
    -------
    The unix timestamp corresponding to `mjd`.

    Notes
    -----
    This does the crudest possible conversion, and does not take into
    account leap seconds (or potentially other discrepancies I don't know about)
    """
    return (mjd - 40587) * 86400

def mjd2djd(mjd):
    """Convert from modified Julian date to Dublin Julian date

    pyephem uses Dublin Julian dates

    Parameters
    ----------
    mjd : float
        The modified Julian date to be converted.

    Returns
    -------
    The Dublin Julian date corresponding to `mjd`.
    """
    # (this function adapted from Peter Yoachim's code)
    doff = 15019.5 # this equals ephem.Date(0)-ephem.Date('1858/11/17')
    return mjd - doff

def djd2mjd(djd):
    """Convert from Dublin Julian date to Modified Julian date

    pyephem uses Dublin Julian dates

    Parameters
    ----------
    djd : float
        The Dublin Julian date to be converted.

    Returns
    -------
    The modified Julian date corresponding to `djd`.
    """
    doff = 15019.5 # this equals ephem.Date(0)-ephem.Date('1858/11/17')
    return djd + doff
