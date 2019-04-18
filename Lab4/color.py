import numpy as np

def color(alt, az):
    """
    Inputs altitude and azimuth arrays and determines which points lie in range for Leuschner.
    Arguments:
    alt, altitude (array-like)
    az, azimuth (array-like)
    Returns:
    col, color code for plotting (array-like)
    ind_bad, indicies of out of range points (array-like)
    ind_good, indicies of in range points (array-like)
    """
    flalt = alt.flatten()
    flaz = az.flatten()

    min_alt, max_alt = 15, 85
    min_azi, max_azi = 5, 350
    
    col = ['g']*len(flalt)
    col = np.array(col)
    col[np.where(flalt<min_alt)] = 'r'
    col[np.where(flaz<min_azi)] = 'r'
    col[np.where(flalt>max_alt)] = 'r'
    col[np.where(flaz>max_azi)] = 'r'
    
    ind_bad = np.where(col == 'r')
    ind_good = np.where(col == 'g')
   
    return col, ind_bad, ind_good
