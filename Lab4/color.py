import numpy as np

def color(alt, az):

    """
    Inputs altitude and azimuth arrays and determines which points lie in range for Leuschner.
    Arguments:
    alt, altitude (array)
    az, azimuth (array)
    Returns:
    col, color code for plotting (array-like)
    ind, indicies of out of range points (array-like)
    """

    flalt = alt.flatten()
    flaz = az.flatten()

    min_alt, max_alt = 15, 85  #Obtained from ugradio.leusch
    min_az, max_az = 5, 350    #Obtained from ugradio.leusch
    
    col = ['b']*len(flalt)
    col = np.array(col)
    col[np.where(flalt<min_alt)] = 'r'
    col[np.where(flaz<min_az)] = 'r'
    col[np.where(flalt>max_alt)] = 'r'
    col[np.where(flaz>max_az)] = 'r'
    
    ind = np.where(col == 'r')
   
    return col, ind
