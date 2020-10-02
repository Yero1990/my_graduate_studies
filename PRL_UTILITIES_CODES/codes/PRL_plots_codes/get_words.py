import numpy as np
import sys

def get_words(w=0, h=0, ncol=0):
    #PRL counts figure sizes as words, for which there is a conversion factor.
    
    fig_w = float(w)  #figure width (inches)
    fig_h = float(h)  #figure height (inches)
    aspect_ratio = fig_w / fig_h;

    if(ncol==1):
        #For single column figures, the word estimate for a figure size is:
        one_column_fig_words = (150./aspect_ratio) + 20.  #units: number of words
        print('(single column) # of words equivalent for figsize: (w=',w,' h=', h,') in. -----> ',one_column_fig_words,' words')
    elif(ncol==2):
        #For figures that occupy both columns of the paper, the word estimate for a figsize is:
        two_column_fig_words = (300./(0.5*aspect_ratio) ) + 40. 
        print('(two columns) # of words equivalent for figsize: (w=',w,' h=', h,') in. -----> ',two_column_fig_words,' words')


if __name__ == '__main__':
    w = float(sys.argv[1])
    h = float(sys.argv[2])
    ncol = int(sys.argv[3])
    
    #print('w = ',w, 'h=',h, 'ncol=',ncol)
    get_words(w, h, ncol)
