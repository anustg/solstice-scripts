import matplotlib.pyplot as plt
import numpy as N


def plot_layout(X,Y, savedir, fmt='simple'):

    '''
    Arguements:
    X,Y coordinates of each heliostat
    savedir - the directory of saving the figure    
    '''
    if fmt=='simple':
 
        plt.figure(figsize=(20,20))
        plt.plot(X,Y,'.')
        plt.grid()
        plt.savefig(open('%s/layout.png'%savedir,'w'), bbox_inches='tight')
        plt.close()          
    

    else:
        
        plt.figure(figsize=(20,20))
        plt.plot(X,Y,'.')

        #plt.xlabel('X (m)', fontsize=24)
        #plt.ylabel('Y (m)', fontsize=24)
        #plt.xticks(fontsize=24)
        #plt.yticks(fontsize=24)
        #plt.axis('equal')
        #plt.axis([-N.max(X)*1.2, N.max(X)*1.2, -N.max(Y)*1.2, N.max(Y)*1.2])
        plt.grid()
        #plt.savefig(open('%s/layout.png'%savedir,'w'), bbox_inches='tight', dpi=500)
        plt.close()          

