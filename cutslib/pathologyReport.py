from __future__ import print_function
from __future__ import absolute_import
from past.builtins import basestring

import moby2
import pandas as pd
import seaborn as sns
sns.set(style="white")
from .pathologies_tools import pathoList, get_pwv
from matplotlib import pyplot as plt


class pathoReport(object):
    """Provide tools to analyze a pathologies report"""

    def __init__( self, filename):
        pl = pathoList(filename)
        self.data = pd.DataFrame.from_dict(pl.data)
        self.data.index = pd.to_datetime(self.data.ctime,unit='s')
        self.header = pl.header
        self.ndata = pl.ndata
        self.filename = filename
        self.reportname = filename.split('/')[-1].split('.')[0]
        self.livekeys = [k for k in list(self.data.keys()) if k[-4:] == 'Live']

    def addPWV(self):
        self.data['PWV'] = get_pwv(self.data.ctime)


    def select_time(self, time):
        """Select a subset of the data by time. The subset can be chosen as:
           - a year: 'YYYY'
           - a month: 'YYYY-MM'
           - a day: 'YYYY-MM-DD'
           - an hour: 'YYYY-MM-DD HH' (etc.)
           - a range between two years/months/days/etc.
        """
        if not hasattr(self,'data_original'): self.data_original = self.data.copy()
        if type(time) is tuple:
            self.data = self.data[time[0]:time[1]]
        else:
            self.data = self.data[time]
        print("%i TODs have been discarded. %i remain." %(self.ndata-self.data.shape[0],
                                                         self.data.shape[0]))
        self.ndata = self.data.shape[0]


    def select_hours(self,hour_min, hour_max):
        """Select a subset of the data by hour of the day.
        """
        if not hasattr(self,'data_original'): self.data_original = self.data.copy()
        self.data = self.data.between_time(hour_min,hour_max)
        print("%i TODs have been discarded. %i remain." %(self.ndata-self.data.shape[0],
                                                         self.data.shape[0]))
        self.ndata = self.data.shape[0]


    def select_tods(self, todlist):
        """Select a subset of the data from a list of TODs
        """
        if not hasattr(self,'data_original'): self.data_original = self.data.copy()
        sel = self.data.todName.isin(todlist)
        self.data = self.data[sel]
        self.ndata = self.data.shape[0]

    def select_condition(self,seldict):
        """Select TODs based on criteria

        Argument should be passed as a dictionnary of the form { crit: ('lt',val) }
        (use 'lt' for < and 'gt' for >)"""
        crit = list(seldict.keys())
        if not hasattr(self,'data_original'): self.data_original = self.data.copy()
        for c in crit:
            if seldict[c][0] == 'lt':
                self.data = self.data[self.data[c] < seldict[c][1]]
            elif seldict[c][0] == 'gt':
                self.data = self.data[self.data[c] > seldict[c][1]]
        print("%i TODs have been discarded. %i remain." %(self.ndata-self.data.shape[0],
                                                         self.data.shape[0]))

        self.ndata = self.data.shape[0]


    def revert_selection(self):
        if not hasattr(self,'data_original'):
            print('No selection has been aplied')
        else:
            self.data = self.data_original.copy()
            self.ndata = self.data.shape[0]
            del self.data_original
        print("%i TODs." %self.ndata)

    def drop_duplicates(self):
        """Remove duplicates"""
        self.data.drop_duplicates(inplace=True)
        self.ndata = self.data.shape[0]

    def seasonplot(self, crit,
                   time_range = None,
                   dets_lim = (0,1056), pwv_lim = (0,5), pwv_max = 3.,
                   figure=None,
                   filename=None, **kwargs):
        if figure is None: fig = plt.figure(figsize=(30,10))
        else: plt.figure(figure)
        if 'PWV' in self.data.columns:
            ax2 = self.data.PWV.plot(marker='.', ls='-', alpha=0.3,
                                     color='grey',secondary_y=True)
            ax2.set_ylim(pwv_lim)
            ax2.set_ylabel('PWV [mm]')
        ax1 = self.data[crit].plot(ls='none', marker='.', **kwargs)
        ax1.set_ylim(dets_lim)
        if time_range is not None:
            ax1.set_xlim( pd.Timestamp(time_range[0]), pd.Timestamp(time_range[1]))
        ax1.set_ylabel('Number of live detectors')
        ax1.set_xlabel('Date')
        plt.title('%s - %s' %(self.reportname,crit))
        if filename is None:
            plt.draw()
        else:
            fig.savefig(filename)
            plt.close()


    def killedbyplot(self, filename=None, dets_lim = (0,1056), type="violin"):
        """Show the distribution of live detectors after
        applying each criteria individually"""
        if not hasattr(self, 'livedata'):
            self.livedata = self.data[self.livekeys]
        fig = plt.figure(figsize=(10,6))
        if type=="violin":
            ax = sns.violinplot(data=self.livedata, color='grey', fliersize=0.15,
                                orient='h', scale='width',inner='box')
        elif type=="box":
            ax = sns.boxplot(data=self.livedata, color='grey', orient='h', fliersize=0.3)
        else:
            raise NotImplemented("plot type not implemented")
        ax.set_title(self.reportname,fontsize='xx-large')
        ax.set_xlim(dets_lim)
        ax.set_xlabel('Number of live detectors')
        fig.add_axes(ax)
        if filename is None:
            plt.draw()
        else:
            fig.savefig(filename)
            plt.close(fig)



    def scatter_plot(self, crit1, crit2, x_lim=None, y_lim=None,filename=None):
        # self.data.plot(x=crit1,y=crit2,
        #                kind='hexbin',gridsize=50,
        #                bins='log',
        #                cmap='bone_r',colorbar=False)
        # plt.figure()
        # fig = plt.figure(figsize=(15,15))
        sns.jointplot(self.data[crit1],self.data[crit2],
                      xlim=x_lim, ylim=y_lim,
                      kind='hex', cmap='bone_r',
                      stat_func=None,bins='log')
        if filename is None:
            plt.draw()
        else:
            plt.savefig(filename)
            plt.close('all')



    def gridlive(self):
        plt.ioff()
        if not hasattr(self, 'livedata'):
            livekeys = [k for k in list(self.data.keys()) if k[-4:] == 'Live']
            self.livedata = self.data[livekeys]
        g = sns.PairGrid(self.livedata,diag_sharey=False,despine=False)
        with sns.color_palette("bone_r"):
            g.map_diag(sns.kdeplot, lw=3)
            g.map_lower(plt.hexbin, cmap="bone_r", bins='log',gridsize=15,
                        edgecolors='none')
            g.map_upper(plt.scatter, marker='.', alpha=0.1)
        plt.savefig('toto.png')
        plt.close()
