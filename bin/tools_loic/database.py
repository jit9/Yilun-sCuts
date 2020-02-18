import MySQLdb
import pandas as pd
import moby2

host = 'morfeo.rolonet'
port = 3306
user = 'act'
passwd = 'act!now'
db = 'actpol_manifest'


class Database:
    """Pandas based class for ACTPol database"""

    def __init__(self, season=2016):
        self.mysql_cn = MySQLdb.connect(
            host=host,
            port=port,
            user =user,
            passwd=passwd,
            db=db)

        self.season = season
        self.query = None
        if self.season == 2013:
            self._prefix = 's1_'
        elif self.season == 2014:
            self._prefix = 's2_'
        elif self.season == 2015:
            self._prefix = 's3_'
        elif self.season == 2016:
            self._prefix = 's4_'
        elif self.season == 2017:
            self._prefix = 's5_'

    def _table(self, table):
        if table == 'tods':
            return 'tod_info'
        elif table == 'acqs':
            return 'mce_data_acq'
        

    def set_season(self, season):
        """Set the season"""
        self.season = season
    
    def describe(self,table):
        """Describe the fields of the database
        Argument is either 'tods' or 'acq'
        """
        if table not in ['tods', 'acqs']:
            print("Argument should be either 'tods' or 'acqs'")
        else:
            print(pd.read_sql('describe %s%s' %(self._prefix,self._table(table)),
                              self.mysql_cn))
    
    def create_query(self, table, column=None):
        """Create a new query from the database

        Argument:
        - table (either 'tods' or 'acqs')

        Option:
        - column (use describe(table) to seee tha available columns
          by default, all columns are loaded
          should be expressed as 'col' or 'col1, col2, col3, ...'
        """
        if column is None:
            column = '*'
        self.query = "select %s from %s%s" %(column, self._prefix, self._table(table))
        self.condition = []
        self.table = table

    def add_condition( self, condition ):
        """Add condition to the existing query

        The condition should be expressed as a tuple:
        (crit, cond, val)
        eg.: (array, =, 'AR1')
             (ctime_start, >, 1441900000)
        """
        if len(self.condition) == 0:
            self.query += " where "
        else:
            self.query += " and "
        self.query += "%s %s '%s' " %(condition[0], condition[1], condition[2])
        self.condition.append(condition)

    def remove_condition( self, condition ):
        """Remove condition to the existing query

        The condition should be 'all' or a tuple:
        (crit, cond, val)
        eg.: (array, =, 'AR1')
             (ctime_start, >, 1441900000)
        """
        if condition == 'all':
            self.condition = []
        else:
            self.condition.remove(condition)

    def make_query(self):
        """Load the data using query"""
        if self.query is None:
            raise RuntimeError("A query needs to be created first")
        self.data = pd.read_sql(self.query, self.mysql_cn)
        self.data = self.data[self.data.datafile_id.notnull()]
        if 'ctime' in list(self.data.keys()):
            k = 'ctime'
        elif 'ctime_start' in list(self.data.keys()):
            k = 'ctime_start'
        self.data.index = pd.to_datetime(self.data[k], unit='s')
        if not hasattr(self.data, 'name'):
            self.load_names()

    def load_names(self):
        """Load tod name for 2013, 2014 and 2015 (not contained in tod_info)"""
        query = 'select id, name from %sdatafiles' %self._prefix
        datafiles = pd.read_sql(query, self.mysql_cn)
        datafiles['name'] = pd.Series([n.strip('.zip') for n in datafiles.name])
        if hasattr(self,'data'):
            self.data = pd.merge(self.data, datafiles,
                                 left_on='datafile_id', right_on='id', how='left')
            if 'ctime' in list(self.data.keys()):
                k = 'ctime'
            elif 'ctime_start' in list(self.data.keys()):
                k = 'ctime_start'
            self.data.index = pd.to_datetime(self.data[k], unit='s')

        if hasattr(self,'tods'):
            self.tods = pd.merge(self.tods, datafiles,
                                 left_on='datafile_id', right_on='id', how='left')

    def load_tods(self):
        """Load all TODs information"""
        request = 'select * from %stod_info;' %self._prefix
        self.tods = pd.read_sql(request, con=self.mysql_cn)
        self.tods = self.tods[self.tods.datafile_id.notnull()]
        self.tods.index = pd.to_datetime(self.tods['ctime_start'], unit='s')
        self.tods.length = self.tods['ctime_end'] - self.tods['ctime_start']
        if not hasattr(self.tods, 'name'):
            self.load_names()

    def load_acqs(self):
        """Load all mce acquisitions information"""
        request = 'select * from %smce_data_acq;' %self._prefix
        self.acqs = pd.read_sql(request, con=self.mysql_cn)
        self.acqs.index = pd.to_datetime(self.acqs['ctime'], unit='s')
        
    def load_datafiles(self):
        """Load all datafiles information"""
        request = 'select * from %sdatafiles;' %self._prefix
        self.datafiles = pd.read_sql(request, con=self.mysql_cn)

        
    def add_column(self,column):
        """Add a (several) column(s) to the database

        Argument:
        ---------
        column: should be expressed as 'col' or 'col1, col2, col3, ...'
        """
        query = "select datafile_id, %s from %s%s" %(
            column, self._prefix, self._table(self.table))
        new_columns = pd.read_sql(query, self.mysql_cn)
        new_columns = new_columns[new_columns.datafile_id.notnull()]
        self.data = pd.merge(self.data, new_columns,
                             on='datafile_id', how='left')

