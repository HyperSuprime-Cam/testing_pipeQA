import os
import lsst.pex.policy as pexPolicy
import time
from lsst.pex.logging import Trace

class DatabaseIdentity:
    """
    Requires file that looks like:
    
    cat ~/.lsst/db-auth.paf
    database: {
        authInfo: {
            host: lsst10.ncsa.uiuc.edu
            user: "XXXX"
            password: "YYYY"
        }
    }
    """
    def __init__(self, mySqlDb):
        self.mySqlDb   = mySqlDb
        self.loadId()

    def loadId(self):
        dbAuth = os.path.join(os.environ["HOME"], ".lsst", "db-auth.paf")
        policy = pexPolicy.Policy(dbAuth)
        authPolicy = policy.get("database").get("authInfo")
        self.mySqlUser = authPolicy.get("user")
        self.mySqlHost = authPolicy.get("host")
        self.mySqlPasswd = authPolicy.get("password")
        

# Base class
class DatabaseInterface():
    def __init__(self):
        pass


# LSST specific interface
class LsstSimDbInterface(DatabaseInterface):
    # Mapping from filter names to database names
    filterMap = { "u" : 0, "g" : 1, "r" : 2, "i" : 3, "z" : 4 }

    def __init__(self, dbId):
        """
        @param dbId  A databaseIdentity object contain connection information
        """
        self.dbId = dbId
        DatabaseInterface.__init__(self)

        self.connect()


    def connect(self):
        import MySQLdb
        self.db     = MySQLdb.connect(
            host   = self.dbId.mySqlHost,
            db     = self.dbId.mySqlDb,
            user   = self.dbId.mySqlUser,
            passwd = self.dbId.mySqlPasswd
            )
        self.cursor = self.db.cursor()


    def execute(self, sql):
        """Execute an sql command

        @param sql Command to be executed.
        """
        Trace("lsst.testing.pipeQA.LsstSimDbInterface", 3, "Executing: %s" % (sql))
        t0 = time.time()

        #print "mysql>", sql
        
        
        # forking to handle plotting the summary figures causes (i think)
        # a disconnection when the child exits.  Need to reconnect when
        connected = True
        try:
            self.cursor.execute(sql)
        except Exception, e:
            connected = False

        if not connected:
            #print "Mysql connection broken.  Reconnecting."
            self.connect()
            self.cursor.execute(sql)

        results = self.cursor.fetchall()
        t1 = time.time()
        #print " (t=%.2fs) " % (t1-t0)
        Trace("lsst.testing.pipeQA.LsstSimDbInterface", 2, "Time for SQL query: %.2f s" % (t1-t0))
        return results

    
