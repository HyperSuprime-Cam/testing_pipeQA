#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
import psycopg2
import lsst.pex.policy as pexPolicy
import time
from lsst.pex.logging import Trace

class DatabaseIdentity:
    """
    Requires file that looks like:
    
    cat ~/.hsc/db-auth.paf
    database: {
        authInfo: {
            host: hsca-01.ipmu.jp
            user: "XXXX"
            password: "YYYY"
        }
    }
    """
    def __init__(self, postgreSqlDb):
        self.sqlSchema   = postgreSqlDb
        self.loadId()

    def loadId(self):
        dbAuth = os.path.join(os.environ["HOME"], ".hsc", "db-auth.paf")
        policy = pexPolicy.Policy(dbAuth)
        authPolicy = policy.get("database").get("authInfo")
        self.sqlDb   = authPolicy.get("db")
        self.sqlUser = authPolicy.get("user")
        self.sqlHost = authPolicy.get("host")
        self.sqlPasswd = authPolicy.get("password")
        self.sqlPort = authPolicy.get("port")        

# Base class
class DatabaseInterface():
    def __init__(self):
        pass


# LSST specific interface
class DbInterface(DatabaseInterface):
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
        self.db     = psycopg2.connect(
            host     = self.dbId.sqlHost,
            database = self.dbId.sqlDb,
            user     = self.dbId.sqlUser,
            password = self.dbId.sqlPasswd,
            port     = self.dbId.sqlPort
            )
        self.cursor = self.db.cursor()
        self.cursor.execute("set search_path to '%s'" % (self.dbId.sqlSchema))


    def execute(self, sql):
        """Execute an sql command

        @param sql Command to be executed.
        """
        Trace("lsst.testing.pipeQA.HscDbInterface", 3, "Executing: %s" % (sql))
        t0 = time.time()

        #print "postgresql>", sql
        
        # forking to handle plotting the summary figures causes (i think)
        # a disconnection when the child exits.  Need to reconnect when
        connected = True
        try:
            self.cursor.execute(sql)
        except Exception, e:
            connected = False

        if not connected:
            self.connect()
            self.cursor.execute(sql)

        results = self.cursor.fetchall()
        t1 = time.time()
        Trace("lsst.testing.pipeQA.HscDbInterface", 2, "Time for SQL query: %.2f s" % (t1-t0))
        return results

    
