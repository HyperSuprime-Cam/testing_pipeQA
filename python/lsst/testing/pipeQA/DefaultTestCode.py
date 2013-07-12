import sys
import traceback
import os
import re
import inspect
import sqlite
import stat


class TestFailError(Exception):
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)
    


class Test(object):
    """A class to verify some condition is met.
    """

    def __init__(self, label, value, limits, comment):
        """
        @param label   A name for this test
        @param value   Value to be tested
        @param limits  A list [min, max] specifying range of acceptable values (inclusive).
        @param comment A comment with extra info about the test
        """
        
        self.label = label
        self.value = value
        self.limits = limits
        self.comment = comment

    def __str__(self):
        return self.label+" "+str(self.evaluate())+" value="+str(self.value)+" limits="+str(self.limits)

    def evaluate(self):
        """Verify that our value is within our limits."""
        
        # grab a traceback for failed tests
        if (self.value < self.limits[0] or self.value > self.limits[1]):
            return False
        else:
            return True


class TestSet(object):
    """A container for Test objects and associated matplotlib figures."""

    def __init__(self, label=None, group="", clean=False):
        """
        @param label  A name for this testSet
        @param group  A category this testSet belongs to
        """
        
        self.testDir = os.path.join("tests", ".tests")
        self.figDir = os.path.join("figures")

        if not os.path.exists(self.testDir):
            os.mkdir(self.testDir)
            
        testfileName = inspect.stack()[-1][1]
        self.testfileName = os.path.split(testfileName)[1]
        if label is not None:
            self.testfileName += "."+label

        self.status = True
        self.passFile = os.path.join(self.testDir, self.testfileName)
        self.failFile = self.passFile + ".failed"
        self.useFile = self.passFile
        
        self.tests = []
 

    def addTest(self, *args):
        """Add a test to this testing suite.

        @param *args  Either a Test object, or the arguments to create one
        """

        if len(args) >= 4:
            test = Test(*args)
        elif len(args) == 1:
            test, = args

        self.tests.append(test)

        # grab a traceback for failed tests
        backtrace = ""
        message = str(test)
        try:
            if not test.evaluate():
                failMessage = "Failed test '"+test.label+"': " + \
                              "value '" + str(test.value) + "' not in range '" + \
                              str(test.limits)+"'."
                raise TestFailError(failMessage)

        except TestFailError, e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            message = "".join(traceback.format_stack()[:-1]) + "\n" + str(e)

            # if we failed a test, change to the .failed file.
            if self.status:
                os.rename(self.useFile, self.failFile)
                self.useFile = self.failFile
            self.status = False
            
        # enter the test in the output file
        fp = open(self.useFile, 'w')
        fp.write(message+"\n")
        fp.close()


    def addMetadata(self, *args):
        """Associate metadata with this TestSet """
        pass
        
    def addTests(self, testList):
        """Add a list of Test objects to this TestSet."""
        for test in testList:
            self.addTest(test)


    def addFigure(self, fig, filename, caption, **kwargs):
        """Add a figure to this test suite.
        
        @param fig      a matplotlib figure
        @param filename The basename of the figure.
        @param caption  text describing the figure
        """
        
        if not os.path.exists(self.figDir):
            os.mkdir(self.figDir)
        path = os.path.join(self.figDir, filename)
        fig.savefig(path)
        captionPath = os.path.join(self.figDir, filename+".caption")
            
        fp = open(captionPath, 'w')
        fp.write(caption+"\n")
        fp.close()
        
    def importExceptionDict(self):
        pass
    def importLogs(self):
        pass
    def importEupsSetups(self):
        pass


