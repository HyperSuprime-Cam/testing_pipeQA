import os, re
import math

import numpy

import lsst.afw.image   as afwImage
import lsst.afw.coord   as afwCoord
import lsst.pex.logging  as pexLog
import lsst.pex.policy  as pexPolicy
import lsst.meas.astrom as measAstrom
import lsst.meas.algorithms.utils as maUtils

def findDataInTestbed(label, raiseOnFailure=True):
    """Scan TESTBED_PATH directories to find a testbed dataset with a given name

    @param label            Directory name to look for in TESTBED_PATH directories.
    @param raiseOnFailure   Raise an exception if nothing is found.
    """
    
    # If label==None -> use TESTBOT_DIR (_DIR, not _PATH ... only one dataset can be run)
    if re.search("^testBot", label):
        testbotDir = os.getenv("TESTBOT_DIR")
        testbedDir, testdataDir  = os.path.split(testbotDir)
        return testbedDir, testbotDir

    
    # otherwise, get a specific test-data set from one of the testbed directories
    testbedPath = os.getenv("TESTBED_PATH")
    if testbedPath is not None:
        testbedDirs = testbedPath.split(":")
    else:
        raise Exception("Must specify environment variable TESTBED_PATH.")
    
    #############################
    # find the label in the testbed path
    testbedDir = None
    testdataDir = None
    for tbDir in testbedDirs:
        path = os.path.join(tbDir, label)
        if os.path.exists(path):
            testbedDir = tbDir
            testdataDir = path
                
    if raiseOnFailure and (testbedDir is None):
        msg = "Testbed %s was not found in any of the TESTBED_PATH directories:\n" % (label)
        msg += "\n".join(testbedDirs) + "\n"
        raise Exception(msg)
    
    return testbedDir, testdataDir


def setSourceBlobsNone(s):
    """Free the blob structures (photometry,astrometry,shape) in Source to reduce object size. """
    s.setPhotometry(None)
    s.setAstrometry(None)
    s.setShape(None)

def setSourceSetBlobsNone(ss):
    """Free the blob structures (photometry,astrometry,shape) for all Source in SourceSet. """
    for s in ss:
        setSourceBlobsNone(s)


def setMatchListBlobsNone(matchList):
    """Free the blob structures (photometry,astrometry,shape) for all Source in MatchList. """
    
    for s1, s2, d in matchList:
        setSourceBlobsNone(s1)
        setSourceBlobsNone(s2)













        
def getSourceSetNameList():
    """Associate Source accessor names to database columns in a list of pairs. """
    
    accessors = [
        #["Id",                           "sourceId" ], #"objectId"                      ],
        ["Ra",                           "ra",                           ],
        ["Dec",                          "decl",                         ],
        #["RaErrForWcs",                  "raSigmaForWcs",                ],
        #["DecErrForWcs",                 "declSigmaForWcs",              ],
        #["RaErrForDetection",            "raSigmaForDetection",          ],
        #["DecErrForDetection",           "declSigmaForDetection",        ],
        #["XFlux",                        "xFlux",                        ],
        #["XFluxErr",                     "xFluxSigma",                   ],
        #["YFlux",                        "yFlux",                        ],
        #["YFluxErr",                     "yFluxSigma",                   ],
        #["RaFlux",                       "raFlux",                       ],
        #["RaFluxErr",                    "raFluxSigma",                  ],
        #["DecFlux",                      "declFlux",                     ],
        #["DecFluxErr",                   "declFluxSigma",                ],
        #["XPeak",                        "xPeak",                        ],
        #["YPeak",                        "yPeak",                        ],
        #["RaPeak",                       "raPeak",                       ],
        #["DecPeak",                      "declPeak",                     ],
        ["XAstrom",                      "x",                      ],
        #["XAstromErr",                   "xAstromSigma",                 ],
        ["YAstrom",                      "y",                      ],
        #["YAstromErr",                   "yAstromSigma",                 ],
        #["RaAstrom",                     "raAstrom",                     ],
        #["RaAstromErr",                  "raAstromSigma",                ],
        #["DecAstrom",                    "declAstrom",                   ],
        #["DecAstromErr",                 "declAstromSigma",              ],
        #["TaiMidPoint",                  "taiMidPoint",                  ],
        #["TaiRange",                     "taiRange",                     ],
        ["PsfFlux",                      "psfFlux",                      ],
        ["PsfFluxErr",                   "psfFluxSigma",                 ],
        ["ApFlux",                       "apFlux",                       ],
        ["ApFluxErr",                    "apFluxSigma",                  ],
        ["ModelFlux",                    "modelFlux",                    ],
        ["ModelFluxErr",                    "modelFluxSigma",                    ],
        #["ModelFlux",                    "flux_ESG",                 ], 
        #["ModelFluxErr",                 "flux_ESG_Sigma",               ],
        ["InstFlux",                     "instFlux",                     ],
        ["InstFluxErr",                     "instFluxSigma",                     ],
        #["InstFlux",                      "flux_Gaussian",                     ], 
        #["InstFluxErr",                  "flux_Gaussian_Sigma",                ],
        #["NonGrayCorrFlux",              "nonGrayCorrFlux",              ],
        #["NonGrayCorrFluxErr",           "nonGrayCorrFluxSigma",         ],
        #["AtmCorrFlux",                  "atmCorrFlux",                  ],
        #["AtmCorrFluxErr",               "atmCorrFluxSigma",             ],
        #["ApDia",                        "apDia",                        ],
        ["Ixx",                          "shapeIxx",                          ],
        #["IxxErr",                       "ixxSigma",                     ],
        ["Iyy",                          "shapeIyy",                          ],
        #["IyyErr",                       "iyySigma",                     ],
        ["Ixy",                          "shapeIxy",                          ],
        #["IxyErr",                       "ixySigma",                     ],
        #["PsfIxx",                       "psfIxx",                       ],
        #["PsfIxxErr",                    "psfIxxSigma",                  ],
        #["PsfIyy",                       "psfIyy",                       ],
        #["PsfIyyErr",                    "psfIyySigma",                  ],
        #["PsfIxy",                       "psfIxy",                       ],
        #["PsfIxyErr",                    "psfIxySigma",                  ],
        #["Resolution",                   "resolution_SG",                ],
        #["E1",                           "e1_SG",                        ],
        #["E1Err",                        "e1_SG_Sigma",                  ],
        #["E2",                           "e2_SG",                        ],
        #["E2Err",                        "e2_SG_Sigma",                  ],
        #["Shear1",                       "shear1_SG",                    ],
        #["Shear1Err",                    "shear1_SG_Sigma",              ],
        #["Shear2",                       "shear2_SG",                    ],
        #["Shear2Err",                    "shear2_SG_Sigma",              ],
        #["Sigma",                        "sourceWidth_SG",               ],
        #["SigmaErr",                     "sourceWidth_SG_Sigma",         ],
        #["ShapeStatus",                  "shapeStatus",                  ],
        #["Snr",                          "snr",                          ],
        #["Chi2",                         "chi2",                         ],
        #["FlagForAssociation",           "flagForAssociation",           ],
        #["FlagForDetection",             "flagForDetection",             ],
        #["FlagForWcs",                   "flagForWcs",                   ],
        ["FlagPixInterpCen",              "flagPixInterpCen",             ],
        ["FlagNegative",                  "flagNegative",             ],
        ["FlagPixEdge",                   "flagPixEdge",             ],
        ["FlagBadCentroid",               "flagBadCentroid",             ],
        ["FlagPixSaturCen",               "flagPixSaturCen",             ],
        ["Extendedness",                  "extendedness",                 ],
        ]
    return accessors


def getSourceSetAccessors():
    """Get a list of all accessor names for Source objects. """
    return zip(*getSourceSetNameList())[0]
def getSourceSetDbNames(replacementDict):
    """Get a list of all database names for Source objects. """
    rawList = list(zip(*getSourceSetNameList())[1])
    for k,v in replacementDict.items():
        matches = [i for i,x in enumerate(rawList) if x == k]
        arg = matches[0]
        rawList[arg] = v
    return rawList
        
def getCalexpNameLookup():
    """Associate calexp_md names to database names."""
    
    nameLookup = {
        'FILTER':           'filterName'           ,
        'RA':                   'ra'               ,
        'DEC':                 'decl'              ,
        'CRPIX1':               'crpix1'           ,
        'CRPIX2':               'crpix2'           ,
        'CRVAL1':               'crval1'           ,
        'CRVAL2':               'crval2'           ,
        'CD1_1':                'cd1_1'            ,
        'CD1_2':                'cd1_2'            ,
        'CD2_1':                'cd2_1'            ,
        'CD2_2':                'cd2_2'            ,
        'FLUXMAG0':             'fluxMag0'         ,
        'FLUXMAG0ERR':        'fluxMag0Sigma'      ,
        'SEEING':                 'fwhm'           ,
        }

    return nameLookup


def getSceNameList(dataIdNames):
    """Associate SourceCcdExposure names to database columns in a list of pairs. """
    
    nameList = [ ['scienceCcdExposureId', 'scienceCcdExposureId' ] ] + dataIdNames
    nameList += [
        #['visit',                'visit'                ],
        #['raft',                 'raft'                 ],
        #['raftName',             'raftName'             ],
        #['ccd',                  'ccd'                  ],
        #['ccdName',              'ccdName'              ],
        ['filterId',             'filterId'             ],
        ['filterName',           'filterName'           ],
        ['ra',                   'ra'                   ],
        ['decl',                 'decl'                 ],
        #['equinox',              'equinox'              ],
        #['raDeSys',              'raDeSys'              ],
        #['ctype1',               'ctype1'               ],
        #['ctype2',               'ctype2'               ],
        ['crpix1',               'crpix1'               ],
        ['crpix2',               'crpix2'               ],
        ['crval1',               'crval1'               ],
        ['crval2',               'crval2'               ],
        ['cd1_1',                'cd1_1'                ],
        ['cd1_2',                'cd1_2'                ],
        ['cd2_1',                'cd2_1'                ],
        ['cd2_2',                'cd2_2'                ],
        #['llcRa',                'llcRa'                ],
        #['llcDecl',              'llcDecl'              ],
        #['ulcRa',                'ulcRa'                ],
        #['ulcDecl',              'ulcDecl'              ],
        #['urcRa',                'urcRa'                ],
        #['urcDecl',              'urcDecl'              ],
        #['lrcRa',                'lrcRa'                ],
        #['lrcDecl',              'lrcDecl'              ],
        #['taiMjd',               'taiMjd'               ],
        #['obsStart',             'obsStart'             ],
        #['expMidpt',             'expMidpt'             ],
        #['expTime',              'expTime'              ],
        #['nCombine',             'nCombine'             ],
        #['binX',                 'binX'                 ],
        #['binY',                 'binY'                 ],
        #['readNoise',            'readNoise'            ],
        #['saturationLimit',      'saturationLimit'      ],
        #['gainEff',              'gainEff'              ],
        ['fluxMag0',             'fluxMag0'             ],
        ['fluxMag0Sigma',        'fluxMag0Sigma'        ],
        ['fwhm',                 'fwhm'                 ],
        ]
    return nameList

def getSceDbNames(dataIdNames):
    """Get SourceCcdExposure database column names."""
    return zip(*getSceNameList(dataIdNames))[1]




def getCalibObjects(butler, filterName, dataId):
    """
    A version of getCalibObjectsImpl that isn't a class method, for use by other code
    
    @param useOutputSrc             # use fluxes from the "src", not "icSrc"
    """

    #################################
    # color calibration
    #
    # this stuff stolen from pipette.
    # - numbers from suprimecam.paf, method in Calibrate.py
    #####################
    primaryLookup = {
        #"B" : "g",
        #"V" : "g",
        #"R" : "r",
        #"I" : "i",
        "g" : "g",
        "g" : "g",
        "r" : "r",
        "i" : "i",
        }
    # The below colour terms are from the last page of
    # http://www.naoj.org/staff/nakata/suprime/illustration/colorterm_report_ver3.pdf
    secondaryLookup = {
        "g" : ["r", [-0.00569, -0.0427]],
        "r" : ["g", [0.00261, 0.0304]],
        "i" : ["r", [0.00586, 0.0827, -0.0118]],
        "z" : ['i', [0.000329, 0.0608, 0.0219]],
        "y" : ['i', [0.000329, 0.0608, 0.0219]],
        }
    # figure out which filters we need for the color term
    primaryFilter   = primaryLookup[filterName]      # Primary band for correction
    secondaryFilter = secondaryLookup[primaryFilter][0]  # Secondary band for correction



    log = pexLog.Log.getDefaultLog()
    psources  = butler.get('icSrc', dataId)
    pmatches  = butler.get('icMatch', dataId)
    calexp_md = butler.get('calexp_md', dataId)

    wcs = afwImage.makeWcs(calexp_md)
    W, H = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")
    calib = afwImage.Calib(calexp_md)

    matches = pmatches.getSourceMatches()
    setMatchListBlobsNone(matches)
    sources = psources.getSources() ## fails
    setSourceSetBlobsNone(sources)

    anid = pmatches.getSourceMatchMetadata().getInt('ANINDID')

    del psources; del pmatches; del calexp_md # cleanup

    useOutputSrc = False
    if useOutputSrc:
        srcs = butler.get('src', dataId).getSources()
        setSourceSetBlobsNone(srcs)
        
        import lsst.afw.detection as afwDetect
        pmMatch = afwDetect.matchXy(sources, srcs, 1.0, True)
        for icSrc, src, d in pmMatch:
            icSrc.setPsfFlux(src.getPsfFlux())

    # ref sources
    xc, yc = 0.5*W, 0.5*H
    radec = wcs.pixelToSky(xc, yc)
    ra = radec.getLongitude(afwCoord.DEGREES)
    dec = radec.getLatitude(afwCoord.DEGREES)
    radius = wcs.pixelScale().asDegrees()*math.hypot(xc, yc)*1.1

    pol = pexPolicy.Policy()
    pol.set('matchThreshold', 30)
    solver = measAstrom.createSolver(pol, log)
    idName = 'id'

    X = solver.getCatalogue(ra, dec, radius, primaryFilter, idName, anid)
    refsources = X.refsources
    setSourceSetBlobsNone(refsources)
    
    inds = X.inds

    referrs, stargal = None, None
    cols = solver.getTagAlongColumns(anid)
    colnames = [c.name for c in cols]

    col = 'starnotgal'
    if col in colnames:
        stargal1 = solver.getTagAlongBool(anid, col, inds)
        stargal = []
        for i in range(len(stargal1)):
            stargal.append(stargal1[i])

    fdict = maUtils.getDetectionFlags()


    # get the data for the secondary
    XX = solver.getCatalogue(ra, dec, radius, secondaryFilter, idName, anid)
    secondaryRefsources = XX.refsources
    setSourceSetBlobsNone(secondaryRefsources)

    polyData = secondaryLookup[primaryFilter][1]   # Polynomial correction
    polyData.reverse()                             # Numpy wants decreasing powers
    polynomial = numpy.poly1d(polyData)


    # We already have the 'primary' magnitudes in the matches
    secondariesDict = dict()
    for s in secondaryRefsources:
        secondariesDict[s.getId()] = (s.getPsfFlux(), s.getPsfFluxErr())
    del secondaryRefsources


    keepref = []
    keepi = []
    for i in xrange(len(refsources)):
        ra, dec = refsources[i].getRa(), refsources[i].getDec() # ra,dec in Rads
        x, y = wcs.skyToPixel(afwCoord.Coord(afwGeom.PointD(numpy.degrees(ra), numpy.degrees(dec))))

        if x < 0 or y < 0 or x > W or y > H:
            continue

        refsources[i].setXAstrom(x)
        refsources[i].setYAstrom(y)
        if stargal[i]:
            refsources[i].setFlagForDetection(refsources[i].getFlagForDetection() | fdict["STAR"])

        # color term
        primaryMag = -2.5*numpy.log10(refsources[i].getPsfFlux())
        secondaryMag = -2.5*numpy.log10(secondariesDict[refsources[i].getId()][0])
        diff = polynomial(primaryMag - secondaryMag)
        refsources[i].setPsfFlux(10.0**(-0.4*(primaryMag+diff)))

        keepref.append(refsources[i])
        keepi.append(i)


    refsources = keepref

    if referrs is not None:
        referrs = [referrs[i] for i in keepi]
    if stargal is not None:
        stargal = [stargal[i] for i in keepi]

    measAstrom.joinMatchList(matches, refsources, first=True, log=log)
    measAstrom.joinMatchList(matches, sources, first=False, log=log)

    return matches, calib, refsources



def calibFluxError(f, df, f0, df0):
    if numpy.isinf(f):
        # issues if f is inf
        return numpy.NaN
    else:
        # or df/f, even if they are finite
        try:
            df/f
        except:
            return numpy.NaN

    if f > 0.0 and f0 > 0.0:
        return (df/f + df0/f0)*f/f0
    else:
        return numpy.NaN

def atEdge(bbox, x, y):

    borderWidth = 18
    x0, y0, x1, y1 = bbox
    imgWidth  = x1 - x0
    imgHeight = y1 - y0
    
    if x < borderWidth or imgWidth - x < borderWidth:
        return True
    if y < borderWidth or imgHeight - y < borderWidth:
        return True
    
    return False



