from QaDataUtils import QaDataUtils


class LsstQaDataUtils(QaDataUtils):

    def __init__(self):
        QaDataUtils.__init__(self)


    def getSourceSetNameList(self):
        """Associate Source accessor names to database columns in a list of pairs. """

        accessors = [
            #["Id",                           "sourceId",                     ],
            ["Ra",                           "ra",                            ],
            ["Dec",                          "decl",                          ],
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
            ["XAstrom",                      "x",                             ],
            #["XAstromErr",                   "xAstromSigma",                 ],
            ["YAstrom",                      "y",                             ],
            #["YAstromErr",                   "yAstromSigma",                 ],
            #["RaAstrom",                     "raAstrom",                     ],
            #["RaAstromErr",                  "raAstromSigma",                ],
            #["DecAstrom",                    "declAstrom",                   ],
            #["DecAstromErr",                 "declAstromSigma",              ],
            #["TaiMidPoint",                  "taiMidPoint",                  ],
            #["TaiRange",                     "taiRange",                     ],
            ["PsfFlux",                      "psfFlux",                       ],
            ["PsfFluxErr",                   "psfFluxSigma",                  ],
            ["ApFlux",                       "apFlux",                        ],
            ["ApFluxErr",                    "apFluxSigma",                   ],
            ["ModelFlux",                    "modelFlux",                     ],
            ["ModelFluxErr",                    "modelFluxSigma",             ],
            #["ModelFlux",                    "flux_ESG",                     ], 
            #["ModelFluxErr",                 "flux_ESG_Sigma",               ],
            ["InstFlux",                     "instFlux",                      ],
            ["InstFluxErr",                     "instFluxSigma",              ],
            #["InstFlux",                      "flux_Gaussian",               ], 
            #["InstFluxErr",                  "flux_Gaussian_Sigma",          ],
            #["NonGrayCorrFlux",              "nonGrayCorrFlux",              ],
            #["NonGrayCorrFluxErr",           "nonGrayCorrFluxSigma",         ],
            #["AtmCorrFlux",                  "atmCorrFlux",                  ],
            #["AtmCorrFluxErr",               "atmCorrFluxSigma",             ],
            #["ApDia",                        "apDia",                        ],
            ["Ixx",                          "shapeIxx",                      ],
            #["IxxErr",                       "ixxSigma",                     ],
            ["Iyy",                          "shapeIyy",                     ],
            #["IyyErr",                       "iyySigma",                    ],
            ["Ixy",                          "shapeIxy",                     ],
            #["IxyErr",                       "ixySigma",                    ],
            #["PsfIxx",                       "psfIxx",                      ],
            #["PsfIxxErr",                    "psfIxxSigma",                 ],
            #["PsfIyy",                       "psfIyy",                      ],
            #["PsfIyyErr",                    "psfIyySigma",                 ],
            #["PsfIxy",                       "psfIxy",                      ],
            #["PsfIxyErr",                    "psfIxySigma",                 ],
            #["Resolution",                   "resolution_SG",               ],
            #["E1",                           "e1_SG",                       ],
            #["E1Err",                        "e1_SG_Sigma",                 ],
            #["E2",                           "e2_SG",                       ],
            #["E2Err",                        "e2_SG_Sigma",                 ],
            #["Shear1",                       "shear1_SG",                   ],
            #["Shear1Err",                    "shear1_SG_Sigma",             ],
            #["Shear2",                       "shear2_SG",                   ],
            #["Shear2Err",                    "shear2_SG_Sigma",             ],
            #["Sigma",                        "sourceWidth_SG",              ],
            #["SigmaErr",                     "sourceWidth_SG_Sigma",        ],
            #["ShapeStatus",                  "shapeStatus",                 ],
            #["Snr",                          "snr",                         ],
            #["Chi2",                         "chi2",                        ],
            #["FlagForAssociation",           "flagForAssociation",          ],
            #["FlagForDetection",             "flagForDetection",            ],
            #["FlagForWcs",                   "flagForWcs",                  ],
            ["FlagPixInterpCen",              "flagPixInterpCen",            ],
            ["FlagNegative",                  "flagNegative",                ],
            ["FlagPixEdge",                   "flagPixEdge",                 ],
            ["FlagBadCentroid",               "flagBadCentroid",             ],
            ["FlagPixSaturCen",               "flagPixSaturCen",             ],
            ["Extendedness",                  "extendedness",                ],
            ]
        return accessors
    

    def getSceNameList(self, dataIdNames, replacements={}):
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
            ['expMidpt',             'expMidpt'             ],
            ['expTime',              'expTime'              ],
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


        for arr in nameList:
            a, b = arr
            if a in replacements:
                arr[1] = replacements[a]

        return nameList
    
