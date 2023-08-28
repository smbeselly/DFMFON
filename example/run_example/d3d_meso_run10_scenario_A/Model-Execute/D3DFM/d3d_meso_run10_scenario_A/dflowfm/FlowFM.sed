[SedimentFileInformation]
    FileCreatedBy         = Deltares, FM-Suite DFlowFM Model Version 4.4.0.0, DFlow FM Version 1.2.110.68456M 
    FileCreationDate      = Sun Jan 09 2022, 11:30:42 
    FileVersion           = 02.00                  
[SedimentOverall]
    Cref                  = 1600                   [kg/m3]   Reference density for hindered settling calculations
[Sediment]
    Name                  = #SedimentMud#                    Name of sediment fraction
    SedTyp                = mud                              Must be "sand", "mud" or "bedload"
    IniSedThick           = 30                     [m]       Initial sediment layer thickness at bed
    FacDss                = 1                      [-]       Factor for suspended sediment diameter
    RhoSol                = 2650                   [kg/m3]   Specific density
    TraFrm                = -3                               Integer selecting the transport formula
    CDryB                 = 500                    [kg/m3]   Dry bed density
    SalMax                = 31                     [ppt]     Salinity for saline settling velocity
    WS0                   = 0.0005                 [m/s]     Settling velocity fresh water
    WSM                   = 0.0005                 [m/s]     Settling velocity saline water
    EroPar                = 5E-05                  [kg/m²s]  Erosion parameter
    TcrSed                = 1000                   [N/m²]    Critical stress for sedimentation
    TcrEro                = 0.3                    [N/m²]    Critical stress for erosion
    TcrFluff              = 4.94065645841247E-324  [N/m²]    Critical stress for fluff layer erosion
