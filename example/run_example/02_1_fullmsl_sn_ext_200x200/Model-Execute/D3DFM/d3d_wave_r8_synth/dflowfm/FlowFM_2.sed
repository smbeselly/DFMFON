[SedimentFileInformation]
    FileCreatedBy         = Deltares, FM-Suite DFlowFM Model Version 4.10.0.3624, DFlow FM Version 1.2.167.141798 
    FileCreationDate      = Tue Jan 17 2023, 20:54:58 
    FileVersion           = 02.00                  
[SedimentOverall]
    Cref                  = 1600                   [kg/m3]   Reference density for hindered settling calculations
    IopSus                = 0                      
[Sediment]
    Name                  = #Mud#                            Name of sediment fraction
    SedTyp                = mud                              Must be "sand", "mud" or "bedload"
    IniSedThick           = 5                      [m]       Initial sediment layer thickness at bed
    FacDss                = 1                      [-]       Factor for suspended sediment diameter
    RhoSol                = 2650                   [kg/m3]   Specific density
    TraFrm                = -3                               Integer selecting the transport formula
    CDryB                 = 500                    [kg/m3]   Dry bed density
    SalMax                = 0                      [ppt]     Salinity for saline settling velocity
    WS0                   = 0.0009                 [m/s]     Settling velocity fresh water
    WSM                   = 0.0009                 [m/s]     Settling velocity saline water
    EroPar                = 0.00007                 [kg/m²s]  Erosion parameter
    TcrSed                = 1000                   [N/m²]    Critical stress for sedimentation
    TcrEro                = 0.24                    [N/m²]    Critical stress for erosion
    TcrFluff              = 4.94065645841247E-324  [N/m²]    Critical stress for fluff layer erosion
