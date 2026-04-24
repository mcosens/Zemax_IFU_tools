'''
make_full_model.py
Written by: Maren Cosens
Date: 4/21/26

Description: uses ZOS-API to generate a full model of MIRMOS IFS mode
-reads in the 5 channel spectrograph design and the 21 slice IFU
-adds configurations and multi-configuration operands that dictate IFU design to generate a model with the full IFU in all spectral channels

To-do: add MOS mode with different field definitions and surface ignores
'''
##import class controlling ZOS-API interface
import zos_pyclass

##run from command line
if __name__ == '__main__':
    zos = zos_pyclass.PythonStandaloneApplication()
    
    # load local variables
    ZOSAPI = zos.ZOSAPI
    IFU_System = zos.TheSystem
    TheApplication = zos.TheApplication
    
    # Setup each file
    IFU_System.LoadFile(r'C:\Users\mcosens\Documents\Zemax\MIRMOS\IFU\virtualIFU_Kband_4-2-26.zmx', False) #old format of lens file
    Full_System=TheApplication.LoadNewSystem(r'C:\Users\mcosens\Documents\Zemax\MIRMOS\IFU\MIRMOS_full.zmx')
    #print("Surfaces in IFU system: " + str(IFU_System.LDE.NumberOfSurfaces)) #for testing that these are not being overwritten
    #print("Surfaces in full system: " + str(Full_System.LDE.NumberOfSurfaces))

    '''
    #add surfaces for a single slice (set parameters that are common now, the rest in the MCE later)
    Full_LDE = Full_System.LDE
    IFS_LDE = IFS_K_System.LDE
    for i in range(0,28):
        Full_LDE.InsertNewSurfaceAt(17+i) #insert 28 new surfaces after Field Lens for IFU
    '''
    #may be easiest to manually add surfaces for 1 configuration then modify multi-configuration operands
    #starting with this for now   
    #add operands to MCE to control IFU
    Full_MCE = Full_System.MCE
    IFS_MCE = IFU_System.MCE
    IFS_nconfigs=IFS_MCE.NumberOfConfigurations
    bands=Full_MCE.NumberOfConfigurations

    # Add and set type for each operand
    start_op_length=Full_MCE.NumberOfOperands
    for i in range(1,IFS_MCE.NumberOfOperands+1):
        Full_MCE.AddOperand()
        current_op=Full_MCE.GetOperandAt(start_op_length+i)
        IFS_op=IFS_MCE.GetOperandAt(i)
        current_op.ChangeType(IFS_op.Type) #set type to match existing operands in K band MCE (i.e. surface parameters, thicknesses, etc.)
        if IFS_op.Param1Enabled:
            #current_op.Param1Enabled = True
            current_op.Param1 = IFS_op.Param1 #set parameters to match existing operands in K band MCE (i.e. surface numbers, etc.)
        if IFS_op.Param2Enabled:
            #current_op.Param2Enabled = True
            current_op.Param2 = IFS_op.Param2
        for j in range(bands):
            current_op.GetOperandCell(j+1).Value = IFS_op.GetOperandCell(1).Value
    
    ##add configurations with parameters for each slice, repeat for each band
    for k in range(bands):
        start_config=1+k*IFS_nconfigs
        for i in range(start_config,start_config+IFS_nconfigs-1): #start at 1 since new config is added after
            Full_MCE.InsertConfiguration(i, False) #insert new configuration with pickup from first configuration
            Full_MCE.SetCurrentConfiguration(i+1)
            IFS_MCE.SetCurrentConfiguration(i-start_config+2)
            Full_MCE.GetOperandAt(start_op_length).GetOperandCell(i+1).Value = ''#remove comment that is automatically copied over (just cleans up editor)
            for j in range(1,IFS_MCE.NumberOfOperands+1):
                IFS_op=IFS_MCE.GetOperandAt(j)
                current_op=Full_MCE.GetOperandAt(start_op_length+j)
                current_op.GetOperandCell(i+1).Value = IFS_op.GetOperandCell(i-start_config+2).Value #set operand values to match existing K band MCE configurations
    
    Full_System.SaveAs("C:\\Users\\mcosens\\Documents\\Zemax\\MIRMOS\\IFU\\MIRMOS_full_IFS.zmx") #save as zmx in order to make plotting functions work in 'make_IFU_plots.py'?

    # close server instance of OpticStudio
    del zos
    zos = None