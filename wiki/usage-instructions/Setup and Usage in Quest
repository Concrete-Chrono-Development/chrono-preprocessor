# 💻 Quest  FreeCAD Singularity Setup and Usage
- **Run FreeCAD using following command:**
>
>singularity exec -B /projects:/projects -B /software:/software -B /home:/home /hpc/software/freecad/0.21.2/freecad_23.10.sif freecad $@.
>
- **After opening FreeCAD, from tools menu find Addon Manager and install Plot Workbench.**
>
>
- **Open python panel inside FreeCAD and run following command:**
> 
 >"App.getUserAppDataDir()"
>
>This command gives you the path where you need to clone Chrono preprocessor.
>	in my case, the path for new module is: "/home/netID/.local/share/FreeCAD/Mod/ 
>
- **After cloning, get and install gmsh software into your home folder:**
>	
>you can directly download exe file from https://gmsh.info/bin/Linux/
>and copy it into your home folder.
>
- **Check tetgen path inside your freecad folder:**
>	in my case i found it in "/home/netID/.local/share/FreeCAD/Mod/chronoConcrete/freecad/chronoWorkbench/tetgen"
>
- **You need to add tetgen and gmsh path into singularity conatiner. I used --env option and modify command for running singularity as follows:**
>
>	"singularity exec --env PATH=/usr/local/bin:/usr/bin:/bin:/home/NetID/.local/share/gmsh/bin:/home/NetID/.local/share/FreeCAD/Mod/chronoConcrete/freecad/chronoWorkbench/tetgen -B /projects:/projects -B /software:/software -B /home/NetID:/home/NetID /hpc/software/freecad/0.21.2/freecad_23.10.sif  gmsh $@"
>
- **If you experiences any problem, you may need to change the file accesibility for tetgen and gmsh exe:**
>	If it is the case, go to folder that executable is inside and run following code from command prompt
>	
>	"chmod -R 777 *"
>	
- **You can run the freecad using following command:**
>
>	"module load singularity"   
>	"singularity exec --env PATH=/usr/local/bin:/usr/bin:/bin:/home/NetID/.local/share/gmsh/bin:/home/NetID/.local/share/FreeCAD/Mod/chronoConcrete/freecad/chronoWorkbench/tetgen -B /projects:/projects -B /software:/software -B /home/NetID:/home/NetID /hpc/software/freecad/0.21.2/freecad_23.10.sif  gmsh $@"
