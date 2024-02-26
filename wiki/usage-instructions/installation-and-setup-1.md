# 💻 Desktop Installation and Setup

Follow the below instructions to get set up with the Chrono Preprocessor. These instructions were written with Windows users in mind. It is also possible to install on Mac or Linux.

Please note the versions of each piece of software. Newer or alternate versions may work, but have not been tested and verified for compatibility.

<details>

<summary>Step 1: Install FreeCAD</summary>

Install the latest version of FreeCAD (use at least version 0.20.2). The download is available for free:

[https://www.freecadweb.org/downloads.php](https://www.freecadweb.org/downloads.php)

</details>

<details>

<summary>Step 2: Install Git Client</summary>

Any Git client can be used to push and pull from the GitHub. We recommend using SourceTree and these instructions will assume you are using that installation. You can download SourceTree for free here:

[https://www.sourcetreeapp.com/](https://www.sourcetreeapp.com/)

</details>

<details>

<summary>Step 3: Pull GitHub</summary>

We recommend pulling the GitHub directly into the FreeCAD workbench directory. Otherwise if you pull to another location then you will need to copy the pulled files to the appropriate directory.

* Open Sourcetree
* Select **File** > **Clone / New...**
* Select "**Remote**" and "**Add an account...**"
* For "Hosting Service" select "**GitHub**". For "Authentication" select "**OAuth**"
* Click on "**Refresh OAuth Token**" and login to GitHub and allow Sourcetree in the browser window that opens
* Click **Ok** in Sourcetree. Then "**chrono-preprocessor**" should populate on the right side of the window. If it doesn't, you may need to click refresh.
* Select "**chrono-preprocessor**" and click "**Clone**"
* When filling out the clone window, it should look like the image below (be sure you are naming it chronoConcrete), and with the appropriate username filled out instead of "**\<usr>**"
* Click "**Clone**"

</details>

<figure><img src="../../.gitbook/assets/clone.png" alt=""><figcaption></figcaption></figure>

<details>

<summary>Step 4: Test Installation</summary>

Verify that everything is installed properly by opening FreeCAD and check if the Chrono Workbench is available in the list of installed workbenches.

</details>


# 💻 Linux Installation and Setup

Follow the below instructions to get set up with the Chrono Preprocessor.  


<details>

<summary>Step 1: Install FreeCAD</summary>

Open a terminal. Run the following commands step by step. 

* apt-get -y update
* apt-get -y install software-properties-common
* add-apt-repository ppa:freecad-maintainers/freecad-stable
* apt-get -y  update
* apt-get -y install mesa-utils libglew-dev freeglut3-dev libgl1-mesa-dri freeca

</details>

<details>

<summary>Step 2: Find FreeCad User Workbench directory. </summary>

* Open FreeCAD.
* To find FreeCad User Workbench directory, run following command in FreeCad python panel:  
               “App.getUserAppDataDir()” 



</details>

<details>

<summary>Step 3: Clone ChronoConcrete Workbench from GitHub to FreeCad user application data directory 
obtained in Step 2. </summary>


* Open a terminal.
* Clone the repisotory 

 **"git clone https://github.com/Concrete-Chrono-Development/chrono-preprocessor chronoConcrete"**

* Check if the Chrono Workbench is available in the list of installed workbenches.

</details>

<details>

<summary>Step 4: Check tetgen Installation</summary>

Verify that **tetgen** is installed properly in FreeCAD. 

* Find **tetgen** under FreeCAD folder.
 You may find it in this folder  " ~/.local/share/FreeCAD/Mod/chronoConcrete/freecad/chronoWorkbench/tetgen"
* Run on a terminal > **./tetgen**

* Check if there is a warning such as > **bash: ./tetgen: Permission denied**
 
 Run the following command > "chmod -R 777 * " to give all permission. 
 

* Add **tetgen** path into the bash script.
Open the file **bashrc** and put the following line and save. 

**export PATH=$PATH:~/.local/share/FreeCAD/Mod/chronoConcrete/freecad/chronoWorkbench/tetgen**

Run the command on a terminal > "**source ~/.bashrc**"

</details>


<details>

<summary>Step 5: gmsh Installation</summary>


* Download **gmsh** from the website and extract into an appropriate folder.

[https://gmsh.info/bin/Linux/ ]

You can select version 4.4.1.

* Open a terminal and go to gmsh folder

You may find it in this folder **~/gmsh-4.4.1-Linux64/bin** 

* Check by running on terminal > **./gmsh**

* Add **gmsh** path in the bash script.
Open the file **bashrc** and put the following line and save

**export PATH=$PATH:~/gmsh-4.4.1-Linux64/bin**

Run the command on a terminal > "**source ~/.bashrc**"

</details>

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






