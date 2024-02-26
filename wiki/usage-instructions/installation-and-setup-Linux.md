# 💻 Linux Installation and Setup

Follow the below instructions to get set up with the Chrono Preprocessor. These instructions were written with Linux users in mind. 

Please note the versions of each piece of software. Newer or alternate versions may work, but have not been tested and verified for compatibility.

<details>

<summary>Step 1: Install FreeCAD</summary>

Install the latest version of FreeCAD. The download is available for free. 
Open a terminal. Run the following commands step by step. 

* apt-get -y update
* apt-get -y install software-properties-common
* add-apt-repository ppa:freecad-maintainers/freecad-stable
* apt-get -y  update
* apt-get -y install mesa-utils libglew-dev freeglut3-dev libgl1-mesa-dri freeca

</details>

<details>

<summary>Step 2: Find FreeCad User Workbench directory </summary>

* Open FreeCAD.
* To find FreeCad User Workbench directory, run following command in FreeCad python panel:  
               “App.getUserAppDataDir()” 



</details>

<details>

<summary>Step 3: Clone ChronoConcrete Workbench from GitHub to FreeCad user application data directory. 
 </summary>


* Open a terminal.
* Clone > "**git clone https://github.com/Concrete-Chrono-Development/chrono-preprocessor chronoConcrete**".
* Check if the Chrono Workbench is available in the list of installed workbenches.

</details>

<details>

<summary>Step 4: Check "tetgen" Installation</summary>

Verify that **tetgen** is installed properly in FreeCAD. 


* Go to " ~/.local/share/FreeCAD/Mod/chronoConcrete/freecad/chronoWorkbench/tetgen"
* Run on a terminal > **./tetgen**

* Check if there is a warning such as > **bash: ./tetgen: Permission denied**
* Run the following command -> "chmod -R 777 *" to give all permission. 
* Run again on a terminal > **./tetgen** 

* Add **tetgen** path in the bash script.
Open the file **bashrc** and put the following line and save. 

**export PATH=$PATH:~/.local/share/FreeCAD/Mod/chronoConcrete/freecad/chronoWorkbench/tetgen**
Run the command on a terminal > "**source ~/.bashrc**"

</details>


<details>

<summary>Step 5: "gmsh" Installation</summary>

* Download **gmsh** from the website 

[https://gmsh.info/bin/Linux/ ]

You can select version 4.4.1.

* Open a terminal and go > **~/gmsh-4.4.1-Linux64/bin** 
* Check by runnung on terminal > **./gmsh**

* Add **gmsh** path in the bash script.
Open the file **bashrc** and put the following line and save.

**export PATH=$PATH:~/gmsh-4.4.1-Linux64/bin**
Run the command on a terminal > "**source ~/.bashrc**"

</details>




