# ☁ Quest Installation and Setup

Follow the below instructions to get set up with Project Chrono on Quest. These instructions assume that you have a Northwestern **NetID** and access to the group project allocation **p31861**. If you do not have either of these, please contact Matthew Troemner at [mtroemner@gmail.com](https://app.gitbook.com/u/xUYe9ATU96ZtZ44QEt0NAkdxV922).

<details>

<summary>Step 1: Install and Configure SSH Client </summary>

Install an SSH Client

* Any SSH client should work, but we recommend PuTTY
* Download and install PuTTY from [https://www.putty.org/](https://www.putty.org/)
* After installation open PuTTY. Within the PuTTY Configuration window:
  * Enter '**quest.northwestern.edu**' for Host Name
  * Enter '**22**' for Port&#x20;
  * Enter '**Quest**' for Saved Sessions
  * Click **Save**
  * Click **Quest** which should not be added to the Session list
  * Click **Open**
* A new SSH window will open. In this window login with your Northwestern NetID and password

</details>

<details>

<summary>Step 2: Install and Configure SFTP Client</summary>

Install an SFTP Client

* Any SFTP client should work, but we recommend FileZilla
* Download and install FileZilla Client from [https://filezilla-project.org/](https://filezilla-project.org/)
* After installation open FileZilla. Within the FileZilla window:
  * Click **File** and then **Site Manager...**
  * In the opened window click **New Site** and enter 'Quest' for the name
  * Enter '**quest.northwestern.edu**' for Host&#x20;
  * Enter '**22**' for Port&#x20;
  * Enter your NetID for **User** and password for **Password**
  * Click **New Bookmark** and enter 'Projects' for the name
  * Choose any Local directory that you want
  * Enter '**/projects/p31861**' for Remote directory
  * Click **OK**
  * Click **File** and then **Site Manager...**
  * Click **Connect**
* The remote site on the right side of your window should automatically connect to the Quest Project Chrono Project and you should see a folder called **Singularity Container**

</details>

<details>

<summary>Step 3: Create Your User Directory</summary>

Create a directory for all of your developments and testing. No files/folders should be created or changed at the top-most '/projects/p31861' directory.

* Within FileZilla, enter the **Users** folder&#x20;
* Right click in the '/projects/p31861/Users' folder and select **Create Directory and Enter It**
* Name the folder with your name in the following format **LastnameFirstname**

</details>

<details>

<summary>Step 4: Copy Singularity File to User Directory</summary>

Copy the SIF file into your User Directory

* In the SSH window run the following command, being sure to replace **LastnameFirstname** with your correct directory name

<pre><code><strong>cp /projects/p31861/SingularityContainer/project-chrono-dependencies.sif /projects/p31861/Users/LastnameFirstname 
</strong></code></pre>

</details>

<details>

<summary>Step 5: Clone Chrono-Concrete</summary>

Clone the Project Chrono GitHub into your User Directory

* In the SSH window cd into your User Directory with the following command, being sure to replace **LastnameFirstname** with your correct directory name&#x20;

<pre><code><strong>cd /projects/p31861/Users/LastnameFirstname
</strong></code></pre>

* Clone the GitHub project here with the following command

```
git clone https://github.com/Concrete-Chrono-Development/chrono-concrete.git
```

* Pull updates to GitHub project

```
cd chrono-concrete
git pull https://github.com/Concrete-Chrono-Development/chrono-concrete.git
git submodule init​
git submodule update
```

</details>

<details>

<summary>Step 6: Build Chrono-Concrete</summary>

Copy example make script, edit, and build Project Chrono

* Copy example make script to User Directory, being sure to replace **LastnameFirstname** with your correct directory name&#x20;

<pre><code><strong>cp /projects/p31861/ExampleScripts/submit_chrono_make.sh /projects/p31861/Users/LastnameFirstname 
</strong></code></pre>

* Navigate to the newly copied 'submit\_chrono\_make.sh' file in FileZilla and double-click on it to edit. Change all instances of **LastnameFirstname** in the file to your appropriate directory and save/upload editted file back to Quest
* In your SSH client navigate to your User Directory and run the following command to submit job

```
sbatch submit_chrono_make.sh
```

You can check the status of your job with the command, being sure to replace **NetID** with your NetID:

```
squeue -u NetID
```

Once the job has completed, proceed to Step 7.

</details>

<details>

<summary>Step 7: Verify Installation</summary>

Verify proper installation of Chrono-Concrete by running a test job with MPI

* Copy example make script to User Directory, being sure to replace **LastnameFirstname** with your correct directory name&#x20;

```
cp /projects/p31861/ExampleScripts/example_submit_mpi.sh /projects/p31861/Users/LastnameFirstname 
```

* Navigate to the newly copied 'submit\_chrono\_make.sh' file in FileZilla and double-click on it to edit. Change all instances of **LastnameFirstname** in the file to your appropriate directory and save/upload editted file back to Quest
* Make an output directory, being sure to replace **LastnameFirstname** with your correct directory name&#x20;

```
mkdir /projects/p31861/Users/LastnameFirstname/outdir
```

* In your SSH client navigate to your User Directory and run the following command to submit job

```
sbatch example_submit_mpi.sh
```

You can check the status of your job with the command, being sure to replace **NetID** with your NetID:

```
squeue -u NetID
```

Once the job has completed, open the outlog file in your User Directory and confirm that the simulation ran. Then navigate to the output directory (./outdir/TestJob) and confirm that several .csv files were created.&#x20;

</details>

<details>

<summary>Step 8: Development</summary>

Code within the chrono-concrete directory can be developed as needed and be pushed/pulled to the GitHub. Please read online about how git works so that you ensure you are properly developing with everyone else.&#x20;

You can modify/copy the example .sh scripts and outdir in your User Directory to help your developments.

Please **do not** edit any files outside of your User Directory.

</details>
