# Build instructions

I use pyinstaller to make a single executable (or a folder with an executable in
it).  Then "inno setup" to turn this into a MSI installer.

To make the pyinstaller outputs as small as possible it helps to use a virtual
environment for python.  I think pyinstaller is supposed to only include
necessary libraries, but I've found that it puts in way more than is necessary.
Using a virtual environment avoids that.

NB: I run the following commands in windows in Git bash shell, if you are using
cmd, you will have to modify appropriately

1. From git bash make and activate virtual environment
   ```
   python -m venv ~/Documents/venvs/andytest
   . ~/Documents/venvs/andytest/Scripts/activate
   ```

2. Now we are in the virtual environment, install pyinstaller and the package requirements
   ```
   python -m pip install --upgrade pip
   pip install pyinstaller
   pip install -r requirements.pip
   ```

3. To build dist/andytest/andytest.exe just run
   ```
   pyinstaller -y andytest.spec
   ```

4. If you want to rebuild the andytest.spec file run
   ```
   pyinstaller -y bin/andytest.pyw --onedir --noconsole --add-data="bin/version.txt;." --add-data="bin/help.html;." --icon="etc/andytest.ico" --add-data="etc/andytest.ico;." --add-binary="../code/cellDynamics.cp39-win_amd64.pyd;." --hidden-import "scipy.integrate"
   ```
 
5. The directory `dist/andytest` can be turned into a proper MSI Windows Setup
   file using "Inno Setup" (download it from the internet) and the setup file [inno.iss](inno.iss)
