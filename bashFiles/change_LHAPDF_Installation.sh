# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  cange_LHAPDF_Installation.sh                                                                             #
#  Author: Alexander Epping: a_eppi01@uni-muenster.de                                                       #
#  6 Oct 2021                                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #    
#  Convenient program to change the currently used installation of LHAPDF.                                  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  - Furthermore it will save the current state of the .bash_aliases file in a backup file                  #
#  - The .bash_aliases file should be constucted in a similar way to the one below to work properly.        #
#    My current .bash_aliases file can be also found in the Github.                                         #
#  - The paths obiously should be changed accordingly.                                                      # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# line in the .bash_aliases file in wich the varibale for Apfel is saved
changeLine=3
programName="LHAPDF"

# set "change" as the default value of the input
# possible inputs are mod and std
var=${1:-change}


# change to the folder with .bash_aliases
cd /home/users/a_eppi01

# bash file to be changed
file='.bash_aliases'

# temporary file to copy contents of 'file' to
temp='temp_bash_aliases.txt'

i=1
while read line; do
    if [ $i == $changeLine ]; then 
        if [ $var == change ]; then

            # find the beginning of the definition of the variable
            j="$(expr index "$line" "=")"
            j=$((j+1))

            # look which installation is active and change it to the other one
            case ${line:j:3} in
                "std") echo ${line/std/mod} >> $temp && i=100;;
                "mod") echo ${line/mod/std} >> $temp && i=101;;
            esac
        
        # if there was a parameter given:
        elif [ $var == mod ]; then echo ${line/std/mod} >> $temp && i=102
        elif [ $var == std ]; then echo ${line/mod/std} >> $temp && i=103
        else echo "Invalid parameter given!" && i=99 # invalid parameter and set i to error-code
        fi
                
    else
        echo $line >> $temp
    fi

    if [ $i -lt $changeLine ]; then i=$((i+1)) 
    fi

done < $file

# backup file of 'file'
fileBackup='backup_bash_aliases'

# make backup from file to fileBackup
cat $file > $fileBackup

if [ $i != 99 ]; then # i=99 is error code

    # write changed contents from temp to file
    cat $temp > $file

    case $i in
        100) echo "Changed current installation of " $programName " from standard to modified!";;
        101) echo "Changed current installation of " $programName " from modified to standard!";;
        102) echo "Changed current installation of " $programName " to modified!";;
        103) echo "Changed current installation of " $programName " to standard!";;
    esac

    echo "Please restart terminal so that changes can come into effect!"
fi

# delete the temporary file
rm $temp 






# backup of file, 06.10.2021 14:28


#std, mod
# apfelInstallation="mod"
# lhapdfInstallation="std"



###################
# apfelxx
###################
# case $apfelInstallation in
#     "std")
#     export LD_LIBRARY_PATH=/local0/a_eppi01/myapfel/lib:$LD_LIBRARY_PATH
#     export PATH=/local0/a_eppi01/myapfel/bin:$PATH
#     export CURRENT_APFEL=/local0/a_eppi01/myapfel ;;

#     "mod")
#     export LD_LIBRARY_PATH=/local0/a_eppi01/myApfelModified/lib:$LD_LIBRARY_PATH
#     export PATH=/local0/a_eppi01/myApfelModified/bin:$PATH
#     export CURRENT_APFEL=/local0/a_eppi01/myApfelModified;;
# esac



###################
# lhapdf
###################
# case $lhapdfInstallation in
#     "std")
#     export LD_LIBRARY_PATH=/local0/a_eppi01/myLHAPDF/lib:$LD_LIBRARY_PATH
#     export PATH=/local0/a_eppi01/myLHAPDF/bin:$PATH
#     export CURRENT_LHAPDF=/local0/a_eppi01/myLHAPDF ;;

#     "mod")
#     export LD_LIBRARY_PATH=/local0/a_eppi01/myLHAPDFModified/lib:$LD_LIBRARY_PATH
#     export PATH=/local0/a_eppi01/myLHAPDFModified/bin:$PATH
#     export CURRENT_LHAPDF=/local0/a_eppi01/myLHAPDFModified;;
# esac
# maybe needed:
# export PYTHONPATH=$PYTHONPATH:/foo/lhapdf/lib/python3.9/site-packages

