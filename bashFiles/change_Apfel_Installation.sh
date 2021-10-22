# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  cange_Apfel_Installation.sh                                                                              #
#  Author: Alexander Epping: a_eppi01@uni-muenster.de                                                       #
#  6 Oct 2021                                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #    
#  Convenient program to change the currently used installation of Apfel++ on my ThinkPad.                  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  - Furthermore it will save the current state of the .bash_aliases_apfel file in a backup file            #
#  - The .bash_aliases_apfel file should be constucted in a similar way to the one found in the Github.     #
#    It should be included into the .bash_aliases file.                                                     #
#  - The paths obiously should be changed accordingly.                                                      # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# line in the .bash_aliases_apfel file in wich the varibale for Apfel is saved
changeLine=2
programName="Apfel++"

# set "change" as the default value of the input
# possible inputs are mod and std
var=${1:-change}


# change to the folder with .bash_aliases_apfel
# cd /home/users/a_eppi01
cd /home/alexander
# bash file to be changed
file='.bash_aliases_apfel'

# temporary file to copy contents of 'file' to
temp='temp_bash_aliases_apfel.txt'

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
fileBackup='backup_bash_aliases_apfel'

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
